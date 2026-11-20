# Adapted from serguppp/simrna-adapter (31be7d6); see LICENSE-SIMRNA-ADAPTER.

"""Evaluate a complete SimRNA mutation scan in one Linux container."""

import csv
import json
import math
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import gemmi

from .frustration import BASES, COMBINATIONS, PAIR_COLUMNS, PairKey
from .simrna_engine import SimRNAResult, _engine_provenance, _native_score, _prepare_input, _sha256
from .simrna_mutation import (
    MutationState,
    mutation_state,
    read_templates,
    replace_base,
    select_contacts,
    write_variant,
)

PAIR_MINIMUM_SD = 2e-6


@dataclass(frozen=True)
class SimRNAOptions:
    cutoff: float = 8.0
    workers: int = 4
    timeout: float = 120.0
    preserve_same_class_ring: bool = False
    contacts: tuple[PairKey, ...] | None = None

    def __post_init__(self) -> None:
        if self.workers < 1 or not all(
            math.isfinite(v) and v > 0 for v in (self.cutoff, self.timeout)
        ):
            raise ValueError("Workers, cutoff and timeout must be positive.")


@dataclass(frozen=True)
class SimRNAConfig:
    engine_dir: Path
    templates_pdb: Path
    options: SimRNAOptions = SimRNAOptions()
    docker_image: str | None = "rna-frustratometer-simrna:local"


@dataclass(frozen=True)
class SimRNAScanResult:
    total_energy: float
    variants_file: Path
    run_file: Path


def run_simrna_frustration(reference: Path, output: Path, config: SimRNAConfig) -> SimRNAScanResult:
    reference = reference.resolve(strict=True)
    templates = config.templates_pdb.resolve(strict=True)
    engine = config.engine_dir.resolve(strict=True)
    output = output.resolve()
    if any(p.suffix.lower() != ".pdb" for p in (reference, templates)):
        raise ValueError("SimRNA reference and templates must be PDB files.")
    for name in ("SimRNA", "config_n0.dat"):
        if not (engine / name).is_file():
            raise FileNotFoundError(engine / name)
    if not (engine / "data").is_dir():
        raise NotADirectoryError(engine / "data")
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError(f"Use a new or empty SimRNA output directory: {output}")

    if config.docker_image is None:
        _scan(reference, templates, output, engine, config.options)
    else:
        _run_container(reference, templates, output, engine, config.docker_image, config.options)
    run_file = output / "run.json"
    metadata = json.loads(run_file.read_text(encoding="utf-8"))
    variants_file = output / "contact_variants.csv"
    if not variants_file.is_file():
        raise FileNotFoundError(variants_file)
    return SimRNAScanResult(float(metadata["E_total"]), variants_file, run_file)


def _run_container(
    reference: Path, templates: Path, output: Path, engine: Path, image: str, options: SimRNAOptions
) -> None:
    with TemporaryDirectory(prefix="simrna-options-") as temporary:
        manifest = Path(temporary) / "options.json"
        manifest.write_text(json.dumps(asdict(options)), encoding="utf-8")
        command = ["docker", "run", "--rm", "--network", "none", "--entrypoint", "python"]
        for source, target in (
            (reference, "/input.pdb"),
            (templates, "/templates.pdb"),
            (engine, "/simrna"),
            (manifest, "/options.json"),
            (Path(__file__).resolve().parent, "/code/rna_frustratometer"),
        ):
            command.extend(["--mount", f"type=bind,source={source},target={target},readonly"])
        command.extend(
            [
                "--env",
                "PYTHONPATH=/code",
                "--mount",
                f"type=bind,source={output},target=/output",
                image,
                "-c",
                "from rna_frustratometer.simrna import _worker; _worker()",
            ]
        )
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        log = output / "docker.log"
        log.write_text(process.stdout + process.stderr, encoding="utf-8")
        if process.returncode:
            raise RuntimeError(f"SimRNA scan exited with {process.returncode}. See {log}.")


def _worker() -> None:
    payload = json.loads(Path("/options.json").read_text(encoding="utf-8"))
    if payload["contacts"] is not None:
        payload["contacts"] = tuple(tuple(key) for key in payload["contacts"])
    _scan(
        Path("/input.pdb"),
        Path("/templates.pdb"),
        Path("/output"),
        Path("/simrna"),
        SimRNAOptions(**payload),
    )


def isolate_pair(structure: gemmi.Structure, pair: tuple[int, int]) -> gemmi.Structure:
    """Move other residues apart while retaining chain membership and adjacency."""
    masked = structure.clone()
    residues = [residue for chain in masked[0] for residue in chain]
    offset = max(1000.0, max(atom.pos.x for residue in residues for atom in residue) + 100.0)
    for index, residue in enumerate(residues):
        if index in pair:
            continue
        anchor = residue[0].pos
        shift = gemmi.Vec3(
            offset + 50 * (index % 40) - anchor.x,
            1000 + 50 * ((index // 40) % 40) - anchor.y,
            1000 + 50 * (index // 1600) - anchor.z,
        )
        for atom in residue:
            position = gemmi.Position(
                atom.pos.x + shift.x, atom.pos.y + shift.y, atom.pos.z + shift.z
            )
            if any(not -999 < value < 9999 for value in (position.x, position.y, position.z)):
                raise ValueError("Pair isolation exceeds PDB coordinate limits.")
            atom.pos = position
    return masked


def _pair_energy(result: SimRNAResult) -> float:
    """Discard artificial backbone strain from the isolated-pair score."""
    energy = result.total_energy - result.components["E_local_geometry"]
    subtotal = result.components["E_base_base"] + result.components["E_base_backbone"]
    if not math.isclose(energy, subtotal, rel_tol=0.0, abs_tol=0.002):
        raise ValueError("Isolated pair energy does not match the nonlocal components.")
    return energy


def _scan(
    pdb: Path, template_path: Path, output: Path, engine: Path, options: SimRNAOptions
) -> None:
    started = perf_counter()
    with TemporaryDirectory(prefix="simrna-scan-") as folder:
        work = Path(folder)
        # Keep the engine and temporary structures on the Linux filesystem for the whole scan.
        local_engine = work / "engine"
        local_engine.mkdir()
        for name in ("SimRNA", "config_n0.dat"):
            shutil.copy2(engine / name, local_engine / name)
        shutil.copytree(engine / "data", local_engine / "data")
        provenance = _engine_provenance(local_engine)

        reference_pdb = work / "reference.pdb"
        _prepare_input(pdb, reference_pdb)
        structure = gemmi.read_structure(str(reference_pdb))
        full_atom = gemmi.read_structure(str(pdb))
        templates, template_ids = read_templates(template_path)
        locations = [
            (ci, ri) for ci, chain in enumerate(structure[0]) for ri, _ in enumerate(chain)
        ]
        residues = [structure[0][ci][ri] for ci, ri in locations]
        identifiers = [
            (structure[0][ci].name, str(residue.seqid))
            for (ci, _), residue in zip(locations, residues, strict=True)
        ]
        lookup = {key: index for index, key in enumerate(identifiers)}
        native = [residue.name for residue in residues]
        contacts = (
            options.contacts
            if options.contacts is not None
            else select_contacts(structure, options.cutoff)
        )
        pairs = [(lookup[ci, ri], lookup[cj, rj]) for ci, ri, cj, rj in contacts]
        if not pairs or any(i == j for i, j in pairs):
            raise ValueError("Select at least one pair of distinct residues.")
        if len({(min(i, j), max(i, j)) for i, j in pairs}) != len(pairs):
            raise ValueError("Duplicate SimRNA contact.")

        replacements = {}
        for index in sorted({index for pair in pairs for index in pair}):
            ci, ri = locations[index]
            for base in BASES:
                if base != native[index]:
                    replacements[index, base] = replace_base(
                        full_atom[0][ci][ri],
                        templates[base],
                        preserve_same_class_ring=options.preserve_same_class_ring,
                    )

        states = sorted(
            {mutation_state(i, j, a, b, native) for i, j in pairs for a, b in COMBINATIONS}
        )
        preparation_seconds = perf_counter() - started
        reference_dir = output / "reference"
        reference_dir.mkdir()
        reference = _native_score(
            reference_pdb, reference_dir, local_engine, 1, options.timeout, provenance
        )

        def evaluate(state: MutationState, pair: tuple[int, int] | None = None) -> SimRNAResult:
            with TemporaryDirectory(dir=work, prefix="variant-") as directory:
                variant_dir = Path(directory)
                variant = variant_dir / "variant.pdb"
                score_dir = variant_dir / "score"
                score_dir.mkdir()
                geometry = structure if pair is None else isolate_pair(structure, pair)
                write_variant(geometry, locations, replacements, state, variant)
                try:
                    return _native_score(
                        variant, score_dir, local_engine, 1, options.timeout, provenance
                    )
                except Exception as error:
                    label = "_".join(f"{index}{base}" for index, base in state) or "native"
                    prefix = "total" if pair is None else f"pair_{pair[0]}_{pair[1]}"
                    failure = output / f"failed_{prefix}_{label}"
                    shutil.copytree(variant_dir, failure)
                    raise RuntimeError(f"SimRNA scoring failed. See {failure}.") from error

        def score_total(state: MutationState) -> tuple[float, dict[str, float]]:
            result = evaluate(state)
            return result.total_energy, result.components

        def score_pair(task: tuple[int, int, str, str]) -> float:
            i, j, a, b = task
            return _pair_energy(evaluate(mutation_state(i, j, a, b, native), (i, j)))

        total_states = [state for state in states if state]
        tasks = [(i, j, a, b) for i, j in pairs for a, b in COMBINATIONS]
        with ThreadPoolExecutor(min(options.workers, os.cpu_count() or 1)) as pool:
            totals = dict(zip(total_states, pool.map(score_total, total_states), strict=True))
            totals[()] = (reference.total_energy, reference.components)
            pair_scores = dict(zip(tasks, pool.map(score_pair, tasks), strict=True))

    with (output / "contact_variants.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                *PAIR_COLUMNS,
                "base_i",
                "base_j",
                "E_total",
                "delta_E_total",
                *reference.components,
                "native_E_total",
                "native_E_pair",
                "E_pair",
                "delta_E_pair",
            ]
        )
        for i, j in pairs:
            native_pair = pair_scores[i, j, native[i], native[j]]
            for a, b in COMBINATIONS:
                energy, components = totals[mutation_state(i, j, a, b, native)]
                pair_energy = pair_scores[i, j, a, b]
                writer.writerow(
                    [
                        *identifiers[i],
                        *identifiers[j],
                        a,
                        b,
                        energy,
                        energy - reference.total_energy,
                        *(components[name] for name in reference.components),
                        reference.total_energy,
                        native_pair,
                        pair_energy,
                        pair_energy - native_pair,
                    ]
                )

    (output / "run.json").write_text(
        json.dumps(
            {
                **provenance,
                "engine": "SimRNA 3.20",
                "iterations": 0,
                "seed": 1,
                "mutation_protocol_version": 4 if options.preserve_same_class_ring else 3,
                "preserve_same_class_ring": options.preserve_same_class_ring,
                "contact_selection": "base_midpoint_distance"
                if options.contacts is None
                else "explicit",
                "contact_cutoff_angstrom": options.cutoff if options.contacts is None else None,
                "reference_sha256": _sha256(pdb),
                "template_sha256": _sha256(template_path),
                "mutation_source_sha256": _sha256(Path(__file__).with_name("simrna_mutation.py")),
                "template_residues": template_ids,
                "E_total": reference.total_energy,
                "contact_count": len(pairs),
                "variant_count": len(tasks),
                "engine_calls": len(states) + len(tasks),
                "total_engine_calls": len(states),
                "pair_engine_calls": len(tasks),
                "energy_strategies": ["total", "pair"],
                "pair_F_minimum_sd": PAIR_MINIMUM_SD,
                "workers": min(options.workers, os.cpu_count() or 1),
                "preparation_seconds": preparation_seconds,
                "worker_seconds": perf_counter() - started,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
