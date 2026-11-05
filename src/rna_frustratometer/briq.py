import csv
import json
import math
import subprocess
from dataclasses import dataclass
from itertools import product
from pathlib import Path

from .mutation import BASES, MutationGenerator

PairKey = tuple[str, str, str, str]

VARIANT_COLUMNS = [
    "chain_i",
    "residue_i",
    "chain_j",
    "residue_j",
    "contact_type",
    "base_i",
    "base_j",
    "native_E_total",
    "E_total",
    "delta_E_total",
]


@dataclass(frozen=True)
class BriqMutationResult:
    total_energy: float
    pairs_file: Path
    variants_file: Path
    run_file: Path


def run_briq_mutations(
    cache: MutationGenerator,
    workspace: Path,
    *,
    image: str = "briq:memory",
) -> BriqMutationResult:
    """Score the reference and all contact variants in one BRiQ process."""
    workspace = workspace.resolve()

    _prepare_inputs(cache, workspace)
    _run_docker(workspace, image)

    return _read_result(workspace / "briq")


def _prepare_inputs(
    cache: MutationGenerator,
    workspace: Path,
) -> None:
    """Write the reference structure and cached variant coordinates."""
    workspace.mkdir(parents=True, exist_ok=True)

    reference_path = workspace / "reference.pdb"
    geometry_path = workspace / "geometry.tsv"
    output_dir = workspace / "briq"

    if any(path.exists() for path in (reference_path, geometry_path, output_dir)):
        raise FileExistsError("Use a new analysis output directory.")

    cache.reference.write_pdb(str(reference_path))
    _write_geometry(cache, geometry_path)


def _write_geometry(
    cache: MutationGenerator,
    path: Path,
) -> None:
    """Export heavy-atom coordinates for every cached residue variant."""
    with path.open("x", encoding="utf-8", newline="") as file:
        writer = csv.writer(
            file,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(["chain", "residue", "base", "atom", "x", "y", "z"])

        for key, residue in sorted(cache.variants.items()):
            chain_id, residue_id, base = key

            if residue.name != base:
                raise ValueError(f"Unexpected cached residue name: {key}")

            for atom in residue:
                if atom.element.name in ("H", "D"):
                    continue

                if atom.altloc != "\x00":
                    raise ValueError(f"Select one atom conformation first: {key}, {atom.name}")

                coordinates = (atom.pos.x, atom.pos.y, atom.pos.z)

                if not all(math.isfinite(value) for value in coordinates):
                    raise ValueError(f"Non-finite coordinate: {key}, {atom.name}")

                writer.writerow(
                    [
                        chain_id,
                        residue_id,
                        base,
                        atom.name,
                        *(f"{value:.3f}" for value in coordinates),
                    ]
                )


def _run_docker(workspace: Path, image: str) -> None:
    """Run BRiQ and write its output to memory.log."""
    command = [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={workspace},target=/work",
        image,
        "BRiQ_Energy",
        "--mutate-total",
        "/work/reference.pdb",
        "/work/geometry.tsv",
        "/work/briq",
    ]

    log_path = workspace / "memory.log"

    try:
        with log_path.open("w", encoding="utf-8") as log:
            subprocess.run(
                command,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"BRiQ failed. See the log: {log_path}") from error


def _read_result(output_dir: Path) -> BriqMutationResult:
    """Read run metadata and validate the generated mutation results."""
    run_file = output_dir / "mutation_run.json"
    pairs_file = output_dir / "pairs.csv"
    variants_file = output_dir / "contact_variants.csv"

    for path in (run_file, pairs_file, variants_file):
        if not path.is_file():
            raise FileNotFoundError(f"Missing BRiQ output: {path}")

    metadata = json.loads(run_file.read_text(encoding="utf-8"))

    native_energy = float(metadata["native_energy"])

    if not math.isfinite(native_energy):
        raise ValueError("Invalid native energy.")

    _validate_variants(
        variants_file,
        native_energy=native_energy,
        expected_variants=metadata["variants"],
        expected_contacts=metadata["contacts"],
    )

    return BriqMutationResult(
        total_energy=native_energy,
        pairs_file=pairs_file,
        variants_file=variants_file,
        run_file=run_file,
    )


def _validate_variants(
    path: Path,
    *,
    native_energy: float,
    expected_variants: int,
    expected_contacts: int,
) -> None:
    """Check energies, variant combinations, and agreement with metadata."""
    groups: dict[PairKey, list[tuple[str, str]]] = {}
    row_count = 0

    with path.open(newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)

        if reader.fieldnames != VARIANT_COLUMNS:
            raise ValueError("Invalid mutation CSV header.")

        for row in reader:
            key = (
                row["chain_i"],
                row["residue_i"],
                row["chain_j"],
                row["residue_j"],
            )

            bases = (row["base_i"], row["base_j"])
            groups.setdefault(key, []).append(bases)

            _validate_energy(row, native_energy)
            row_count += 1

    expected_combinations = list(product(BASES, repeat=2))

    for key, combinations in groups.items():
        if combinations != expected_combinations:
            raise ValueError(f"Expected 16 ordered variants for contact {key}.")

    if row_count != expected_variants:
        raise ValueError("Variant count does not match run metadata.")

    if len(groups) != expected_contacts:
        raise ValueError("Contact count does not match run metadata.")


def _validate_energy(
    row: dict[str, str],
    native_energy: float,
) -> None:
    """Check that a CSV row contains finite and consistent energies."""
    energy = float(row["E_total"])
    reference_energy = float(row["native_E_total"])
    delta = float(row["delta_E_total"])

    if not all(math.isfinite(value) for value in (energy, reference_energy, delta)):
        raise ValueError("Invalid mutant energy.")

    if abs(reference_energy - native_energy) > 1e-9:
        raise ValueError("Inconsistent reference energy.")

    if abs(delta - (energy - native_energy)) > 1e-9:
        raise ValueError("Inconsistent energy difference.")
