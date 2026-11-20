"""Prepare one RNA, score its variants with one potential, and compute F."""

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from typing import Literal

from gemmi import Structure

from .cgrnasp import CgrnaspMutationResult, run_cgrnasp_mutations
from .frustration import (
    STRATEGIES,
    UNIFORM_WEIGHTS,
    BaseFrequencies,
    Weighting,
    energy_statistics,
    read_contact_energies,
    save_frustration_scores,
)
from .simrna import PAIR_MINIMUM_SD, SimRNAConfig, SimRNAScanResult, run_simrna_frustration
from .structures import load_structure_from_file_path, retain_rna_only

PairKey = tuple[str, str, str, str]


@dataclass(frozen=True)
class ContactResult:
    chain_i: str
    residue_i: str
    chain_j: str
    residue_j: str
    native_energy: float
    variant_energies: tuple[float, ...]
    variant_weights: tuple[float, ...] = UNIFORM_WEIGHTS
    minimum_sd: float = 0.0

    @property
    def frustration(self) -> float | None:
        """Configured score; RNAfrustratometer supplies composition weights by default."""
        return energy_statistics(
            self.native_energy,
            self.variant_energies,
            self.variant_weights,
            minimum_sd=self.minimum_sd,
        ).frustration

    @property
    def baseline_frustration(self) -> float | None:
        return energy_statistics(
            self.native_energy, self.variant_energies, minimum_sd=self.minimum_sd
        ).frustration


@dataclass(frozen=True)
class AnalysisResult:
    native_energy: float
    contacts: tuple[ContactResult, ...]
    pair_contacts: tuple[ContactResult, ...] = ()
    background: BaseFrequencies | None = None
    weighting: Weighting = "composition"
    variants_file: Path | None = None
    scores_file: Path | None = None


class RNAfrustratometer:
    def __init__(
        self,
        input_dir: Path,
        output_dir: Path,
        *,
        executable: Path | None = None,
        data_dir: Path | None = None,
        simrna: SimRNAConfig | None = None,
        contacts_file: Path | None = None,
        sequence_separation: Literal["pdb", "ordinal"] = "ordinal",
        weighting: Weighting = "composition",
        background: BaseFrequencies | None = None,
    ) -> None:
        if weighting not in ("composition", "uniform"):
            raise ValueError("Weighting must be composition or uniform.")
        if simrna is None and (executable is None or data_dir is None):
            raise ValueError("cgRNASP requires executable and data_dir.")
        if simrna is not None and any(
            arg is not None for arg in (executable, data_dir, contacts_file)
        ):
            raise ValueError("Use SimRNAOptions to configure SimRNA contacts.")
        self.input_path = Path(input_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.contacts_file = Path(contacts_file).resolve() if contacts_file is not None else None
        self.executable = Path(executable).resolve() if executable is not None else None
        self.data_dir = Path(data_dir).resolve() if data_dir is not None else None
        self.simrna = simrna
        self.sequence_separation = sequence_separation
        self.weighting = weighting
        self.background = background

    def run(self) -> AnalysisResult:
        structure = self._load_structure()
        if self.output_dir.exists() and any(self.output_dir.iterdir()):
            raise FileExistsError("Use a new or empty analysis output directory.")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        reference = self.output_dir / "reference.pdb"
        structure.write_pdb(str(reference))
        (self.output_dir / "input.json").write_text(
            json.dumps(
                {
                    "source": str(self.input_path),
                    "sha256": hashlib.sha256(self.input_path.read_bytes()).hexdigest(),
                    "reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        background = self.background or BaseFrequencies.from_pdb(reference)
        result: CgrnaspMutationResult | SimRNAScanResult
        if self.simrna is None:
            assert self.executable is not None and self.data_dir is not None
            result = run_cgrnasp_mutations(
                reference_path=reference,
                output_dir=self.output_dir / "cgrnasp",
                contacts_file=self.contacts_file,
                executable=self.executable,
                data_dir=self.data_dir,
                sequence_separation=self.sequence_separation,
            )
        else:
            result = run_simrna_frustration(reference, self.output_dir / "simrna", self.simrna)
        variants_file = self.output_dir / "contact_variants.csv"
        copyfile(result.variants_file, variants_file)
        pair_minimum_sd = PAIR_MINIMUM_SD if self.simrna is not None else 0.0
        groups = read_contact_energies(variants_file)
        scores_file = save_frustration_scores(
            variants_file,
            self.output_dir,
            background,
            weighting=self.weighting,
            pair_minimum_sd=pair_minimum_sd,
            energies=groups,
        )
        contacts: dict[str, list[ContactResult]] = {"total": [], "pair": []}
        weights = background.weights(mode=self.weighting)
        for key, energies in groups.items():
            for strategy in STRATEGIES:
                value = energies[strategy]
                if strategy == "total" and not math.isclose(
                    value.native, result.total_energy, abs_tol=1e-8
                ):
                    raise ValueError(
                        "Native energy in the variant table differs from the scorer result."
                    )
                contacts[strategy].append(
                    ContactResult(
                        *key,
                        value.native,
                        value.variants,
                        weights,
                        pair_minimum_sd if strategy == "pair" else 0.0,
                    )
                )
        analysis = AnalysisResult(
            result.total_energy,
            tuple(contacts["total"]),
            tuple(contacts["pair"]),
            background,
            self.weighting,
            variants_file,
            scores_file,
        )
        print(f"Original energy: {analysis.native_energy:.6f}")
        print(f"Contacts: {len(analysis.contacts)} · weighting: {self.weighting}")
        print(f"Saved: {self.output_dir / 'frustration_scores.csv'}")
        return analysis

    def _load_structure(self) -> Structure:
        if not self.input_path.is_file():
            raise FileNotFoundError(self.input_path)
        if self.input_path.suffix.lower() != ".pdb":
            raise ValueError("Input must be a PDB file.")
        structure = load_structure_from_file_path(self.input_path)
        if not len(structure):
            raise ValueError("Input contains no structural models.")
        for index in range(len(structure) - 1, 0, -1):
            del structure[index]
        retain_rna_only(structure)
        if self.simrna is None and self.sequence_separation == "ordinal":
            self._validate_sequence_continuity(structure)
        return structure

    @staticmethod
    def _validate_sequence_continuity(structure: Structure) -> None:
        for chain in structure[0]:
            for index in range(1, len(chain)):
                left, right = chain[index - 1], chain[index]
                oxygen = next((a for a in left if a.name == "O3'"), None)
                phosphorus = next((a for a in right if a.name == "P"), None)
                label = f"{chain.name}:{left.seqid}-{right.seqid}"
                if oxygen is None or phosphorus is None:
                    raise ValueError(f"Cannot validate backbone continuity at {label}.")
                distance = oxygen.pos.dist(phosphorus.pos)
                if not 1.2 <= distance <= 2.2:
                    raise ValueError(f"Backbone break at {label}: O3'-P = {distance:.3f} A.")
