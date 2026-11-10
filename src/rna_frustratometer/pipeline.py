import csv
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from statistics import fmean, pstdev
from typing import Literal

from gemmi import Structure

from .cgrnasp import CgrnaspMutationResult, run_cgrnasp_mutations
from .structures import (
    load_structure_from_file_path,
    retain_rna_only,
)

PairKey = tuple[str, str, str, str]  # (chain_i, residue_i, chain_j, residue_j)
Contact = tuple[
    PairKey, dict[str, str]
]  # PairKey; pairs.csv row e.g. {"contact_type": "basepair", "base_i": "A", "base_j": "U"}


@dataclass(frozen=True)
class ContactResult:
    chain_i: str
    residue_i: str
    chain_j: str
    residue_j: str
    native_energy: float
    variant_energies: tuple[float, ...]

    @property
    def frustration(self) -> float | None:
        if len(self.variant_energies) != 16:
            raise ValueError("Expected exactly 16 variant energies.")

        sigma = pstdev(self.variant_energies)
        if sigma == 0:
            return None

        return (fmean(self.variant_energies) - self.native_energy) / sigma


@dataclass(frozen=True)
class AnalysisResult:
    native_energy: float
    contacts: tuple[ContactResult, ...]


class RNAfrustratometer:
    def __init__(
        self,
        input_dir: Path,
        output_dir: Path,
        *,
        contacts_file: Path | None = None,
        executable: Path,
        data_dir: Path,
        sequence_separation: Literal["pdb", "ordinal"] = "ordinal",
    ) -> None:
        self.input_path = Path(input_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.contacts_file = Path(contacts_file).resolve() if contacts_file is not None else None
        self.executable = Path(executable).resolve()
        self.data_dir = Path(data_dir).resolve()
        self.sequence_separation = sequence_separation

    def run(self) -> AnalysisResult:
        structure = self._load_structure()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        reference_path = self.output_dir / "reference.pdb"
        if reference_path.exists():
            raise FileExistsError("Use a new analysis output directory.")
        structure.write_pdb(str(reference_path))

        cgrnasp_result = run_cgrnasp_mutations(
            reference_path=reference_path,
            contacts_file=self.contacts_file,
            output_dir=self.output_dir / "cgrnasp",
            executable=self.executable,
            data_dir=self.data_dir,
            sequence_separation=self.sequence_separation,  # type: ignore
        )

        analysis = self._build_analysis(cgrnasp_result)
        self._save_summary(cgrnasp_result, analysis)

        return analysis

    def _load_structure(self) -> Structure:
        if not self.input_path.is_file():
            raise FileNotFoundError(self.input_path)

        if self.input_path.suffix.lower() != ".pdb":
            raise ValueError("Input must be a PDB file.")

        structure = load_structure_from_file_path(self.input_path)

        if len(structure) == 0:
            raise ValueError("Input contains no structural models.")

        for index in range(len(structure) - 1, 0, -1):
            del structure[index]

        # remove_water_and_metal_ions(structure)
        retain_rna_only(structure)
        if self.sequence_separation == "ordinal":
            self._validate_sequence_continuity(structure)

        return structure

    @staticmethod
    def _validate_sequence_continuity(structure: Structure) -> None:
        """Reject ordinal indexing across unresolved or broken RNA backbones."""
        for chain in structure[0]:
            for index in range(1, len(chain)):
                left, right = chain[index - 1], chain[index]
                oxygen = next((atom for atom in left if atom.name == "O3'"), None)
                phosphorus = next((atom for atom in right if atom.name == "P"), None)
                label = f"{chain.name}:{left.seqid}-{right.seqid}"
                if oxygen is None or phosphorus is None:
                    raise ValueError(f"Cannot validate backbone continuity at {label}.")
                distance = oxygen.pos.dist(phosphorus.pos)
                if not 1.2 <= distance <= 2.2:
                    raise ValueError(
                        f"Backbone break at {label}: O3'-P = {distance:.3f} A. "
                        "Verify the construct sequence before choosing a separation protocol."
                    )

    def _read_contacts(self, variants_file: Path) -> list[PairKey]:
        with variants_file.open(newline="", encoding="utf-8") as file:
            return list(
                dict.fromkeys(
                    (row["chain_i"], row["residue_i"], row["chain_j"], row["residue_j"])
                    for row in csv.DictReader(file)
                )
            )

    def _build_analysis(self, result: CgrnaspMutationResult) -> AnalysisResult:
        contacts = self._read_contacts(result.variants_file)
        energies: dict[PairKey, list[float]] = {key: [] for key in contacts}

        with result.variants_file.open(newline="", encoding="utf-8") as file:
            for row in csv.DictReader(file):
                key = (
                    row["chain_i"],
                    row["residue_i"],
                    row["chain_j"],
                    row["residue_j"],
                )

                if key not in energies:
                    raise ValueError(f"Variant refers to an unknown contact: {key}")

                energies[key].append(float(row["E_total"]))

        contact_results: list[ContactResult] = []

        for key in contacts:
            variant_energies = energies[key]

            if len(variant_energies) != 16:
                raise ValueError(
                    f"Expected 16 energies for contact {key}, got {len(variant_energies)}."
                )

            chain_i, residue_i, chain_j, residue_j = key

            contact_results.append(
                ContactResult(
                    chain_i=chain_i,
                    residue_i=residue_i,
                    chain_j=chain_j,
                    residue_j=residue_j,
                    native_energy=result.total_energy,
                    variant_energies=tuple(variant_energies),
                )
            )

        return AnalysisResult(
            native_energy=result.total_energy,
            contacts=tuple(contact_results),
        )

    def _save_summary(
        self,
        cgrnasp_result: CgrnaspMutationResult,
        analysis: AnalysisResult,
    ) -> None:
        summary_path = self.output_dir / "contact_variants.csv"

        if cgrnasp_result.variants_file.resolve() != summary_path.resolve():
            copyfile(cgrnasp_result.variants_file, summary_path)

        print(f"Original energy: {analysis.native_energy:.6f}")
        print(f"Contacts: {len(analysis.contacts)}")
        print(f"Saved: {summary_path}")
