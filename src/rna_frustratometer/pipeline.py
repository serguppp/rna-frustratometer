import csv
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from statistics import fmean, pstdev

from gemmi import Structure

from .briq import BriqMutationResult, run_briq_mutations
from .mutation import MutationGenerator
from .structures import (
    load_structure_from_file_path,
    remove_water_and_metal_ions,
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
    def __init__(self, input_dir: Path, output_dir: Path) -> None:
        self.input_path = Path(input_dir).resolve()
        self.output_dir = Path(output_dir).resolve()

    def run(self) -> AnalysisResult:
        structure = self._load_structure()
        generator = self._prepare_variants(structure)

        briq_result = run_briq_mutations(
            generator,
            workspace=self.output_dir,
        )

        analysis = self._build_analysis(briq_result)
        self._save_summary(briq_result, analysis)

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

        remove_water_and_metal_ions(structure)
        return structure

    def _prepare_variants(self, structure: Structure) -> MutationGenerator:
        targets = {(chain.name, str(residue.seqid)) for chain in structure[0] for residue in chain}

        generator = MutationGenerator(structure)
        generator.prepare_residues_variants(targets)
        return generator

    def _read_contacts(self, pairs_file: Path) -> list[PairKey]:
        with pairs_file.open(newline="", encoding="utf-8") as file:
            return [
                (
                    row["chain_i"],
                    row["residue_i"],
                    row["chain_j"],
                    row["residue_j"],
                )
                for row in csv.DictReader(file)
                if row["contact_type"] not in ("", "unannotated")
            ]

    def _build_analysis(self, result: BriqMutationResult) -> AnalysisResult:
        contacts = self._read_contacts(result.pairs_file)
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
        briq_result: BriqMutationResult,
        analysis: AnalysisResult,
    ) -> None:
        summary_path = self.output_dir / "contact_variants.csv"

        if briq_result.variants_file.resolve() != summary_path.resolve():
            copyfile(briq_result.variants_file, summary_path)

        print(f"Original energy: {analysis.native_energy:.6f}")
        print(f"Contacts: {len(analysis.contacts)}")
        print(f"Saved: {summary_path}")
