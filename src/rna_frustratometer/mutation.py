from copy import deepcopy
from io import StringIO
from itertools import product
from pathlib import Path

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue as BioResidue
from gemmi import Residue, Structure, read_pdb_string
from rna_tools.tools.mini_moderna3.moderna import (  # type: ignore[import-untyped]
    exchange_single_base,
    load_model,
)

BASES = ("A", "C", "G", "U")

ResidueKey = tuple[str, str]  # (chain_id, residue_id)
VariantKey = tuple[str, str, str]  # (chain_id, residue_id, base)


class MutationGenerator:
    def __init__(self, structure: Structure) -> None:
        if len(structure) != 1:
            raise ValueError("Select exactly one structural model.")

        self.reference = structure.clone()
        self.variants: dict[VariantKey, Residue] = {}

        parser = PDBParser(QUIET=True)  # type: ignore[no-untyped-call]
        self.bio_structure = parser.get_structure(  # type: ignore[no-untyped-call]
            "rna",
            StringIO(self.reference.make_pdb_string()),
        )

    def prepare_residues_variants(self, targets: set[ResidueKey]) -> None:
        """Prepare and cache A/C/G/U variants for the requested residues."""
        chain_ids = sorted({chain_id for chain_id, _ in targets})

        for chain_id in chain_ids:
            residue_ids = sorted(residue_id for chain, residue_id in targets if chain == chain_id)

            self._prepare_chain_variants(chain_id, residue_ids)

    def _prepare_chain_variants(
        self,
        chain_id: str,
        residue_ids: list[str],
    ) -> None:
        bio_structure = self.bio_structure
        if bio_structure is None:
            raise RuntimeError("BioPython structure not loaded.")

        model = load_model(
            bio_structure[0][chain_id],
            chain_id,
            data_type="chain",
        )

        if model is None:
            raise RuntimeError(f"Could not load chain {chain_id}.")

        reference_chain = self.reference[0][chain_id]

        for residue_id in residue_ids:
            matches = [residue for residue in reference_chain if str(residue.seqid) == residue_id]

            if len(matches) != 1:
                raise ValueError(f"Expected one reference residue: {chain_id}:{residue_id}")

            reference_residue = matches[0]

            for base in BASES:
                key = (chain_id, residue_id, base)

                if key in self.variants:
                    continue

                if base == reference_residue.name:
                    self.variants[key] = reference_residue.clone()
                    continue

                source = model[residue_id]

                if not isinstance(source, BioResidue):
                    raise TypeError(
                        f"Expected a single residue at {chain_id}:{residue_id}, "
                        f"got {type(source).__name__}."
                    )

                variant = self._build_variant(source, base)

                if str(variant.seqid) != residue_id:
                    raise ValueError(f"Residue identifier changed: {key}")

                self.variants[key] = variant

    def _build_variant(self, source: BioResidue, base: str) -> Residue:
        """Rebuild a copy of the residue and return it as a Gemmi residue."""
        residue = deepcopy(source)
        exchange_single_base(residue, base)

        return self._to_gemmi(residue)

    def _to_gemmi(self, residue: BioResidue) -> Residue:
        """Convert a BioPython residue to Gemmi through an in-memory PDB."""
        buffer = StringIO()

        writer = PDBIO()  # type: ignore[no-untyped-call]
        writer.set_structure(residue)  # type: ignore[no-untyped-call]
        writer.save(buffer)  # type: ignore[no-untyped-call]

        structure = read_pdb_string(buffer.getvalue())
        return structure[0][0][0].clone()

    def write_pair_mutants(
        self,
        chain_i: str,
        residue_i: str,
        chain_j: str,
        residue_j: str,
        output_dir: Path,
    ) -> list[tuple[str, str, Path]]:
        """Write all 16 combinations using the prepared residue cache."""
        if (chain_i, residue_i) == (chain_j, residue_j):
            raise ValueError("Contact must contain two distinct residues.")

        output_dir.mkdir(parents=True, exist_ok=True)
        results: list[tuple[str, str, Path]] = []

        for base_i, base_j in product(BASES, repeat=2):
            mutant = self.reference.clone()

            self._replace_residue(mutant, chain_i, residue_i, base_i)
            self._replace_residue(mutant, chain_j, residue_j, base_j)

            output_path = output_dir / f"{base_i}{base_j}.pdb"
            mutant.write_pdb(str(output_path))

            results.append((base_i, base_j, output_path))

        return results

    def _replace_residue(
        self,
        structure: Structure,
        chain_id: str,
        residue_id: str,
        base: str,
    ) -> None:
        """Replace one residue with a copy of its cached variant."""
        chain = structure[0][chain_id]

        indices = [index for index, residue in enumerate(chain) if str(residue.seqid) == residue_id]

        if len(indices) != 1:
            raise ValueError(f"Expected one residue: {chain_id}:{residue_id}")

        index = indices[0]
        replacement = self.variants[(chain_id, residue_id, base)].clone()

        del chain[index]
        chain.add_residue(replacement, pos=index)
