from pathlib import Path

import gemmi
import pytest

from rna_frustratometer.simrna_mutation import atom_position, glycosidic_frame, replace_base


def test_cross_class_mutation_preserves_native_anchor_and_glycosidic_frame() -> None:
    reference = gemmi.read_structure(str(Path(__file__).parent / "data/2KOC.pdb"))
    residues = [residue for chain in reference[0] for residue in chain]
    guanine = next(r for r in residues if r.name == "G")
    cytosine = next(r for r in residues if r.name == "C")
    variant = replace_base(guanine, cytosine)

    assert variant.name == "C"
    assert atom_position(variant, "N1").dist(atom_position(guanine, "N9")) == 0
    for name in ("C4'", "C1'", "O4'"):
        assert atom_position(variant, name).dist(atom_position(guanine, name)) == 0
    for native_axis, mutant_axis in zip(
        glycosidic_frame(guanine)[1], glycosidic_frame(variant)[1], strict=True
    ):
        assert native_axis.dot(mutant_axis) == pytest.approx(1.0)
    assert atom_position(variant, "C2").dist(atom_position(variant, "N1")) == pytest.approx(
        atom_position(cytosine, "C2").dist(atom_position(cytosine, "N1"))
    )
