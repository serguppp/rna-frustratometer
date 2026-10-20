from pathlib import Path

import gemmi
import pytest


def test_gemmi_reads_rna_residue_from_pdb(tmp_path: Path) -> None:
    pdb = "ATOM      1  P     A A  10B      1.000   2.000   3.000  1.00 20.00           P  \nEND\n"
    path = tmp_path / "rna.pdb"
    path.write_text(pdb, encoding="ascii")

    structure = gemmi.read_structure(str(path))

    assert len(structure) == 1

    chain = structure[0][0]
    assert chain.name == "A"

    residue = chain[0]
    assert residue.name == "A"
    assert residue.seqid.num == 10
    assert residue.seqid.icode == "B"

    atom = residue[0]
    assert atom.name == "P"
    assert (atom.pos.x, atom.pos.y, atom.pos.z) == pytest.approx((1.0, 2.0, 3.0))
