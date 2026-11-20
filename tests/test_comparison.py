import csv
from pathlib import Path

from rna_frustratometer.comparison import write_shared_scores


def test_shared_scores_match_reversed_pair_and_keep_undefined_f(tmp_path: Path) -> None:
    cgr = tmp_path / "cg.csv"
    sim = tmp_path / "sim.csv"
    header = (
        "chain_i,residue_i,chain_j,residue_j,"
        "F_pair_composition,F_pair_uniform,F_total_composition,F_total_uniform\n"
    )
    cgr.write_text(header + "A,1,A,2,1,2,3,4\nA,2,A,3,5,6,7,8\n", encoding="utf-8")
    sim.write_text(header + "A,2,A,1,,2,3,4\nA,1,A,3,5,6,7,8\n", encoding="utf-8")
    output = tmp_path / "shared.csv"
    assert write_shared_scores(cgr, sim, output) == 1
    with output.open(newline="", encoding="utf-8") as stream:
        row = next(csv.DictReader(stream))
    assert row["cgrnasp_F_pair_composition"] == "1"
    assert row["simrna_F_pair_composition"] == ""
    assert row["simrna_F_total_composition"] == "3"
