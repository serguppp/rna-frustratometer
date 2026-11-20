"""Edit the paths below and run from the editor (or use uv run scripts/analysis.py)."""

import os
from datetime import datetime
from pathlib import Path

from rna_frustratometer import RNAfrustratometer, SimRNAConfig, SimRNAOptions
from rna_frustratometer.comparison import POTENTIALS, Potential, read_scores, write_shared_scores
from rna_frustratometer.frustration import STRATEGIES
from rna_frustratometer.visualization import contact_maps

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data/sim"
TOOLS = Path.home() / "Desktop/tools"
CGRNASP = TOOLS / "cgRNASP"
SIMRNA = TOOLS / "SimRNA-3.20/installation_dynamic/SimRNA_64bitIntel_Linux"
TEMPLATES = DATA / "2KOC.pdb"
INPUTS = (DATA / "2KOC.pdb",)


def run(pdb: Path, output: Path, *, potentials: tuple[Potential, ...] = POTENTIALS) -> None:
    if not potentials or len(set(potentials)) != len(potentials):
        raise ValueError("Choose cgRNASP, SimRNA or both.")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Use a new or empty output directory: {output}")
    scores = {}
    for potential in potentials:
        destination = output / potential
        if potential == "cgrnasp":
            analysis = RNAfrustratometer(
                pdb,
                destination,
                executable=CGRNASP
                / (
                    "build/cgRNASP_frustration.exe"
                    if os.name == "nt"
                    else "build/cgRNASP_frustration"
                ),
                data_dir=CGRNASP / "cgRNASP/data",
            )
        else:
            analysis = RNAfrustratometer(
                pdb,
                destination,
                simrna=SimRNAConfig(SIMRNA, TEMPLATES, SimRNAOptions()),
            )
        result = analysis.run()
        assert result.scores_file is not None
        scores[potential] = result.scores_file

    if set(scores) == set(POTENTIALS):
        n = write_shared_scores(scores["cgrnasp"], scores["simrna"], output / "shared_scores.csv")
        print(f"Shared mutation targets: {n}")
    tables = {name: read_scores(path) for name, path in scores.items()}
    for strategy in STRATEGIES:
        column = f"F_{strategy}"
        values = {
            name: {key: float(row[column]) for key, row in table.items() if row[column]}
            for name, table in tables.items()
        }
        reference = output / potentials[0] / "reference.pdb"
        contact_maps(reference, values, output / f"{column}.png")


if __name__ == "__main__":
    batch = ROOT / "results/scans" / datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    for structure in INPUTS:
        run(structure, batch / structure.stem)
