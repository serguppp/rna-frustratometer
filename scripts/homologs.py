"""Compare the existing homolog panel from the editor."""

from datetime import datetime
from pathlib import Path

from rna_frustratometer.homologs import analyze_homologs

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
EXPERIMENT = RESULTS / "homolog_validation_20261114"

if __name__ == "__main__":
    analyze_homologs(
        RESULTS / "homologs_20261114/primary_panel.json",
        EXPERIMENT / "alignment_mapping.csv",
        EXPERIMENT / "scores",
        RESULTS / "homolog_comparison" / datetime.now().strftime("%Y%m%d-%H%M%S-%f"),
        plot_families=("RF00059",),
    )
