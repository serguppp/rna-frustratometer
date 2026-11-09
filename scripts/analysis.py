from pathlib import Path

from rna_frustratometer.pipeline import RNAfrustratometer
from rna_frustratometer.visualization import create_frustration_plots

data = ["2K96"]
for pdb_id in data:
    analysis = RNAfrustratometer(
        input_dir=Path(f"tests/data/{pdb_id}.pdb"),
        output_dir=Path(f"results/{pdb_id}_cache_01"),
    )

    result = analysis.run()

    figures = create_frustration_plots(
        result,
        analysis.output_dir / "reference.pdb",
        components_path=analysis.output_dir / "briq" / "energy_components.csv",
        radius=8.0,
    )

    for name, figure in figures.items():
        figure.write_html(str(analysis.output_dir / f"{name}.html"))
