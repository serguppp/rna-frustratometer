from pathlib import Path

from rna_frustratometer.pipeline import RNAfrustratometer
from rna_frustratometer.visualization import create_frustration_plots

output = Path("results/1ANR_cgrnasp")
tool = Path.home() / "Desktop/tools/cgRNASP"
analysis = RNAfrustratometer(
    input_dir=Path("tests/data/1ANR.pdb"),
    output_dir=output,
    executable=tool / "build/cgRNASP_frustration.exe",
    data_dir=tool / "cgRNASP/data",
    sequence_separation="ordinal",
)
result = analysis.run()
figures = create_frustration_plots(
    result,
    output / "reference.pdb",
    components_path=output / "cgrnasp/energy_components.csv",
)
for name, figure in figures.items():
    figure.write_html(str(output / f"{name}.html"))
