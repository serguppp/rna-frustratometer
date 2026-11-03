from pathlib import Path
from tempfile import TemporaryDirectory

from gemmi import Structure, cif, make_structure_from_block, read_pdb_string, read_structure
from rna_tools.tools.mini_moderna3.moderna import exchange_single_base, load_model, write_model  # type: ignore[import-untyped]

from .formats import SupportedFormats


def load_structure_from_file_path(path: str | Path) -> Structure:
    return read_structure(str(path))


def load_structure_from_file_content(file_content: str, file_format: SupportedFormats) -> Structure:
    if SupportedFormats.normalize_format(file_format) == SupportedFormats.MMCIF:
        cif_document = cif.read_string(file_content)
        structure = make_structure_from_block(cif_document.sole_block())
    else:
        structure = read_pdb_string(file_content)

    return structure


def remove_water_and_metal_ions(structure: Structure) -> None:
    structure.remove_waters()

    for model in structure:
        for chain in model:
            for index in range(len(chain) - 1, -1, -1):
                residue = chain[index]

                if len(residue) == 1 and residue[0].element.is_metal:
                    del chain[index]


def mutate_rna(path: str | Path) -> None:
    structure = load_structure_from_file_path(path)
    remove_water_and_metal_ions(structure)

    with TemporaryDirectory() as directory:
        cleaned_path = Path(directory) / "cleaned.pdb"
        structure.write_pdb(str(cleaned_path))

        model = load_model(str(cleaned_path), "A")
        exchange_single_base(model["1"], "A")
        write_model(model, "1EHZ_A1A.pdb")


def frustrate_rna(structure: Structure) -> float:
    energy = 0.0
    return energy
