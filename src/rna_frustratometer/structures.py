from pathlib import Path

from gemmi import Structure, cif, make_structure_from_block, read_pdb_string, read_structure

from .formats import SupportedFormats

BASES = ["A", "C", "G", "U"]


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
