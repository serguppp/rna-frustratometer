from pathlib import Path

from gemmi import Structure, cif, make_structure_from_block, read_pdb_string, read_structure

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


def frustrate_rna(structure: Structure) -> float:
    energy = 0.0
    return energy
