# Adapted from serguppp/simrna-adapter (31be7d6); see LICENSE-SIMRNA-ADAPTER.

"""Construct SimRNA sequence variants in the reference glycosidic frame."""

from pathlib import Path

import gemmi

from .frustration import BASES, PairKey

Location = tuple[int, int]
MutationState = tuple[tuple[int, str], ...]
Frame = tuple[gemmi.Position, tuple[gemmi.Vec3, gemmi.Vec3, gemmi.Vec3]]

BASE_BEADS = {
    "A": ("N9", "C2", "C6"),
    "G": ("N9", "C2", "C6"),
    "C": ("N1", "C2", "C4"),
    "U": ("N1", "C2", "C4"),
}
BACKBONE_ATOMS = frozenset(
    (
        "P",
        "OP1",
        "OP2",
        "OP3",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "O2'",
        "C1'",
    )
)
PURINE_RING = frozenset(("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"))
PYRIMIDINE_RING = frozenset(("N1", "C2", "N3", "C4", "C5", "C6"))
BASE_ATOMS = {
    "A": PURINE_RING | {"N6"},
    "G": PURINE_RING | {"O6", "N2"},
    "C": PYRIMIDINE_RING | {"O2", "N4"},
    "U": PYRIMIDINE_RING | {"O2", "O4"},
}
EXOCYCLIC_PARENT = {"N6": "C6", "O6": "C6", "N2": "C2", "O2": "C2", "N4": "C4", "O4": "C4"}


def atom_position(residue: gemmi.Residue, name: str) -> gemmi.Position:
    atoms = [atom for atom in residue if atom.name == name]
    if len(atoms) != 1 or atoms[0].has_altloc() or atoms[0].occ <= 0:
        raise ValueError(f"Expected one occupied {name} without altloc at {residue.seqid}.")
    return atoms[0].pos


def glycosidic_frame(residue: gemmi.Residue) -> Frame:
    """Origin at N9/N1. The sugar sets x and the χ-defining atom fixes the plane."""
    nitrogen, edge = ("N9", "C4") if residue.name in ("A", "G") else ("N1", "C2")
    origin = atom_position(residue, nitrogen)
    bond = atom_position(residue, "C1'") - origin
    if bond.length() < 1e-6:
        raise ValueError(f"Degenerate glycosidic bond at {residue.seqid}.")
    x = bond.normalized()
    normal = x.cross(atom_position(residue, edge) - origin)
    if normal.length() < 1e-6:
        raise ValueError(f"Degenerate glycosidic frame at {residue.seqid}.")
    z = normal.normalized()
    return origin, (x, z.cross(x), z)


def replace_base(
    residue: gemmi.Residue, template: gemmi.Residue, *, preserve_same_class_ring: bool = False
) -> gemmi.Residue:
    replacement = residue.clone()
    if residue.name == template.name:
        return replacement

    origin, axes = glycosidic_frame(residue)
    template_origin, template_axes = glycosidic_frame(template)
    for index in reversed(range(len(replacement))):
        if replacement[index].name not in BACKBONE_ATOMS:
            del replacement[index]
    replacement.name = template.name
    for template_atom in template:
        if template_atom.name not in BASE_ATOMS[template.name]:
            continue
        displacement = template_atom.pos - template_origin
        xyz = gemmi.Vec3(origin.x, origin.y, origin.z)
        for axis, template_axis in zip(axes, template_axes, strict=True):
            xyz = xyz + axis * displacement.dot(template_axis)
        atom = template_atom.clone()
        atom.pos = gemmi.Position(xyz.x, xyz.y, xyz.z)
        replacement.add_atom(atom)

    native_purine = residue.name in ("A", "G")
    if preserve_same_class_ring and native_purine == (template.name in ("A", "G")):
        ring = PURINE_RING if native_purine else PYRIMIDINE_RING
        aligned = {atom.name: atom.pos for atom in replacement}
        native = {atom.name: atom.pos for atom in residue}
        for atom in replacement:
            if atom.name in ring and atom.name in native:
                atom.pos = native[atom.name]
            elif atom.name in EXOCYCLIC_PARENT:
                parent = EXOCYCLIC_PARENT[atom.name]
                atom.pos = native[parent] + (aligned[atom.name] - aligned[parent])
    return replacement


def read_templates(path: Path) -> tuple[dict[str, gemmi.Residue], dict[str, tuple[str, str]]]:
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("The SimRNA template must contain one structural model.")
    templates = {}
    identifiers = {}
    for chain in structure[0]:
        for residue in chain:
            base = residue.name
            if base in BASES and base not in templates:
                templates[base] = residue.clone()
                identifiers[base] = (chain.name, str(residue.seqid))
    if set(templates) != set(BASES):
        raise ValueError("The SimRNA template must contain A, C, G and U.")
    for residue in templates.values():
        glycosidic_frame(residue)
        for name in BASE_BEADS[residue.name]:
            atom_position(residue, name)
    return templates, identifiers


def select_contacts(structure: gemmi.Structure, cutoff: float) -> list[PairKey]:
    """Choose mutation targets by distance between native base midpoints."""
    residues = [(chain.name, str(r.seqid), r) for chain in structure[0] for r in chain]
    centers = []
    for _, _, residue in residues:
        start, _, end = BASE_BEADS[residue.name]
        centers.append((atom_position(residue, start) + atom_position(residue, end)) / 2.0)
    return [
        (chain_i, residue_i, chain_j, residue_j)
        for i, (chain_i, residue_i, _) in enumerate(residues)
        for j, (chain_j, residue_j, _) in enumerate(residues[i + 1 :], i + 1)
        if centers[i].dist(centers[j]) <= cutoff
    ]


def mutation_state(i: int, j: int, a: str, b: str, native: list[str]) -> MutationState:
    """Identify the same single or double mutation wherever it occurs."""
    return tuple(sorted((index, base) for index, base in ((i, a), (j, b)) if base != native[index]))


def write_variant(
    structure: gemmi.Structure,
    locations: list[Location],
    replacements: dict[tuple[int, str], gemmi.Residue],
    state: MutationState,
    path: Path,
) -> None:
    mutant = structure.clone()
    for index, base in state:
        chain_index, residue_index = locations[index]
        chain = mutant[0][chain_index]
        del chain[residue_index]
        chain.add_residue(replacements[index, base].clone(), pos=residue_index)
    mutant.write_pdb(str(path))
