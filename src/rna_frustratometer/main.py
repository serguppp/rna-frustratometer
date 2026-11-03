import argparse
from pathlib import Path

from .structures import (
    frustrate_rna,
    load_structure_from_file_path,
    mutate_rna,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="RNA Frustratometer")
    parser.add_argument("file", type=Path, help="Path to a PDB or mmCIF file")
    parser.add_argument("--mutate", action="store_true", help="Generate a mutant")
    parser.add_argument(
        "--frustrate",
        action="store_true",
        help="Analyze the input structure",
    )
    args = parser.parse_args()

    path = args.file.resolve()
    if not path.is_file():
        parser.error(f"File does not exist: {path}")

    structure = load_structure_from_file_path(path)
    print(f"Loaded: {path.name}")
    print(f"Models: {len(structure)}")

    if args.mutate:
        if path.suffix.lower() != ".pdb":
            parser.error("The current mutation function requires a PDB file.")

        mutate_rna(path)
        print("Mutant saved to: 1EHZ_A1A.pdb")

    if args.frustrate:
        energy = frustrate_rna(structure)
        print(f"Input structure score: {energy}")


if __name__ == "__main__":
    main()
