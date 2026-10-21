import argparse

from rna_frustratometer.structures import read_structure


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyse RNA frustration")
    parser.add_argument("path", help="Path to a PDB or mmCIF file")
    args = parser.parse_args()

    structure = read_structure(args.path)
    for model in structure:
        for chain in model:
            print(f"Chain {chain.name}: {len(chain)} residues")


if __name__ == "__main__":
    main()
