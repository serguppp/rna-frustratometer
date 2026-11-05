from pathlib import Path

from rna_frustratometer.pipeline import RNAfrustratometer

data = ["1ANR", "1NC0", "2KOC"]
for pdb_id in data:
    analysis = RNAfrustratometer(
        input_dir=Path(f"tests/data/{pdb_id}.pdb"),
        output_dir=Path(f"results/{pdb_id}_cache_01"),
    )

    result = analysis.run()

    for contact in result.contacts:
        print(
            f"{contact.chain_i}:{contact.residue_i} {contact.chain_j}:{contact.residue_j}",
            contact.frustration,
        )
