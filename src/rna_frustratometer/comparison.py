"""Compare potentials at the level of scored residue pairs."""

import csv
from pathlib import Path
from typing import Literal

from .frustration import PAIR_COLUMNS, STRATEGIES, PairKey, canonical_pair

Potential = Literal["cgrnasp", "simrna"]
POTENTIALS: tuple[Potential, ...] = ("cgrnasp", "simrna")


def read_scores(path: Path) -> dict[PairKey, dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        records: dict[PairKey, dict[str, str]] = {}
        for row in csv.DictReader(stream):
            key = canonical_pair(
                (row["chain_i"], row["residue_i"], row["chain_j"], row["residue_j"])
            )
            if key in records:
                raise ValueError(f"Duplicate pair {key} in {path}.")
            records[key] = row
    return records


def write_shared_scores(cgrnasp: Path, simrna: Path, output: Path) -> int:
    """Join the two score tables by pair identifier. Blank F stays undefined."""
    cgr, sim = read_scores(cgrnasp), read_scores(simrna)
    shared = sorted(cgr.keys() & sim.keys())
    columns = [
        f"F_{strategy}_{weight}" for strategy in STRATEGIES for weight in ("composition", "uniform")
    ]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([*PAIR_COLUMNS, *(f"{p}_{c}" for p in POTENTIALS for c in columns)])
        for key in shared:
            writer.writerow([*key, *(table[key][c] for table in (cgr, sim) for c in columns)])
    return len(shared)
