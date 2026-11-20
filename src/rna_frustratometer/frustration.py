"""Exact population moments over all 16 RNA identity combinations."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import gemmi

BASES = ("A", "C", "G", "U")
COMBINATIONS: tuple[tuple[str, str], ...] = tuple((a, b) for a in BASES for b in BASES)
UNIFORM_WEIGHTS = (1 / 16,) * 16
PAIR_COLUMNS = ("chain_i", "residue_i", "chain_j", "residue_j")
PairKey = tuple[str, str, str, str]
Strategy = Literal["pair", "total"]
STRATEGIES: tuple[Strategy, ...] = ("pair", "total")
Weighting = Literal["composition", "uniform"]
WEIGHTINGS: tuple[Weighting, ...] = ("composition", "uniform")


@dataclass(frozen=True)
class BaseFrequencies:
    counts: tuple[int, int, int, int]
    source: str
    pseudocount: float = 0.0
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        if len(self.counts) != 4 or any(n < 0 for n in self.counts) or sum(self.counts) == 0:
            raise ValueError("Expected nonnegative A/C/G/U counts and at least one nucleotide.")
        if not math.isfinite(self.pseudocount) or self.pseudocount < 0:
            raise ValueError("Pseudocount must be finite and nonnegative.")

    @classmethod
    def from_sequence(
        cls, sequence: str, *, source: str = "sequence", pseudocount: float = 0.0
    ) -> BaseFrequencies:
        sequence = "".join(sequence.split()).upper()
        if not sequence or set(sequence) - set(BASES):
            raise ValueError("Provide a nonempty, unambiguous A/C/G/U sequence.")
        counts = Counter(sequence)
        return cls((counts["A"], counts["C"], counts["G"], counts["U"]), source, pseudocount)

    @classmethod
    def from_pdb(cls, path: Path, *, pseudocount: float = 0.0) -> BaseFrequencies:
        path = Path(path).resolve(strict=True)
        structure = gemmi.read_structure(str(path))
        if not len(structure):
            raise ValueError("Weight reference has no structural models.")
        sequence = []
        seen = set()
        # Count model index 0 once; do not count atoms, contacts or NMR copies.
        for chain in structure[0]:
            for residue in chain:
                if residue.name not in BASES:
                    raise ValueError(
                        f"Weight reference contains noncanonical residue {residue.name}."
                    )
                key = (chain.name, str(residue.seqid))
                if key in seen:
                    raise ValueError(f"Duplicate residue in weight reference: {key}.")
                seen.add(key)
                sequence.append(residue.name)
        background = cls.from_sequence("".join(sequence), source=str(path), pseudocount=pseudocount)
        return cls(
            background.counts,
            background.source,
            pseudocount,
            hashlib.sha256(path.read_bytes()).hexdigest(),
        )

    @property
    def probabilities(self) -> tuple[float, ...]:
        total = sum(self.counts) + 4 * self.pseudocount
        return tuple((n + self.pseudocount) / total for n in self.counts)

    def weights(
        self,
        combinations: Sequence[tuple[str, str]] = COMBINATIONS,
        *,
        mode: Weighting = "composition",
    ) -> tuple[float, ...]:
        if len(combinations) != 16 or set(combinations) != set(COMBINATIONS):
            raise ValueError("Expected exactly the 16 distinct A/C/G/U pairs.")
        if mode == "uniform":
            return UNIFORM_WEIGHTS
        if mode != "composition":
            raise ValueError("Weighting must be composition or uniform.")
        p = dict(zip(BASES, self.probabilities, strict=True))
        return tuple(p[a] * p[b] for a, b in combinations)

    def metadata(self) -> dict[str, object]:
        return {
            "source": self.source,
            "source_sha256": self.source_sha256,
            "model_index": 0,
            "length": sum(self.counts),
            "counts": dict(zip(BASES, self.counts, strict=True)),
            "probabilities": dict(zip(BASES, self.probabilities, strict=True)),
            "pseudocount": self.pseudocount,
            "combinations": [
                {"base_i": a, "base_j": b, "composition": w, "uniform": 1 / 16}
                for (a, b), w in zip(COMBINATIONS, self.weights(), strict=True)
            ],
        }


@dataclass(frozen=True)
class EnergyStatistics:
    native: float
    mean: float
    sd: float
    frustration: float | None


@dataclass(frozen=True)
class ContactEnergies:
    native: float
    variants: tuple[float, ...]

    def statistics(
        self, weights: Sequence[float] = UNIFORM_WEIGHTS, *, minimum_sd: float = 0.0
    ) -> EnergyStatistics:
        return energy_statistics(self.native, self.variants, weights, minimum_sd=minimum_sd)


def energy_statistics(
    native: float,
    energies: Sequence[float],
    weights: Sequence[float] = UNIFORM_WEIGHTS,
    *,
    minimum_sd: float = 0.0,
) -> EnergyStatistics:
    """Population variance, without Bessel or effective-sample-size correction."""
    if len(energies) != 16 or len(weights) != 16:
        raise ValueError("Expected 16 energies and 16 weights.")
    if not all(math.isfinite(v) for v in (native, *energies, *weights, minimum_sd)):
        raise ValueError("Energies, weights and minimum SD must be finite.")
    if minimum_sd < 0 or any(w < 0 for w in weights):
        raise ValueError("Weights and minimum SD must be nonnegative.")
    if not math.isclose(math.fsum(weights), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("Weights must sum to one.")
    deltas = [e - native for e in energies]
    mean_delta = math.fsum(w * d for w, d in zip(weights, deltas, strict=True))
    variance = math.fsum(w * (d - mean_delta) ** 2 for w, d in zip(weights, deltas, strict=True))
    sigma = math.sqrt(variance)
    score = mean_delta / sigma if sigma > minimum_sd else None
    return EnergyStatistics(native, native + mean_delta, sigma, score)


def read_variant_groups(path: Path) -> dict[PairKey, list[dict[str, str]]]:
    groups: dict[PairKey, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            key = (row["chain_i"], row["residue_i"], row["chain_j"], row["residue_j"])
            groups.setdefault(key, []).append(row)
    if not groups:
        raise ValueError(f"No variant energies in {path}.")
    seen = set()
    for key, rows in groups.items():
        left, right = sorted(((key[0], key[1]), (key[2], key[3])))
        if left == right or (left, right) in seen:
            raise ValueError(f"Duplicate/reversed/self contact: {key}.")
        seen.add((left, right))
        combinations = [(r["base_i"], r["base_j"]) for r in rows]
        if len(rows) != 16 or set(combinations) != set(COMBINATIONS):
            raise ValueError(f"Expected all 16 distinct variants for {key}.")
        rows.sort(key=lambda r: COMBINATIONS.index((r["base_i"], r["base_j"])))
    return groups


def read_contact_energies(path: Path) -> dict[PairKey, dict[Strategy, ContactEnergies]]:
    contacts: dict[PairKey, dict[Strategy, ContactEnergies]] = {}
    for key, rows in read_variant_groups(path).items():
        strategies: dict[Strategy, ContactEnergies] = {}
        for strategy in STRATEGIES:
            try:
                native_values = {float(r[f"native_E_{strategy}"]) for r in rows}
                energies = tuple(float(r[f"E_{strategy}"]) for r in rows)
            except (KeyError, ValueError) as error:
                raise ValueError(
                    f"Missing or invalid {strategy} energies for {key} in {path}."
                ) from error
            if len(native_values) != 1:
                raise ValueError(f"Inconsistent native {strategy} energy for {key}.")
            strategies[strategy] = ContactEnergies(native_values.pop(), energies)
        contacts[key] = strategies
    return contacts


def canonical_pair(key: PairKey) -> PairKey:
    left, right = sorted((key[:2], key[2:]))
    return left[0], left[1], right[0], right[1]


def save_frustration_scores(
    variants_file: Path,
    output_dir: Path,
    background: BaseFrequencies,
    *,
    weighting: Weighting = "composition",
    pair_minimum_sd: float = 0.0,
    energies: dict[PairKey, dict[Strategy, ContactEnergies]] | None = None,
) -> Path:
    """Reweight saved pair/total energies without invoking either scoring engine."""
    if weighting not in WEIGHTINGS:
        raise ValueError("Weighting must be composition or uniform.")
    if energies is None:
        energies = read_contact_energies(variants_file)
    output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(PAIR_COLUMNS)
    for strategy in STRATEGIES:
        fields.extend(
            [f"native_E_{strategy}", f"mean_E_{strategy}", f"sd_E_{strategy}", f"F_{strategy}"]
        )
        for mode in WEIGHTINGS:
            fields.extend(
                [f"mean_E_{strategy}_{mode}", f"sd_E_{strategy}_{mode}", f"F_{strategy}_{mode}"]
            )
    output = output_dir / "frustration_scores.csv"
    records = []
    for key, strategies in energies.items():
        record: dict[str, str | float] = dict(zip(PAIR_COLUMNS, key, strict=True))
        for strategy, contact in strategies.items():
            record[f"native_E_{strategy}"] = contact.native
            for mode in WEIGHTINGS:
                stats = contact.statistics(
                    background.weights(mode=mode),
                    minimum_sd=pair_minimum_sd if strategy == "pair" else 0,
                )
                values: dict[str, str | float] = {
                    "mean_E": stats.mean,
                    "sd_E": stats.sd,
                    "F": stats.frustration if stats.frustration is not None else "",
                }
                for name, value in values.items():
                    record[f"{name}_{strategy}_{mode}"] = value
                    if mode == weighting:
                        record[f"{name}_{strategy}"] = value
        records.append(record)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    (output_dir / "decoy_weights.json").write_text(
        json.dumps(
            {
                **background.metadata(),
                "primary_weighting": weighting,
                "variance": "weighted population variance",
                "native_combination_included": True,
                "pair_minimum_sd": pair_minimum_sd,
                "undefined_F": "empty CSV field",
                "energies_file": str(variants_file.resolve()),
                "energies_sha256": hashlib.sha256(variants_file.read_bytes()).hexdigest(),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return output
