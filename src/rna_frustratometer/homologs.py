"""Compare homologs using saved energies and a precomputed Rfam/3D alignment."""

import csv
import hashlib
import json
from collections import defaultdict
from collections.abc import Collection
from itertools import groupby
from pathlib import Path
from statistics import StatisticsError, correlation, median

from .comparison import POTENTIALS, Potential
from .frustration import (
    STRATEGIES,
    UNIFORM_WEIGHTS,
    BaseFrequencies,
    ContactEnergies,
    PairKey,
    Strategy,
    canonical_pair,
    read_contact_energies,
)
from .simrna import PAIR_MINIMUM_SD
from .visualization import homolog_profiles

Scan = dict[PairKey, dict[Strategy, ContactEnergies]]


def _spearman(a: list[float], b: list[float]) -> float | None:
    try:
        return correlation(a, b, method="ranked")
    except StatisticsError:
        return None


def _load_scans(
    directory: Path, members: dict[str, dict[str, str]]
) -> dict[str, dict[Potential, Scan]]:
    scans: dict[str, dict[Potential, Scan]] = defaultdict(dict)
    for name, record in members.items():
        folder = directory / name
        cg_dir = folder / "cgrnasp"
        if not cg_dir.is_dir():  # Layout of the existing homolog experiment.
            cg_dir = folder / "independent/cgrnasp"
            folder = folder / "independent"
        sim_dir = folder / "simrna"
        sim_run = sim_dir / "simrna/run.json"  # New layout.
        if not sim_run.is_file():
            sim_run = sim_dir / "run.json"
        metadata = json.loads(sim_run.read_text(encoding="utf-8"))
        reference_hashes = {metadata["reference_sha256"]}
        reference = sim_dir / "reference.pdb"
        if reference.is_file():
            if hashlib.sha256(reference.read_bytes()).hexdigest() not in reference_hashes:
                raise ValueError(f"SimRNA reference differs from the recorded input: {name}.")
            input_metadata = sim_dir / "input.json"
            if input_metadata.is_file():
                reference_hashes.add(json.loads(input_metadata.read_text())["sha256"])
        if record["prepared_sha256"] not in reference_hashes:
            raise ValueError(f"Saved energies do not match the panel input: {name}.")
        for potential, subdir in zip(POTENTIALS, (cg_dir, sim_dir), strict=True):
            scans[name][potential] = {
                canonical_pair(key): value
                for key, value in read_contact_energies(subdir / "contact_variants.csv").items()
            }
    return scans


def _aligned_contacts(
    scans: dict[Potential, Scan], alignment: dict[str, int]
) -> dict[tuple[int, int], PairKey]:
    shared = scans["cgrnasp"].keys() & scans["simrna"].keys()
    if any(ci != "A" or cj != "A" for ci, _, cj, _ in shared):
        raise ValueError("This panel alignment expects prepared, single-chain PDBs (chain A).")
    result = {}
    for key in shared:
        _, ri, _, rj = key
        if ri in alignment and rj in alignment:
            i, j = sorted((alignment[ri], alignment[rj]))
            result[i, j] = key
    return result


def analyze_homologs(
    panel_file: Path,
    alignment_file: Path,
    scans_dir: Path,
    output: Path,
    *,
    plot_families: Collection[str] = (),
) -> None:
    """Save matched contact values and descriptive statistics for each homolog pair."""
    members = {row["name_id"]: row for row in json.loads(panel_file.read_text(encoding="utf-8"))}
    scans = _load_scans(scans_dir, members)
    sequences: dict[str, list[str]] = defaultdict(list)
    for record in members.values():
        sequences[record["rfam"]].append(record["sequence_observed"])
    family_weights = {
        family: BaseFrequencies.from_sequence("".join(seq)).weights()
        for family, seq in sequences.items()
    }
    with alignment_file.open(newline="", encoding="utf-8-sig") as stream:
        alignment = sorted(
            (
                row
                for row in csv.DictReader(stream)
                if row["structure_a"] in members and row["structure_b"] in members
            ),
            key=lambda row: row["comparison"],
        )
    if not alignment:
        raise ValueError("No panel members match the supplied alignment.")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Use a new or empty homolog output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    (output / "inputs.json").write_text(
        json.dumps(
            {
                "panel": str(panel_file.resolve()),
                "panel_sha256": hashlib.sha256(panel_file.read_bytes()).hexdigest(),
                "alignment": str(alignment_file.resolve()),
                "alignment_sha256": hashlib.sha256(alignment_file.read_bytes()).hexdigest(),
                "scans": str(scans_dir.resolve()),
                "family_weights_acgu_cross_acgu": family_weights,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    labels = [
        "family",
        "structure_a",
        "structure_b",
        "scope",
        "potential",
        "strategy",
        "background",
    ]
    with (
        (output / "matched_contacts.csv").open("w", newline="", encoding="utf-8") as detail_file,
        (output / "homolog_statistics.csv").open("w", newline="", encoding="utf-8") as summary_file,
    ):
        details, summary = csv.writer(detail_file), csv.writer(summary_file)
        details.writerow(
            [*labels, "column_i", "column_j", "F_a", "F_b", "native_E_pair_a", "native_E_pair_b"]
        )
        summary.writerow(
            [
                *labels,
                "matched_pairs",
                "defined_pairs",
                "rho_F",
                "median_abs_delta_F",
                "rho_native_E_pair",
            ]
        )
        for comparison, group in groupby(alignment, key=lambda row: row["comparison"]):
            rows = list(group)
            first = rows[0]
            a, b, family = first["structure_a"], first["structure_b"], first["rfam"]
            scope_rows = {
                "all_mapped": rows,
                "reliable": [r for r in rows if r["USalign_agrees"] == "True"],
                "core_3A": [
                    r
                    for r in rows
                    if r["family_core3A"] == "True" and r["USalign_agrees"] == "True"
                ],
                "conserved": [
                    r
                    for r in rows
                    if r["identical_base"] == "True" and r["USalign_agrees"] == "True"
                ],
            }
            for scope, subset in scope_rows.items():
                mapping_a = {r["a_residue"]: int(r["alignment_column"]) for r in subset}
                mapping_b = {r["b_residue"]: int(r["alignment_column"]) for r in subset}
                ca = _aligned_contacts(scans[a], mapping_a)
                cb = _aligned_contacts(scans[b], mapping_b)
                shared = sorted(ca.keys() & cb.keys())
                for strategy in STRATEGIES:
                    # Exclude undefined comparisons before splitting by potential or weights.
                    backgrounds = {"family": family_weights[family], "uniform": UNIFORM_WEIGHTS}
                    pairs = []
                    for columns in shared:
                        scores_by_background = {}
                        for pot in POTENTIALS:
                            floor = (
                                PAIR_MINIMUM_SD if pot == "simrna" and strategy == "pair" else 0.0
                            )
                            left, right = scans[a][pot][ca[columns]], scans[b][pot][cb[columns]]
                            for background, weights in backgrounds.items():
                                fa = (
                                    left[strategy].statistics(weights, minimum_sd=floor).frustration
                                )
                                fb = (
                                    right[strategy]
                                    .statistics(weights, minimum_sd=floor)
                                    .frustration
                                )
                                if fa is not None and fb is not None:
                                    scores_by_background[background, pot] = (
                                        *columns,
                                        fa,
                                        fb,
                                        left["pair"].native,
                                        right["pair"].native,
                                    )
                        if len(scores_by_background) == len(POTENTIALS) * len(backgrounds):
                            pairs.append(scores_by_background)
                    for background in backgrounds:
                        for pot in POTENTIALS:
                            values = [pair[background, pot] for pair in pairs]
                            label = (family, a, b, scope, pot, strategy, background)
                            details.writerows((*label, *row) for row in values)
                            left_f, right_f = [r[2] for r in values], [r[3] for r in values]
                            rho = _spearman(left_f, right_f)
                            summary.writerow(
                                [
                                    *label,
                                    len(shared),
                                    len(values),
                                    rho,
                                    median(abs(x - y) for x, y in zip(left_f, right_f, strict=True))
                                    if values
                                    else None,
                                    _spearman([r[4] for r in values], [r[5] for r in values]),
                                ]
                            )
                            if (
                                family in plot_families
                                and scope == "reliable"
                                and background == "family"
                            ):
                                homolog_profiles(
                                    ((i, j, fa, fb) for i, j, fa, fb, _, _ in values),
                                    (a, b),
                                    output / f"{comparison}_{pot}_{strategy}.png",
                                    title=f"{family} · {pot} · F_{strategy} · {len(values)} pairs",
                                )
