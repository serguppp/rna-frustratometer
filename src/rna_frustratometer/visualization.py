"""Static figures for comparing frustration scores."""

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import gemmi
import matplotlib.pyplot as plt
import numpy as np

from .frustration import PairKey

plt.switch_backend("Agg")


def contact_maps(
    reference: Path, scores: Mapping[str, Mapping[PairKey, float]], output: Path
) -> None:
    """Show pairs with defined F in both potentials, using a common color scale."""
    structure = gemmi.read_structure(str(reference))
    positions = {
        (chain, residue): i
        for i, (chain, residue) in enumerate(
            (chain.name, str(r.seqid)) for chain in structure[0] for r in chain
        )
    }
    common = set.intersection(*(set(values) for values in scores.values()))
    if not common:
        return
    limit = max(abs(values[key]) for values in scores.values() for key in common) or 1
    fig, axes = plt.subplots(
        1, len(scores), squeeze=False, figsize=(4 * len(scores) + 1, 4), layout="constrained"
    )
    for axis, (name, values) in zip(axes[0], scores.items(), strict=True):
        matrix = np.full((len(positions), len(positions)), np.nan)
        for ci, ri, cj, rj in common:
            i, j = sorted((positions[ci, ri], positions[cj, rj]))
            matrix[j, i] = values[ci, ri, cj, rj]
        image = axis.imshow(
            matrix,
            origin="lower",
            interpolation="none",
            cmap="RdBu",
            vmin=-limit,
            vmax=limit,
            extent=(0.5, len(positions) + 0.5, 0.5, len(positions) + 0.5),
        )
        axis.set(title=name, xlabel="Nucleotide index i", ylabel="Nucleotide index j")
    fig.colorbar(image, ax=list(axes[0]), label="Frustration index F")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def homolog_profiles(
    rows: Iterable[tuple[int, int, float, float]],
    labels: Sequence[str],
    output: Path,
    *,
    title: str,
) -> None:
    """Compare mean incident F per alignment column over identical matched contacts."""
    profiles: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for i, j, fa, fb in rows:
        for column in (i, j):
            profiles[column].append((fa, fb))
    if not profiles:
        return
    columns = sorted(profiles)
    means = [np.mean(profiles[c], axis=0) for c in columns]
    fig, (axis, support) = plt.subplots(
        2,
        1,
        figsize=(7, 4.2),
        sharex=True,
        layout="constrained",
        gridspec_kw={"height_ratios": [3, 1]},
    )
    for index, label in enumerate(labels):
        axis.plot(columns, [mean[index] for mean in means], ".", label=label)
    axis.axhline(0, color="0.65", linewidth=0.7)
    axis.set(title=title, ylabel="Mean incident F")
    axis.legend(frameon=False)
    support.bar(columns, [len(profiles[c]) for c in columns], color="0.5")
    support.set(xlabel="Rfam alignment column", ylabel="Contacts")
    fig.savefig(output, dpi=200)
    plt.close(fig)
