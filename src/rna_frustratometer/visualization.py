"""Visualize RNA contact frustration, spatial class profiles, and energy components."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from statistics import fmean, pstdev
from typing import TYPE_CHECKING

import numpy as np
import plotly.graph_objects as go
from gemmi import read_structure
from numpy.typing import NDArray
from plotly.subplots import make_subplots

if TYPE_CHECKING:
    from rna_frustratometer.pipeline import AnalysisResult

ResidueKey = tuple[str, str]
CLASS_NAMES = ("Minimally frustrated", "Neutral", "Highly frustrated")
CLASS_COLORS = ("#1b7837", "#c4c9d0", "#b2182b")


@dataclass(frozen=True)
class FrustrationPlotData:
    """Reference order, scored contacts, and spatial contact counts."""

    keys: tuple[ResidueKey, ...]
    labels: tuple[str, ...]
    edges: NDArray[np.int64]
    scores: NDArray[np.float64]
    counts: NDArray[np.int64]
    covered_contacts: int
    radius: float
    anchor_atom: str
    high_cutoff: float
    minimal_cutoff: float


def prepare_frustration_plot_data(
    result: AnalysisResult,
    reference_path: Path,
    *,
    radius: float = 5.0,
    anchor_atom: str = "C4'",
    high_cutoff: float = -1.0,
    minimal_cutoff: float = 0.78,
) -> FrustrationPlotData:
    """Count contact midpoints within radius of each RNA anchor atom.

    Thresholds are borrowed from protein frustration, not calibrated for RNA.
    Each unordered contact is counted once in each spatial neighborhood.
    Undefined F values are rejected rather than assigned a class.
    """
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("Radius must be finite and positive.")
    if not all(math.isfinite(x) for x in (high_cutoff, minimal_cutoff)):
        raise ValueError("Cutoffs must be finite.")
    if high_cutoff >= minimal_cutoff:
        raise ValueError("Highly frustrated cutoff must be below minimal cutoff.")

    structure = read_structure(str(reference_path))
    if len(structure) != 1:
        raise ValueError("Use the single-model reference.pdb from the analysis.")

    keys: list[ResidueKey] = []
    labels: list[str] = []
    coordinates: list[tuple[float, float, float]] = []

    for chain in structure[0]:
        for residue in chain:
            key = (chain.name, str(residue.seqid))
            atoms = [atom for atom in residue if atom.name == anchor_atom]
            if len(atoms) != 1:
                raise ValueError(f"Expected one {anchor_atom} atom at {key}.")
            atom = atoms[0]
            if atom.altloc != "\x00" or atom.occ <= 0:
                raise ValueError(f"Select an occupied atom conformation at {key}.")
            position = (atom.pos.x, atom.pos.y, atom.pos.z)
            if not all(math.isfinite(value) for value in position):
                raise ValueError(f"Non-finite coordinates at {key}.")
            keys.append(key)
            labels.append(f"{chain.name}:{residue.seqid} {residue.name}")
            coordinates.append(position)

    if not keys:
        raise ValueError("Reference contains no residues.")
    index = {key: i for i, key in enumerate(keys)}
    if len(index) != len(keys):
        raise ValueError("Duplicate residue identifiers in reference.")

    edges: list[tuple[int, int]] = []
    scores: list[float] = []
    seen: set[tuple[int, int]] = set()
    for contact in result.contacts:
        i = index[(contact.chain_i, contact.residue_i)]
        j = index[(contact.chain_j, contact.residue_j)]
        edge = (min(i, j), max(i, j))
        if i == j or edge in seen:
            raise ValueError(f"Self-contact or duplicate contact: {edge}.")
        score = contact.frustration
        if score is None or not math.isfinite(score):
            raise ValueError(f"Cannot classify undefined frustration: {edge}.")
        seen.add(edge)
        edges.append(edge)
        scores.append(score)

    xyz: NDArray[np.float64] = np.asarray(coordinates, dtype=np.float64)
    edge_array: NDArray[np.int64] = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    score_array: NDArray[np.float64] = np.asarray(scores, dtype=np.float64)
    midpoints = (xyz[edge_array[:, 0]] + xyz[edge_array[:, 1]]) / 2
    squared_distances = np.sum((xyz[:, None, :] - midpoints[None, :, :]) ** 2, axis=2)
    nearby = squared_distances < radius**2

    masks = (
        score_array >= minimal_cutoff,
        (score_array > high_cutoff) & (score_array < minimal_cutoff),
        score_array <= high_cutoff,
    )
    counts = np.column_stack([np.sum(nearby[:, mask], axis=1, dtype=np.int64) for mask in masks])
    return FrustrationPlotData(
        keys=tuple(keys),
        labels=tuple(labels),
        edges=edge_array,
        scores=score_array,
        counts=counts,
        covered_contacts=int(np.sum(np.any(nearby, axis=0))),
        radius=radius,
        anchor_atom=anchor_atom,
        high_cutoff=high_cutoff,
        minimal_cutoff=minimal_cutoff,
    )


def _apply_style(
    figure: go.Figure, data: FrustrationPlotData, title: str, *, height: int = 450
) -> None:
    step = max(1, math.ceil(len(data.keys) / 24))
    ticks = list(range(0, len(data.keys), step))
    figure.update_layout(
        title={"text": title, "x": 0.04, "font": {"size": 19}},
        template="plotly_white",
        paper_bgcolor="white",
        plot_bgcolor="white",
        font={"family": "Arial, sans-serif", "size": 13, "color": "#2a3f5f"},
        width=1050,
        height=height,
        margin={"l": 80, "r": 35, "t": 110, "b": 140},
        legend={"orientation": "h", "x": 0, "y": 1.16, "font": {"size": 12}},
        hoverlabel={"bgcolor": "white", "font_size": 13},
    )
    figure.update_xaxes(
        title="Residue (PDB order)",
        tickmode="array",
        tickvals=[i + 1 for i in ticks],
        ticktext=[data.labels[i].split()[0] for i in ticks],
        tickangle=-45,
        range=[0.5, len(data.keys) + 0.5],
        gridcolor="#edf1f5",
        zeroline=False,
    )
    figure.update_yaxes(gridcolor="#edf1f5", zerolinecolor="#dbe3ed")
    for i in range(1, len(data.keys)):
        if data.keys[i][0] != data.keys[i - 1][0]:
            figure.add_vline(x=i + 0.5, line_color="#94a3b8", line_dash="dash")


def _add_profile_note(figure: go.Figure, data: FrustrationPlotData) -> None:
    figure.add_annotation(
        x=0,
        y=-0.38,
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        align="left",
        font={"size": 11, "color": "#64748b"},
        text=(
            f"{data.radius:g} Å around {data.anchor_atom}; contact midpoints. "
            f"Coverage: {data.covered_contacts}/{len(data.scores)} contacts."
            f"<br>Protein cutoffs: F ≤ {data.high_cutoff:g}, "
            f"F ≥ {data.minimal_cutoff:g}; exploratory for RNA."
        ),
    )


def _profile_x(data: FrustrationPlotData) -> list[int | None]:
    """Insert gaps so line traces do not join separate chains."""
    positions: list[int | None] = []
    for i in range(len(data.keys)):
        if i > 0 and data.keys[i][0] != data.keys[i - 1][0]:
            positions.append(None)
        positions.append(i + 1)
    return positions


def plot_contact_map_classic(data: FrustrationPlotData, *, color_limit: float = 3.0) -> go.Figure:
    """Keep the triangular map and thresholds, using the existing Plotly style."""
    if not math.isfinite(color_limit) or color_limit <= 0:
        raise ValueError("Color limit must be finite and positive.")
    if not -color_limit < data.high_cutoff < data.minimal_cutoff < color_limit:
        raise ValueError("Color range must contain both classification cutoffs.")

    size = len(data.keys)
    matrix: list[list[float | None]] = [[None] * size for _ in range(size)]
    text = [[""] * size for _ in range(size)]
    for (i, j), score in zip(data.edges, data.scores, strict=True):
        matrix[j][i] = float(score)
        text[j][i] = f"{data.labels[i]} ↔ {data.labels[j]}"

    lower = (data.high_cutoff + color_limit) / (2 * color_limit)
    upper = (data.minimal_cutoff + color_limit) / (2 * color_limit)
    figure = go.Figure(
        go.Heatmap(
            x=list(range(1, size + 1)),
            y=list(range(1, size + 1)),
            z=matrix,
            text=text,
            zmin=-color_limit,
            zmax=color_limit,
            zauto=False,
            colorscale=[
                [0, CLASS_COLORS[2]],
                [lower, CLASS_COLORS[1]],
                [upper, CLASS_COLORS[1]],
                [1, CLASS_COLORS[0]],
            ],
            colorbar={"title": "F (BRiQ)", "thickness": 14, "outlinewidth": 0},
            xgap=1,
            ygap=1,
            hoverongaps=False,
            hovertemplate="%{text}<br>F=%{z:.3f}<extra></extra>",
        )
    )
    _apply_style(figure, data, "Mutational frustration — contact map", height=800)
    figure.update_layout(width=850, margin={"l": 85, "r": 95, "t": 85, "b": 135})
    figure.update_xaxes(title="Residue i", showgrid=False)
    step = max(1, math.ceil(size / 24))
    tick_indices = list(range(0, size, step))

    figure.update_yaxes(
        title="Residue j",
        scaleanchor="x",
        scaleratio=1,
        tickmode="array",
        tickvals=[i + 1 for i in tick_indices],
        ticktext=[data.labels[i].split()[0] for i in tick_indices],
        range=[0.5, size + 0.5],
        showgrid=False,
    )
    figure.add_annotation(
        x=0,
        y=-0.2,
        xref="paper",
        yref="paper",
        xanchor="left",
        showarrow=False,
        font={"size": 11, "color": "#64748b"},
        text="White: no scored contact. Gray: neutral range under exploratory protein cutoffs.",
    )
    return figure


def plot_contact_density(data: FrustrationPlotData) -> go.Figure:
    """Plot exactly the same spatial counts with interactive residue labels."""
    figure = go.Figure()
    total = data.counts.sum(axis=1)
    x = _profile_x(data)
    for column, (name, color) in enumerate(
        zip((*CLASS_NAMES, "Total"), (*CLASS_COLORS, "#334155"), strict=True)
    ):
        y: list[int | None] = []
        details: list[list[str | int] | None] = []
        for i, label in enumerate(data.labels):
            if i > 0 and data.keys[i][0] != data.keys[i - 1][0]:
                y.append(None)
                details.append(None)
            y.append(int(total[i] if column == 3 else data.counts[i, column]))
            details.append([label, int(total[i])])
        figure.add_trace(
            go.Scatter(
                x=x,
                y=y,
                customdata=details,
                mode="lines+markers",
                name=name,
                line={"color": color, "width": 2, "dash": "dot" if column == 3 else "solid"},
                marker={"color": color, "size": 5},
                connectgaps=False,
                hovertemplate=(
                    "%{customdata[0]}<br>Count=%{y}<br>Total in sphere=%{customdata[1]}"
                    "<extra>%{fullData.name}</extra>"
                ),
            )
        )
    _apply_style(figure, data, "Contact density along sequence", height=480)
    figure.update_yaxes(
        title=f"Contacts within {data.radius:g} Å sphere", dtick=1, rangemode="tozero"
    )
    _add_profile_note(figure, data)
    return figure


def plot_contact_fractions(data: FrustrationPlotData) -> go.Figure:
    """Show the same class fractions with counts available on hover."""
    figure = go.Figure()
    total = data.counts.sum(axis=1)
    for column, (name, color) in enumerate(zip(CLASS_NAMES, CLASS_COLORS, strict=True)):
        fractions = [
            float(count / n) if n else None
            for count, n in zip(data.counts[:, column], total, strict=True)
        ]
        details = [
            [label, int(count), int(n)]
            for label, count, n in zip(data.labels, data.counts[:, column], total, strict=True)
        ]
        figure.add_trace(
            go.Bar(
                x=list(range(1, len(data.keys) + 1)),
                y=fractions,
                name=name,
                marker_color=color,
                marker_line_width=0,
                width=0.96,
                customdata=details,
                hovertemplate=(
                    "%{customdata[0]}<br>Fraction=%{y:.1%}"
                    "<br>Contacts=%{customdata[1]} / %{customdata[2]}"
                    "<extra>%{fullData.name}</extra>"
                ),
            )
        )
    _apply_style(figure, data, "Fraction of contacts in each frustration class", height=480)
    figure.update_layout(barmode="stack", bargap=0)
    figure.update_yaxes(title="Fraction of contacts", range=[0, 1], tickformat=".0%", dtick=0.2)
    _add_profile_note(figure, data)
    return figure


PairKey = tuple[str, str, str, str]
COMPONENTS = (
    "E_bb",
    "E_bo",
    "E_oo",
    "E_clash",
    "E_base_ribose",
    "E_base_phosphate",
    "E_ribose_phosphate",
    "E_rotamer",
    "E_phosphate_local",
)


def frustration(native: float, variants: list[float]) -> float | None:
    sigma = pstdev(variants)
    return None if sigma == 0 else (fmean(variants) - native) / sigma


def plot_energy_components(path: Path) -> go.Figure:
    """Show reference energies, contributions to F, and clash sensitivity."""
    with path.open(newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))

    references = [row for row in rows if row["state"] == "reference"]
    if len(references) != 1:
        raise ValueError("Expected exactly one reference energy row.")

    reference = references[0]
    native = float(reference["E_total"])

    for row in rows:
        values = [float(row[name]) for name in (*COMPONENTS, "E_total")]

        if not all(math.isfinite(value) for value in values):
            raise ValueError("Non-finite energy component.")

        if not math.isclose(sum(values[:-1]), values[-1], rel_tol=1e-10, abs_tol=1e-6):
            raise ValueError("Components do not sum to E_total.")

    groups: dict[PairKey, list[dict[str, str]]] = defaultdict(list)

    for row in rows:
        if row["state"] == "variant":
            key = (
                row["chain_i"],
                row["residue_i"],
                row["chain_j"],
                row["residue_j"],
            )
            groups[key].append(row)

    labels: list[str] = []
    total_scores: list[float | None] = []
    no_clash_scores: list[float | None] = []
    contributions: dict[str, list[float | None]] = {name: [] for name in COMPONENTS}
    mean_deltas: dict[str, list[float]] = {name: [] for name in COMPONENTS}

    expected = set(product("ACGU", repeat=2))

    for key, variants in groups.items():
        combinations = {(row["base_i"], row["base_j"]) for row in variants}
        if len(variants) != 16 or combinations != expected:
            raise ValueError(f"Expected all 16 combinations for contact {key}.")

        chain_i, residue_i, chain_j, residue_j = key
        labels.append(f"{chain_i}:{residue_i} – {chain_j}:{residue_j}")

        totals = [float(row["E_total"]) for row in variants]
        sigma = pstdev(totals)
        total_scores.append(frustration(native, totals))

        for name in COMPONENTS:
            delta = fmean(float(row[name]) for row in variants) - float(reference[name])
            mean_deltas[name].append(delta)
            contributions[name].append(None if sigma == 0 else delta / sigma)

        totals_without_clash = [float(row["E_total"]) - float(row["E_clash"]) for row in variants]
        native_without_clash = native - float(reference["E_clash"])
        no_clash_scores.append(frustration(native_without_clash, totals_without_clash))

    figure = make_subplots(
        rows=2,
        cols=1,
        vertical_spacing=0.18,
        subplot_titles=(
            f"Reference energy components — E_total = {native:.3f}",
            "Contributions to F and sensitivity to E_clash",
        ),
    )

    figure.add_trace(
        go.Bar(
            x=list(COMPONENTS),
            y=[float(reference[name]) for name in COMPONENTS],
            marker_color=["#d62728" if name == "E_clash" else "#4078a8" for name in COMPONENTS],
            showlegend=False,
            hovertemplate="%{x}<br>E=%{y:.4f}<extra></extra>",
        ),
        row=1,
        col=1,
    )

    for name in COMPONENTS:
        figure.add_trace(
            go.Bar(
                x=labels,
                y=contributions[name],
                name=name,
                customdata=mean_deltas[name],
                hovertemplate=(
                    "%{x}<br>Contribution to F=%{y:.4f}"
                    "<br>Mean ΔE=%{customdata:.4f}<extra>%{fullData.name}</extra>"
                ),
            ),
            row=2,
            col=1,
        )

    for name, scores, color, symbol in (
        ("F: full energy", total_scores, "black", "diamond"),
        ("F: excluding E_clash", no_clash_scores, "#9333ea", "x"),
    ):
        figure.add_trace(
            go.Scatter(
                x=labels,
                y=scores,
                mode="markers",
                marker={"color": color, "symbol": symbol, "size": 9},
                name=name,
                hovertemplate="%{x}<br>F=%{y:.4f}<extra>%{fullData.name}</extra>",
            ),
            row=2,
            col=1,
        )

    figure.update_layout(
        template="plotly_white",
        barmode="relative",
        height=950,
        width=1300,
    )
    figure.update_yaxes(title_text="Energy [BRiQ units]", row=1, col=1)
    figure.update_yaxes(title_text="Contribution / F", row=2, col=1)
    figure.update_xaxes(title_text="Contact", row=2, col=1)
    return figure


def plot_contact_map(
    result: AnalysisResult,
    reference_path: Path,
    *,
    color_limit: float = 3.0,
) -> go.Figure:
    """Keep the existing plot_contact_map(result, reference_path) interface."""
    data = prepare_frustration_plot_data(result, reference_path)
    return plot_contact_map_classic(data, color_limit=color_limit)


def create_frustration_plots(
    result: AnalysisResult,
    reference_path: Path,
    *,
    components_path: Path | None = None,
    radius: float = 5.0,
    high_cutoff: float = -1.0,
    minimal_cutoff: float = 0.78,
    color_limit: float = 3.0,
) -> dict[str, go.Figure]:
    """Prepare data once and create the contact plots and optional energy report.

    Plotting never reruns BRiQ. Radius affects only the spatial profiles.
    Frustration values remain those returned by the analysis pipeline.
    """
    data = prepare_frustration_plot_data(
        result,
        reference_path,
        radius=radius,
        high_cutoff=high_cutoff,
        minimal_cutoff=minimal_cutoff,
    )
    figures = {
        "contact_map": plot_contact_map_classic(data, color_limit=color_limit),
        "contact_density": plot_contact_density(data),
        "contact_fractions": plot_contact_fractions(data),
    }
    if components_path is not None:
        figures["energy_components"] = plot_energy_components(components_path)
    return figures
