"""Plotting helpers for multi-definition CloneDisp results."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .multi import MultiCloneDispResult


MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
FOCUS_COLORS = ["#E8A45A", "#6A4C93", "#4C9F70", "#D1495B"]
OTHER_COLOR = "#C5CBD3"
INK = "#1F1F1F"


def plot_multi_definition(
    result: MultiCloneDispResult,
    output_prefix: str | Path,
    *,
    state_order: Sequence[str] | None = None,
    focus_states: Sequence[str] = (),
    title: str = "Clone dispersion across definitions",
    figsize: tuple[float, float] = (4.3, 3.1),
    dpi: int = 400,
) -> plt.Figure:
    """Plot one patient-level mean CDI point per state and definition."""

    if result.definition_summary.empty:
        raise ValueError("No estimable definition-state metrics to plot")

    definitions = list(result.parameters["definitions"])
    if state_order is None:
        state_order = result.state_summary.sort_values(
            ["rank_median_cdi", "state"]
        )["state"].tolist()
    state_order = [state for state in state_order if state in set(result.state_summary["state"])]
    if not state_order:
        raise ValueError("state_order contains no states present in the result")

    lookup = result.definition_summary.set_index(["definition", "state"])
    medians = result.state_summary.set_index("state")["median_definition_cdi"]
    offsets = (
        np.linspace(-0.22, 0.22, len(definitions))
        if len(definitions) > 1
        else np.array([0.0])
    )

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": dpi,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.34, top=0.86)

    for state_index, state in enumerate(state_order):
        focal = state in focus_states
        color = (
            FOCUS_COLORS[list(focus_states).index(state) % len(FOCUS_COLORS)]
            if focal
            else OTHER_COLOR
        )
        for definition_index, definition in enumerate(definitions):
            key = (definition, state)
            if key not in lookup.index:
                continue
            value = float(lookup.loc[key, "mean_clone_dispersion"])
            ax.scatter(
                state_index + offsets[definition_index],
                value,
                marker=MARKERS[definition_index % len(MARKERS)],
                s=27 if focal else 21,
                facecolor=color,
                edgecolor=INK if focal else "white",
                linewidth=0.55 if focal else 0.3,
                zorder=4 if focal else 3,
            )
        median = float(medians.loc[state])
        ax.plot(
            [state_index - 0.26, state_index + 0.26],
            [median, median],
            color=INK,
            linewidth=1.2,
            solid_capstyle="round",
            zorder=5,
        )

    ax.set_xticks(np.arange(len(state_order)))
    ax.set_xticklabels(state_order)
    for label, state in zip(ax.get_xticklabels(), state_order):
        if state in focus_states:
            label.set_color(
                FOCUS_COLORS[list(focus_states).index(state) % len(FOCUS_COLORS)]
            )
            label.set_fontweight("bold")
    ax.set_xlim(-0.55, len(state_order) - 0.45)
    ax.set_xlabel("State")
    ax.set_ylabel("Clone-agnostic spread (CDI)")
    ax.set_title(title, loc="left", pad=5)
    ax.margins(x=0, y=0.08)

    handles = [
        Line2D(
            [0],
            [0],
            linestyle="none",
            marker=MARKERS[index % len(MARKERS)],
            markersize=4.2,
            markerfacecolor="#7A7F87",
            markeredgecolor="#7A7F87",
            label=definition,
        )
        for index, definition in enumerate(definitions)
    ]
    handles.append(Line2D([0], [0], color=INK, linewidth=1.2, label="median"))
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.20),
        ncol=min(4, len(handles)),
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.0,
        borderaxespad=0,
    )

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(prefix.with_suffix(".pdf"), dpi=dpi)
    fig.savefig(prefix.with_suffix(".png"), dpi=dpi)
    return fig
