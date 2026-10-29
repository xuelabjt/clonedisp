"""Multi-definition effective-clone analysis."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Hashable

import numpy as np
import pandas as pd

from .core import (
    ALGORITHM_NAME,
    METRIC_NAME,
    CloneDispResult,
    analyze_clone_dispersion,
)


@dataclass(frozen=True)
class MultiCloneDispResult:
    """Results from independent clone definitions."""

    metrics: pd.DataFrame
    definition_summary: pd.DataFrame
    state_summary: pd.DataFrame
    skipped: pd.DataFrame
    definition_results: dict[str, CloneDispResult]
    parameters: dict[str, object]


def _definition_state_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "definition",
        "state",
        "n_patients",
        "sum_definition_cells",
        "mean_effective_clones",
        "mean_clone_dispersion",
        "median_clone_dispersion",
    ]
    if metrics.empty:
        return pd.DataFrame(columns=columns)
    summary = (
        metrics.groupby(["definition", "state"], sort=True, observed=True)
        .agg(
            n_patients=("sample", "nunique"),
            sum_definition_cells=("n_cells_in_valid_clones", "sum"),
            mean_effective_clones=("effective_clones", "mean"),
            mean_clone_dispersion=("clone_dispersion", "mean"),
            median_clone_dispersion=("clone_dispersion", "median"),
        )
        .reset_index()
    )
    return summary[columns]


def _cross_definition_summary(
    definition_summary: pd.DataFrame,
) -> pd.DataFrame:
    columns = [
        "state",
        "n_definitions",
        "sum_definition_cells",
        "mean_definition_cdi",
        "median_definition_cdi",
        "sd_definition_cdi",
        "min_definition_cdi",
        "max_definition_cdi",
        "mean_definition_effective_clones",
        "median_definition_effective_clones",
        "rank_median_cdi",
    ]
    if definition_summary.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, object]] = []
    for state, subset in definition_summary.groupby("state", sort=True):
        cdi = subset["mean_clone_dispersion"].to_numpy(dtype=float)
        effective = subset["mean_effective_clones"].to_numpy(dtype=float)
        rows.append(
            {
                "state": str(state),
                "n_definitions": int(len(subset)),
                "sum_definition_cells": int(subset["sum_definition_cells"].sum()),
                "mean_definition_cdi": float(np.mean(cdi)),
                "median_definition_cdi": float(np.median(cdi)),
                "sd_definition_cdi": float(np.std(cdi, ddof=1))
                if cdi.size > 1
                else np.nan,
                "min_definition_cdi": float(np.min(cdi)),
                "max_definition_cdi": float(np.max(cdi)),
                "mean_definition_effective_clones": float(np.mean(effective)),
                "median_definition_effective_clones": float(np.median(effective)),
            }
        )

    summary = pd.DataFrame(rows)
    summary["rank_median_cdi"] = (
        summary["median_definition_cdi"]
        .rank(ascending=False, method="min")
        .astype(int)
    )
    return summary[columns]


def analyze_multi_definition(
    definitions: Mapping[str, pd.DataFrame],
    *,
    sample_col: Hashable = "sample",
    state_col: Hashable = "state",
    clone_col: Hashable = "clone",
    cell_id_col: Hashable | None = None,
    min_clone_cells: int = 5,
    min_state_cells: int = 20,
    min_valid_clones: int = 1,
    n_boot: int = 5000,
    seed: int = 0,
) -> MultiCloneDispResult:
    """Calculate CloneDisp independently for each clone definition.

    Definitions are never pooled at the cell-label level. The returned state
    summary contains one patient-level mean CDI per definition and summarizes
    those values across definitions.
    """

    if not definitions:
        raise ValueError("definitions must contain at least one clone definition")

    definition_results: dict[str, CloneDispResult] = {}
    metric_tables: list[pd.DataFrame] = []
    skipped_tables: list[pd.DataFrame] = []

    for definition_index, (definition, frame) in enumerate(definitions.items()):
        definition = str(definition).strip()
        if not definition:
            raise ValueError("definition names must be non-empty")
        if definition in definition_results:
            raise ValueError(f"duplicate definition name: {definition}")

        result = analyze_clone_dispersion(
            frame,
            sample_col=sample_col,
            state_col=state_col,
            clone_col=clone_col,
            cell_id_col=cell_id_col,
            min_clone_cells=min_clone_cells,
            min_state_cells=min_state_cells,
            min_valid_clones=min_valid_clones,
            n_boot=n_boot,
            seed=seed + definition_index,
        )
        definition_results[definition] = result

        metrics = result.metrics.copy()
        if not metrics.empty:
            metrics.insert(0, "definition", definition)
            metric_tables.append(metrics)

        skipped = result.skipped.copy()
        if not skipped.empty:
            skipped.insert(0, "definition", definition)
            skipped_tables.append(skipped)

    metrics = (
        pd.concat(metric_tables, ignore_index=True)
        if metric_tables
        else definition_results[next(iter(definition_results))].metrics.copy()
    )
    if metrics.empty and "definition" not in metrics.columns:
        metrics.insert(0, "definition", pd.Series(dtype="string"))

    skipped = (
        pd.concat(skipped_tables, ignore_index=True)
        if skipped_tables
        else pd.DataFrame(
            columns=[
                "definition",
                "sample",
                "state",
                "n_cells_total",
                "n_cells_in_valid_clones",
                "n_valid_clones",
                "reason",
            ]
        )
    )
    definition_summary = _definition_state_summary(metrics)
    state_summary = _cross_definition_summary(definition_summary)
    parameters: dict[str, object] = {
        "algorithm": ALGORITHM_NAME,
        "metric": METRIC_NAME,
        "analysis_mode": "multi_definition",
        "definitions": list(definition_results),
        "sample_col": str(sample_col),
        "state_col": str(state_col),
        "clone_col": str(clone_col),
        "cell_id_col": str(cell_id_col) if cell_id_col is not None else None,
        "min_clone_cells": int(min_clone_cells),
        "min_state_cells": int(min_state_cells),
        "min_valid_clones": int(min_valid_clones),
        "bootstrap_replicates": int(n_boot),
        "seed": int(seed),
    }
    return MultiCloneDispResult(
        metrics=metrics,
        definition_summary=definition_summary,
        state_summary=state_summary,
        skipped=skipped,
        definition_results=definition_results,
        parameters=parameters,
    )
