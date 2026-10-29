"""Core effective-clone and clone-dispersion calculations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Hashable

import numpy as np
import pandas as pd


ALGORITHM_NAME = "Effective Clone Dispersion (ECD)"
METRIC_NAME = "Clonal Dispersion Index (CDI)"


@dataclass(frozen=True)
class CloneDispResult:
    """Container for CloneDisp metrics and audit tables."""

    metrics: pd.DataFrame
    clone_validity: pd.DataFrame
    skipped: pd.DataFrame
    state_summary: pd.DataFrame
    parameters: dict[str, object]


def read_label_table(path: str | Path, sep: str | None = None) -> pd.DataFrame:
    """Read a CSV/TSV cell-label table while preserving label columns as text."""

    table_path = Path(path)
    if not table_path.exists():
        raise FileNotFoundError(f"Input table does not exist: {table_path}")
    if table_path.is_dir():
        raise IsADirectoryError(f"Expected an input file, got a directory: {table_path}")

    if sep is None:
        sep = "\t" if table_path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(table_path, sep=sep, dtype="string")


def _validate_inputs(
    data: pd.DataFrame,
    sample_col: Hashable,
    state_col: Hashable,
    clone_col: Hashable,
    cell_id_col: Hashable | None,
    min_clone_cells: int,
    min_state_cells: int,
    min_valid_clones: int,
) -> None:
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas.DataFrame")
    if data.empty:
        raise ValueError("data contains no cell rows")

    required = [sample_col, state_col, clone_col]
    if cell_id_col is not None:
        required.append(cell_id_col)
    missing = [str(column) for column in required if column not in data.columns]
    if missing:
        raise ValueError("Missing required column(s): " + ", ".join(missing))

    for column in required:
        if data[column].isna().any():
            count = int(data[column].isna().sum())
            raise ValueError(f"Column {column!r} contains {count} missing value(s)")

    for name, value in (
        ("min_clone_cells", min_clone_cells),
        ("min_state_cells", min_state_cells),
        ("min_valid_clones", min_valid_clones),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"{name} must be an integer")
        if int(value) < 1:
            raise ValueError(f"{name} must be at least 1")

    if cell_id_col is not None and data[cell_id_col].duplicated().any():
        examples = data.loc[data[cell_id_col].duplicated(), cell_id_col].head(5)
        examples_text = ", ".join(map(str, examples))
        raise ValueError(
            f"Column {cell_id_col!r} must uniquely identify cells; "
            f"duplicated examples: {examples_text}"
        )


def _colon_separated_labels(values: list[str]) -> str:
    """Return a stable, CSV-safe display for clone identifiers."""

    return ";".join(values)


def _bootstrap_mean_ci(
    values: np.ndarray,
    n_boot: int,
    seed: int,
) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, values.size, size=(n_boot, values.size))
    means = values[indices].mean(axis=1)
    low, median, high = np.quantile(means, [0.025, 0.5, 0.975])
    return float(low), float(median), float(high)


def summarize_by_state(
    metrics: pd.DataFrame,
    n_boot: int = 5000,
    seed: int = 0,
) -> pd.DataFrame:
    """Summarize patient-state metrics across patients with bootstrap CIs."""

    if n_boot < 1:
        raise ValueError("n_boot must be at least 1")

    columns = [
        "state",
        "n_patients",
        "total_cells",
        "mean_effective_clones",
        "mean_clone_dispersion",
        "median_clone_dispersion",
        "sd_clone_dispersion",
        "min_clone_dispersion",
        "max_clone_dispersion",
        "dispersion_ci_low",
        "dispersion_ci_high",
        "mean_effective_clone_ci_low",
        "mean_effective_clone_ci_high",
    ]
    if metrics.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, object]] = []
    states = sorted(metrics["state"].astype(str).unique())
    for state_index, state in enumerate(states):
        subset = metrics.loc[metrics["state"].astype(str).eq(state)]
        dispersion = subset["clone_dispersion"].to_numpy(dtype=float)
        effective = subset["effective_clones"].to_numpy(dtype=float)
        dispersion_low, _, dispersion_high = _bootstrap_mean_ci(
            dispersion,
            n_boot=n_boot,
            seed=seed + 2 * state_index,
        )
        effective_low, _, effective_high = _bootstrap_mean_ci(
            effective,
            n_boot=n_boot,
            seed=seed + 2 * state_index + 1,
        )
        rows.append(
            {
                "state": state,
                "n_patients": int(subset["sample"].nunique()),
                "total_cells": int(subset["n_cells_in_valid_clones"].sum()),
                "mean_effective_clones": float(np.mean(effective)),
                "mean_clone_dispersion": float(np.mean(dispersion)),
                "median_clone_dispersion": float(np.median(dispersion)),
                "sd_clone_dispersion": float(np.std(dispersion, ddof=1))
                if dispersion.size > 1
                else np.nan,
                "min_clone_dispersion": float(np.min(dispersion)),
                "max_clone_dispersion": float(np.max(dispersion)),
                "dispersion_ci_low": dispersion_low,
                "dispersion_ci_high": dispersion_high,
                "mean_effective_clone_ci_low": effective_low,
                "mean_effective_clone_ci_high": effective_high,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def analyze_clone_dispersion(
    data: pd.DataFrame,
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
) -> CloneDispResult:
    """Calculate effective clone number and Clonal Dispersion Index.

    Each row represents one cell. Clones below ``min_clone_cells`` are removed
    before state-level eligible cell counts are evaluated. A state also needs
    at least ``min_state_cells`` cells in the retained clones and at least
    ``min_valid_clones`` retained clones in its patient.
    """

    _validate_inputs(
        data,
        sample_col,
        state_col,
        clone_col,
        cell_id_col,
        min_clone_cells,
        min_state_cells,
        min_valid_clones,
    )

    labels = data[[sample_col, state_col, clone_col]].copy()
    labels.columns = ["sample", "state", "clone"]
    labels = labels.astype(str)

    clone_counts = (
        labels.groupby(["sample", "clone"], sort=True, observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    clone_counts["valid"] = clone_counts["n_cells"].ge(min_clone_cells)

    metric_rows: list[dict[str, object]] = []
    validity_rows: list[dict[str, object]] = []
    skipped_rows: list[dict[str, object]] = []

    clone_frame = labels[["sample", "clone"]].drop_duplicates()
    clone_frame = clone_frame.merge(
        clone_counts,
        on=["sample", "clone"],
        how="left",
        validate="one_to_one",
    )
    clone_frame["valid"] = clone_frame["n_cells"].ge(min_clone_cells)

    for sample in sorted(labels["sample"].unique()):
        sample_labels = labels.loc[labels["sample"].eq(sample)].copy()
        sample_clones = clone_frame.loc[clone_frame["sample"].eq(sample)].copy()
        valid_clones = sorted(
            sample_clones.loc[sample_clones["valid"], "clone"].astype(str).tolist()
        )
        excluded_clones = sorted(
            sample_clones.loc[~sample_clones["valid"], "clone"].astype(str).tolist()
        )
        valid_clone_set = set(valid_clones)
        valid_cells = int(
            sample_clones.loc[sample_clones["valid"], "n_cells"].sum()
        )
        observed_cells = int(sample_clones["n_cells"].sum())
        excluded_cells = observed_cells - valid_cells
        n_valid_clones = len(valid_clones)

        validity_rows.append(
            {
                "sample": sample,
                "n_cells_observed": observed_cells,
                "n_cells_in_valid_clones": valid_cells,
                "n_cells_excluded": excluded_cells,
                "valid_cell_fraction": valid_cells / observed_cells
                if observed_cells
                else np.nan,
                "n_clones_observed": len(sample_clones),
                "n_valid_clones": n_valid_clones,
                "valid_clone_fraction": n_valid_clones / len(sample_clones)
                if len(sample_clones)
                else np.nan,
                "valid_clone_ids": _colon_separated_labels(valid_clones),
                "excluded_clone_ids": _colon_separated_labels(excluded_clones),
            }
        )

        states = sorted(sample_labels["state"].unique())
        for state in states:
            state_labels = sample_labels.loc[sample_labels["state"].eq(state)]
            eligible_labels = state_labels.loc[
                state_labels["clone"].isin(valid_clone_set)
            ]
            n_cells_total = int(len(state_labels))
            n_cells_eligible = int(len(eligible_labels))

            reason = ""
            if n_valid_clones == 0:
                reason = "no_valid_clones"
            elif n_valid_clones < min_valid_clones:
                reason = "insufficient_valid_clones"
            elif n_cells_eligible < min_state_cells:
                reason = "insufficient_state_cells"

            if reason:
                skipped_rows.append(
                    {
                        "sample": sample,
                        "state": state,
                        "n_cells_total": n_cells_total,
                        "n_cells_in_valid_clones": n_cells_eligible,
                        "n_valid_clones": n_valid_clones,
                        "reason": reason,
                    }
                )
                continue

            clone_state_counts = (
                eligible_labels["clone"]
                .value_counts()
                .reindex(valid_clones, fill_value=0)
                .astype(float)
            )
            distribution = clone_state_counts.to_numpy(dtype=float, copy=True)
            distribution /= distribution.sum()
            effective_clones = 1.0 / float(np.sum(distribution**2))
            dispersion = effective_clones / n_valid_clones
            dominant_index = int(np.argmax(distribution))

            metric_rows.append(
                {
                    "sample": sample,
                    "state": state,
                    "n_cells_total": n_cells_total,
                    "n_cells_in_valid_clones": n_cells_eligible,
                    "n_valid_clones": n_valid_clones,
                    "effective_clones": effective_clones,
                    "clone_dispersion": dispersion,
                    "theoretical_minimum": 1.0 / n_valid_clones,
                    "max_clone_fraction": float(np.max(distribution)),
                    "dominant_clone": valid_clones[dominant_index],
                }
            )

    metrics = pd.DataFrame(
        metric_rows,
        columns=[
            "sample",
            "state",
            "n_cells_total",
            "n_cells_in_valid_clones",
            "n_valid_clones",
            "effective_clones",
            "clone_dispersion",
            "theoretical_minimum",
            "max_clone_fraction",
            "dominant_clone",
        ],
    )
    clone_validity = pd.DataFrame(
        validity_rows,
        columns=[
            "sample",
            "n_cells_observed",
            "n_cells_in_valid_clones",
            "n_cells_excluded",
            "valid_cell_fraction",
            "n_clones_observed",
            "n_valid_clones",
            "valid_clone_fraction",
            "valid_clone_ids",
            "excluded_clone_ids",
        ],
    )
    skipped = pd.DataFrame(
        skipped_rows,
        columns=[
            "sample",
            "state",
            "n_cells_total",
            "n_cells_in_valid_clones",
            "n_valid_clones",
            "reason",
        ],
    )
    state_summary = summarize_by_state(metrics, n_boot=n_boot, seed=seed)
    parameters: dict[str, object] = {
        "algorithm": ALGORITHM_NAME,
        "metric": METRIC_NAME,
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
    return CloneDispResult(
        metrics=metrics,
        clone_validity=clone_validity,
        skipped=skipped,
        state_summary=state_summary,
        parameters=parameters,
    )
