from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from clonedisp import analyze_clone_dispersion


def make_labels(counts: dict[tuple[str, str, str], int]) -> pd.DataFrame:
    rows = []
    for (sample, state, clone), count in counts.items():
        rows.extend(
            {
                "sample": sample,
                "state": state,
                "clone": clone,
            }
            for _ in range(count)
        )
    return pd.DataFrame(rows)


def metric_row(result, sample: str, state: str):
    query = result.metrics.query("sample == @sample and state == @state")
    assert len(query) == 1
    return query.iloc[0]


def test_state_restricted_to_one_of_four_clones():
    labels = make_labels(
        {
            ("P1", "restricted", "C1"): 10,
            ("P1", "restricted", "C2"): 0,
            ("P1", "restricted", "C3"): 0,
            ("P1", "restricted", "C4"): 0,
            ("P1", "other", "C1"): 1,
            ("P1", "other", "C2"): 10,
            ("P1", "other", "C3"): 10,
            ("P1", "other", "C4"): 10,
        }
    )
    result = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=10,
        n_boot=20,
    )
    row = metric_row(result, "P1", "restricted")
    assert row["n_valid_clones"] == 4
    assert row["effective_clones"] == pytest.approx(1.0)
    assert row["clone_dispersion"] == pytest.approx(0.25)
    assert row["theoretical_minimum"] == pytest.approx(0.25)
    assert row["dominant_clone"] == "C1"


def test_state_evenly_distributed_across_four_clones():
    labels = make_labels(
        {
            ("P1", "uniform", "C1"): 5,
            ("P1", "uniform", "C2"): 5,
            ("P1", "uniform", "C3"): 5,
            ("P1", "uniform", "C4"): 5,
        }
    )
    result = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=20,
        n_boot=20,
    )
    row = metric_row(result, "P1", "uniform")
    assert row["effective_clones"] == pytest.approx(4.0)
    assert row["clone_dispersion"] == pytest.approx(1.0)
    assert row["max_clone_fraction"] == pytest.approx(0.25)


def test_unequal_distribution_matches_manual_calculation():
    labels = make_labels(
        {
            ("P1", "mixed", "C1"): 6,
            ("P1", "mixed", "C2"): 3,
            ("P1", "mixed", "C3"): 1,
            ("P1", "filler", "C1"): 1,
            ("P1", "filler", "C2"): 4,
            ("P1", "filler", "C3"): 6,
        }
    )
    result = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=5,
        n_boot=20,
    )
    row = metric_row(result, "P1", "mixed")
    expected_effective = 1.0 / (0.6**2 + 0.3**2 + 0.1**2)
    assert row["effective_clones"] == pytest.approx(expected_effective)
    assert row["clone_dispersion"] == pytest.approx(expected_effective / 3)
    assert row["max_clone_fraction"] == pytest.approx(0.6)


def test_small_clone_is_excluded_before_state_filtering():
    labels = make_labels(
        {
            ("P1", "target", "C1"): 8,
            ("P1", "target", "C2"): 5,
            ("P1", "target", "C3_tiny"): 4,
        }
    )
    result = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=5,
        n_boot=20,
    )
    row = metric_row(result, "P1", "target")
    validity = result.clone_validity.set_index("sample").loc["P1"]
    assert row["n_cells_total"] == 17
    assert row["n_cells_in_valid_clones"] == 13
    assert row["n_valid_clones"] == 2
    assert validity["n_cells_excluded"] == 4
    assert validity["valid_clone_ids"] == "C1;C2"
    assert validity["excluded_clone_ids"] == "C3_tiny"


def test_state_below_minimum_is_reported_as_skipped():
    labels = make_labels(
        {
            ("P1", "sparse", "C1"): 10,
            ("P1", "sparse", "C2"): 9,
        }
    )
    result = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=20,
        n_boot=20,
    )
    assert result.metrics.empty
    skipped = result.skipped.iloc[0]
    assert skipped["reason"] == "insufficient_state_cells"
    assert skipped["n_cells_in_valid_clones"] == 19


def test_patient_without_valid_clone_is_audited():
    labels = make_labels(
        {
            ("P1", "s0", "C1"): 4,
            ("P1", "s0", "C2"): 4,
        }
    )
    result = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=1,
        n_boot=20,
    )
    assert result.skipped.iloc[0]["reason"] == "no_valid_clones"
    assert result.clone_validity.iloc[0]["n_valid_clones"] == 0


def test_min_valid_clones_can_enforce_comparability():
    labels = make_labels(
        {
            ("P1", "s0", "C1"): 20,
            ("P1", "s0", "C2"): 1,
        }
    )
    result = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=5,
        min_valid_clones=2,
        n_boot=20,
    )
    assert result.metrics.empty
    assert result.skipped.iloc[0]["reason"] == "insufficient_valid_clones"


def test_state_summary_deduplicates_patients_and_is_deterministic():
    labels = make_labels(
        {
            ("P1", "s0", "C1"): 5,
            ("P1", "s0", "C2"): 5,
            ("P2", "s0", "C1"): 10,
            ("P2", "s0", "C2"): 10,
            ("P3", "s0", "C1"): 20,
            ("P3", "s0", "C2"): 5,
        }
    )
    first = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=5,
        n_boot=50,
        seed=7,
    )
    second = analyze_clone_dispersion(
        labels,
        min_clone_cells=5,
        min_state_cells=5,
        n_boot=50,
        seed=7,
    )
    row = first.state_summary.iloc[0]
    p3_effective = 1.0 / (0.8**2 + 0.2**2)
    expected_mean = np.mean([1.0, 1.0, p3_effective / 2.0])
    assert row["n_patients"] == 3
    assert row["total_cells"] == 55
    assert row["mean_clone_dispersion"] == pytest.approx(expected_mean)
    pd.testing.assert_frame_equal(
        first.state_summary,
        second.state_summary,
        check_exact=True,
    )


def test_missing_required_column_raises_clear_error():
    labels = pd.DataFrame({"sample": ["P1"], "clone": ["C1"]})
    with pytest.raises(ValueError, match="state"):
        analyze_clone_dispersion(labels)


def test_missing_label_value_is_rejected():
    labels = pd.DataFrame(
        {
            "sample": ["P1", "P1"],
            "state": ["s0", None],
            "clone": ["C1", "C1"],
        }
    )
    with pytest.raises(ValueError, match="missing value"):
        analyze_clone_dispersion(labels)


def test_duplicate_cell_ids_are_rejected_when_supplied():
    labels = pd.DataFrame(
        {
            "sample": ["P1", "P1"],
            "state": ["s0", "s0"],
            "clone": ["C1", "C1"],
            "cell_id": ["cell-1", "cell-1"],
        }
    )
    with pytest.raises(ValueError, match="uniquely identify"):
        analyze_clone_dispersion(labels, cell_id_col="cell_id")
