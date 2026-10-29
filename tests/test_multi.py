from __future__ import annotations

import pandas as pd
import pytest

from clonedisp import analyze_multi_definition


def event_frame(state_values: dict[str, list[int]]) -> pd.DataFrame:
    rows = []
    for state, clone_counts in state_values.items():
        for clone, count in enumerate(clone_counts, start=1):
            rows.extend(
                {"sample": "P1", "state": state, "clone": f"C{clone}"}
                for _ in range(count)
            )
    return pd.DataFrame(rows)


def test_multi_definition_computes_each_definition_separately():
    definitions = {
        "res 0.5": event_frame(
            {
                "e1": [10, 10],
                "e2": [20, 0],
            }
        ),
        "res 1.0": event_frame(
            {
                "e1": [10, 0],
                "e2": [20, 0],
                "filler": [0, 10],
            }
        ),
    }
    result = analyze_multi_definition(
        definitions,
        min_clone_cells=5,
        min_state_cells=10,
        n_boot=20,
    )
    e1 = result.metrics.loc[result.metrics["state"].eq("e1")].set_index("definition")
    assert e1.loc["res 0.5", "clone_dispersion"] == pytest.approx(1.0)
    assert e1.loc["res 1.0", "clone_dispersion"] == pytest.approx(0.5)
    summary = result.state_summary.set_index("state")
    assert summary.loc["e1", "median_definition_cdi"] == pytest.approx(0.75)


def test_multi_definition_summary_detects_no_cross_definition_pooling():
    definitions = {
        "a": event_frame({"e1": [20, 0], "filler": [0, 10]}),
        "b": event_frame({"e1": [10, 10]}),
        "c": event_frame({"e1": [15, 5]}),
    }
    result = analyze_multi_definition(
        definitions,
        min_clone_cells=5,
        min_state_cells=10,
        n_boot=20,
    )
    assert result.metrics.loc[result.metrics["state"].eq("e1")].shape[0] == 3
    assert result.state_summary.iloc[0]["n_definitions"] == 3
    e1_summary = result.state_summary.set_index("state").loc["e1"]
    assert e1_summary["median_definition_cdi"] == pytest.approx(0.8)
    assert e1_summary["min_definition_cdi"] == pytest.approx(0.5)
    assert e1_summary["max_definition_cdi"] == pytest.approx(1.0)


def test_multi_definition_rejects_empty_mapping():
    with pytest.raises(ValueError, match="at least one"):
        analyze_multi_definition({})
