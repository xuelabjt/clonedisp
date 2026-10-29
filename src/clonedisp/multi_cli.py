"""CLI for multi-definition CloneDisp analysis."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import read_label_table
from .multi import analyze_multi_definition


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="clonedisp-multi",
        description=(
            "Calculate effective clone number and CDI independently across "
            "multiple clone definitions, then summarize their sensitivity."
        ),
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--sample", default="sample")
    parser.add_argument("--state", default="state")
    parser.add_argument("--clone", default="clone")
    parser.add_argument("--cell-id", default=None)
    parser.add_argument(
        "--sample-regex",
        default=None,
        help="Derive sample from --cell-id using the first capture group.",
    )
    parser.add_argument(
        "--state-replace",
        default=None,
        metavar="OLD=NEW",
        help="Replace a state-label prefix, for example newepi=e.",
    )
    parser.add_argument("--min-clone-cells", type=int, default=5)
    parser.add_argument("--min-state-cells", type=int, default=20)
    parser.add_argument("--min-valid-clones", type=int, default=1)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument(
        "--state-order",
        default=None,
        help="Comma-separated state order for plotting.",
    )
    parser.add_argument(
        "--focus",
        action="append",
        default=[],
        help="State to highlight in the plot; repeatable.",
    )
    parser.add_argument("--title", default="Clone dispersion across definitions")
    parser.add_argument("--fig-width", type=float, default=4.3)
    parser.add_argument("--fig-height", type=float, default=3.1)
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def _read_manifest(path: Path) -> list[tuple[str, Path]]:
    separator = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    manifest = pd.read_csv(path, sep=separator, dtype="string")
    required = {"definition", "path"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(
            "Manifest requires columns: definition,path; missing "
            + ", ".join(sorted(missing))
        )
    if manifest.isna().any(axis=None):
        raise ValueError("Manifest contains missing definition or path values")
    entries = []
    for row in manifest.itertuples(index=False):
        definition = str(row.definition)
        table_path = Path(str(row.path)).expanduser()
        if not table_path.is_absolute():
            table_path = path.parent / table_path
        entries.append((definition, table_path))
    if not entries:
        raise ValueError("Manifest contains no definition rows")
    return entries


def _prepare_table(
    table: pd.DataFrame,
    *,
    sample_col: str,
    state_col: str,
    cell_id_col: str | None,
    sample_regex: str | None,
    state_replace: str | None,
) -> pd.DataFrame:
    table = table.copy()
    if sample_col not in table.columns:
        if cell_id_col is None or sample_regex is None:
            raise ValueError(
                f"Column {sample_col!r} is absent; provide --cell-id and --sample-regex"
            )
        if cell_id_col not in table.columns:
            raise ValueError(f"Cell ID column not found: {cell_id_col}")
        extracted = table[cell_id_col].astype(str).str.extract(
            sample_regex,
            expand=False,
        )
        if isinstance(extracted, pd.DataFrame):
            raise ValueError("--sample-regex must contain exactly one capture group")
        if extracted.isna().any():
            raise ValueError("--sample-regex did not match every cell identifier")
        table[sample_col] = extracted

    if state_replace:
        if "=" not in state_replace:
            raise ValueError("--state-replace must use OLD=NEW")
        old, new = state_replace.split("=", 1)
        table[state_col] = table[state_col].astype(str).str.replace(
            old,
            new,
            regex=False,
        )
    return table


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        entries = _read_manifest(args.manifest.expanduser().resolve())
        definitions = {}
        for definition, table_path in entries:
            if definition in definitions:
                raise ValueError(f"Duplicate definition in manifest: {definition}")
            table = read_label_table(table_path)
            definitions[definition] = _prepare_table(
                table,
                sample_col=args.sample,
                state_col=args.state,
                cell_id_col=args.cell_id,
                sample_regex=args.sample_regex,
                state_replace=args.state_replace,
            )
        result = analyze_multi_definition(
            definitions,
            sample_col=args.sample,
            state_col=args.state,
            clone_col=args.clone,
            cell_id_col=args.cell_id,
            min_clone_cells=args.min_clone_cells,
            min_state_cells=args.min_state_cells,
            min_valid_clones=args.min_valid_clones,
            n_boot=args.bootstrap,
            seed=args.seed,
        )
    except (OSError, TypeError, ValueError) as exc:
        parser.exit(2, f"clonedisp-multi: error: {exc}\n")

    prefix = args.output_prefix.expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    output_paths = {
        "definition_patient_metrics": Path(f"{prefix}_definition_patient_metrics.csv"),
        "definition_summary": Path(f"{prefix}_definition_summary.csv"),
        "state_summary": Path(f"{prefix}_state_summary.csv"),
        "skipped": Path(f"{prefix}_skipped.csv"),
        "run_metadata": Path(f"{prefix}_run.json"),
    }
    result.metrics.to_csv(output_paths["definition_patient_metrics"], index=False)
    result.definition_summary.to_csv(output_paths["definition_summary"], index=False)
    result.state_summary.to_csv(output_paths["state_summary"], index=False)
    result.skipped.to_csv(output_paths["skipped"], index=False)

    if args.plot:
        try:
            from .plotting import plot_multi_definition
        except ImportError as exc:
            parser.exit(
                2,
                "clonedisp-multi: error: plotting requires matplotlib; "
                "install with `pip install .[plot]`\n",
            )
        order = (
            [state.strip() for state in args.state_order.split(",")]
            if args.state_order
            else None
        )
        plot_multi_definition(
            result,
            prefix,
            state_order=order,
            focus_states=args.focus,
            title=args.title,
            figsize=(args.fig_width, args.fig_height),
        )

    metadata = {
        "software": "CloneDisp",
        "version": __version__,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(args.manifest.expanduser().resolve()),
        "parameters": result.parameters,
        "counts": {
            "definition_patient_metrics": int(len(result.metrics)),
            "definition_summary": int(len(result.definition_summary)),
            "state_summary": int(len(result.state_summary)),
            "skipped": int(len(result.skipped)),
        },
        "outputs": {name: str(path) for name, path in output_paths.items()},
    }
    with output_paths["run_metadata"].open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(
        f"CloneDisp multi: analyzed {len(definitions)} definitions and "
        f"{result.state_summary['state'].nunique()} states"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
