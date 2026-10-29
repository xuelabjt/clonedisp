"""Command-line interface for CloneDisp."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .core import ALGORITHM_NAME, METRIC_NAME, analyze_clone_dispersion, read_label_table


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="clonedisp",
        description=(
            "Calculate the inverse-Simpson effective clone number and the "
            "normalized Clonal Dispersion Index from cell-level labels."
        ),
    )
    parser.add_argument("--input", required=True, type=Path, help="Input CSV/TSV")
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--sample", default="sample", help="Sample column name")
    parser.add_argument("--state", default="state", help="State column name")
    parser.add_argument("--clone", default="clone", help="Clone column name")
    parser.add_argument(
        "--cell-id",
        default=None,
        help="Optional unique cell identifier column",
    )
    parser.add_argument("--min-clone-cells", type=int, default=5)
    parser.add_argument("--min-state-cells", type=int, default=20)
    parser.add_argument("--min-valid-clones", type=int, default=1)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--sep",
        default=None,
        help="Input delimiter. Default: comma for CSV, tab for TSV/TXT.",
    )
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    output_prefix = args.output_prefix.expanduser().resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        labels = read_label_table(args.input, sep=args.sep)
        result = analyze_clone_dispersion(
            labels,
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
        parser.exit(2, f"clonedisp: error: {exc}\n")

    outputs = {
        "metrics": Path(f"{output_prefix}_metrics.csv"),
        "clone_validity": Path(f"{output_prefix}_clone_validity.csv"),
        "state_summary": Path(f"{output_prefix}_state_summary.csv"),
        "skipped": Path(f"{output_prefix}_skipped.csv"),
        "run_metadata": Path(f"{output_prefix}_run.json"),
    }
    result.metrics.to_csv(outputs["metrics"], index=False)
    result.clone_validity.to_csv(outputs["clone_validity"], index=False)
    result.state_summary.to_csv(outputs["state_summary"], index=False)
    result.skipped.to_csv(outputs["skipped"], index=False)

    metadata = {
        "algorithm": ALGORITHM_NAME,
        "metric": METRIC_NAME,
        "software": "CloneDisp",
        "version": __version__,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(args.input.expanduser().resolve()),
        "input_sha256": _sha256(args.input),
        "formula": {
            "effective_clones": "1 / sum(q_c^2)",
            "clone_dispersion": "effective_clones / n_valid_clones",
        },
        "parameters": result.parameters,
        "counts": {
            "metric_rows": int(len(result.metrics)),
            "clone_validity_rows": int(len(result.clone_validity)),
            "state_summary_rows": int(len(result.state_summary)),
            "skipped_rows": int(len(result.skipped)),
        },
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    with outputs["run_metadata"].open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    if not args.quiet:
        print(
            f"CloneDisp: wrote {len(result.metrics)} patient-state metric rows "
            f"for {result.metrics['state'].nunique() if not result.metrics.empty else 0} "
            f"states to {output_prefix}_*.csv"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
