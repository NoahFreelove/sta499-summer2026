#!/usr/bin/env python3
"""Validate coding-agent shards and finalize a blind LOT benchmark run."""

from __future__ import annotations

import argparse
from pathlib import Path

from blind_lot.agent_handoff import (
    finalize_agent_run,
    import_worker_bundle_outputs,
)
from lot_data import find_project_root
from run_blind_lot_benchmark import format_metric


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--project-root", type=Path)
    parser.add_argument(
        "--worker-bundle-dir",
        type=Path,
        help="verify and import worker outputs before finalization",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.worker_bundle_dir is not None:
        import_worker_bundle_outputs(args.run_dir, args.worker_bundle_dir)
    result = finalize_agent_run(args.run_dir, find_project_root(args.project_root))
    metrics = result["public"]["metrics"]
    print(
        f"Evaluated {metrics['eligible_patients']} patients against reviewer consensus."
    )
    print(
        format_metric(
            "Generated-total exact accuracy",
            metrics["generated_total"]["exact_accuracy"],
        )
    )
    print(format_metric("Usable AI vote coverage", metrics["selective_ai"]["coverage"]))
    print(
        format_metric(
            "Non-abstained exact accuracy", metrics["selective_ai"]["exact_accuracy"]
        )
    )
    print(
        format_metric(
            "Three-way consensus coverage", metrics["three_way_consensus"]["agreement"]
        )
    )
    print("Python API provider calls: 0")
    print(f"Aggregate output: {result['aggregate_output']}")


if __name__ == "__main__":
    main()
