#!/usr/bin/env python3
"""Build an aggregate-only leaderboard across blind LOT model runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

from lot_data import assert_aggregate_only_report


FROZEN_FIELDS = (
    "folds",
    "retrieval_training_folds",
    "input_artifact_sha256",
    "run_purpose",
)
CONDITION_FIELDS = (
    "provider",
    "model",
    "reasoning_effort",
    "temperature",
    "retrieval_k",
    "prompt_version",
    "knowledge_version",
)


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _cohort(path: Path) -> tuple[str, int]:
    keys = sorted(
        json.loads(line)["case_key"]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    if len(keys) != len(set(keys)):
        raise ValueError(f"duplicate case key in {path}")
    encoded = json.dumps(keys, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest(), len(keys)


def _metric(report: dict[str, Any], *path: str) -> Any:
    value: Any = report
    for key in path:
        value = value[key]
    return value


def _load_run(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    metadata = _json(run_dir / "restricted" / "experiment_metadata.json")
    report = _json(run_dir / "public" / "aggregate_evaluation.json")
    cohort_hash, cohort_count = _cohort(
        run_dir / "restricted" / "joined_evaluation.jsonl"
    )
    return {
        "run_dir": str(run_dir),
        "metadata": metadata,
        "report": report,
        "cohort_sha256": cohort_hash,
        "cohort_count": cohort_count,
    }


def _condition_key(metadata: dict[str, Any]) -> tuple[Any, ...]:
    return tuple(metadata.get(field) for field in CONDITION_FIELDS)


def build_leaderboard(
    run_dirs: list[Path], *, phase: str, expected_repeats: int | None
) -> dict[str, Any]:
    if len(run_dirs) < 2:
        raise ValueError("provide at least two completed runs")
    runs = [_load_run(path) for path in run_dirs]
    reference = runs[0]
    expected_purpose = phase
    for run in runs:
        metadata = run["metadata"]
        if metadata.get("run_purpose") != expected_purpose:
            raise ValueError(f"{run['run_dir']} is not a {expected_purpose} run")
        mismatches = [
            field
            for field in FROZEN_FIELDS
            if metadata.get(field) != reference["metadata"].get(field)
        ]
        if run["cohort_sha256"] != reference["cohort_sha256"]:
            mismatches.append("evaluation_cohort")
        if mismatches:
            raise ValueError(
                f"{run['run_dir']} differs on frozen fields: {sorted(set(mismatches))}"
            )

    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for run in runs:
        grouped[_condition_key(run["metadata"])].append(run)

    rows = []
    for condition, condition_runs in grouped.items():
        repeats = [int(run["metadata"]["repeat_index"]) for run in condition_runs]
        if len(repeats) != len(set(repeats)):
            raise ValueError(f"duplicate repeat index for condition {condition}")
        if expected_repeats is not None and sorted(repeats) != list(
            range(1, expected_repeats + 1)
        ):
            raise ValueError(
                f"condition {condition} has repeats {sorted(repeats)}, "
                f"expected 1..{expected_repeats}"
            )
        exact = [
            _metric(
                run["report"], "metrics", "generated_total", "exact_accuracy", "rate"
            )
            for run in condition_runs
        ]
        mae = [
            _metric(run["report"], "metrics", "generated_total", "mean_absolute_error")
            for run in condition_runs
        ]
        within_one = [
            _metric(
                run["report"],
                "metrics",
                "generated_total",
                "within_one_accuracy",
                "rate",
            )
            for run in condition_runs
        ]
        usable = [
            _metric(run["report"], "metrics", "selective_ai", "coverage", "rate")
            for run in condition_runs
        ]
        three_way_coverage = [
            _metric(
                run["report"], "metrics", "three_way_consensus", "agreement", "rate"
            )
            for run in condition_runs
        ]
        three_way_accuracy = [
            _metric(
                run["report"],
                "metrics",
                "three_way_consensus",
                "reviewer_accuracy",
                "rate",
            )
            for run in condition_runs
        ]
        row = {
            **dict(zip(CONDITION_FIELDS, condition)),
            "repeat_count": len(condition_runs),
            "repeat_indices": sorted(repeats),
            "eligible_patients_per_repeat": reference["cohort_count"],
            "generated_exact_accuracy_mean": round(statistics.mean(exact), 6),
            "generated_exact_accuracy_sd": (
                round(statistics.stdev(exact), 6) if len(exact) > 1 else 0.0
            ),
            "generated_mae_mean": round(statistics.mean(mae), 6),
            "generated_within_one_mean": round(statistics.mean(within_one), 6),
            "usable_coverage_mean": round(statistics.mean(usable), 6),
            "three_way_coverage_mean": round(statistics.mean(three_way_coverage), 6),
            "three_way_accuracy_mean": (
                round(
                    statistics.mean(
                        value for value in three_way_accuracy if value is not None
                    ),
                    6,
                )
                if any(value is not None for value in three_way_accuracy)
                else None
            ),
            "run_ids": sorted(run["metadata"]["run_id"] for run in condition_runs),
        }
        rows.append(row)

    rows.sort(
        key=lambda row: (
            -row["generated_exact_accuracy_mean"],
            row["generated_mae_mean"],
            -row["generated_within_one_mean"],
            str(tuple(row.get(field) for field in CONDITION_FIELDS)),
        )
    )
    report = {
        "schema_version": "blind-lot-model-leaderboard-1.0.0",
        "report_scope": "aggregate_only",
        "phase": phase,
        "primary_selection_metric": "mean generated-total exact accuracy",
        "tie_breakers": [
            "lower mean absolute error",
            "higher mean within-one accuracy",
        ],
        "frozen_configuration": {
            field: reference["metadata"].get(field) for field in FROZEN_FIELDS
        },
        "evaluation_cohort_sha256": reference["cohort_sha256"],
        "evaluation_cohort_count": reference["cohort_count"],
        "conditions": rows,
        "selected_condition": rows[0] if phase != "final" else None,
        "final_phase_warning": (
            None
            if phase != "final"
            else "Final results are report-only; do not select or retune a condition on them."
        ),
    }
    assert_aggregate_only_report(report)
    return report


def write_leaderboard(report: dict[str, Any], output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=False)
    json_path = output_dir / "model_leaderboard.json"
    csv_path = output_dir / "model_leaderboard.csv"
    json_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        rows = report["conditions"]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return {"json": json_path, "csv": csv_path}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", action="append", required=True, type=Path)
    parser.add_argument(
        "--phase", choices=("screening", "development", "final"), required=True
    )
    parser.add_argument("--expected-repeats", type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.expected_repeats is not None and args.expected_repeats < 1:
        parser.error("--expected-repeats must be positive")
    report = build_leaderboard(
        args.run_dir,
        phase=args.phase,
        expected_repeats=args.expected_repeats,
    )
    outputs = write_leaderboard(report, args.output)
    best = report["conditions"][0]
    print(
        "Top condition: "
        f"{best['provider']} / {best['model']} / k={best['retrieval_k']} / "
        f"accuracy={best['generated_exact_accuracy_mean']:.1%}"
    )
    print(f"Leaderboard: {outputs['json']}")


if __name__ == "__main__":
    main()
