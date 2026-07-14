#!/usr/bin/env python3
"""Compare policy analyses from completed blind LOT runs without API calls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from blind_lot.reaggregation import REAGGREGATED_AGGREGATE_NAME, REAGGREGATED_JOINED_NAME


FROZEN_METADATA_FIELDS = (
    "folds", "model", "reasoning_effort", "prompt_version", "knowledge_version",
)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _cohort_hash(joined_path: Path) -> tuple[str, int]:
    case_keys = sorted(
        json.loads(line)["case_key"]
        for line in joined_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    payload = json.dumps(case_keys, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest(), len(case_keys)


def load_runs(run_dirs: list[Path]) -> list[dict[str, Any]]:
    if len(run_dirs) < 2:
        raise ValueError("provide at least two --run-dir values")
    runs = []
    for run_dir in run_dirs:
        run_dir = run_dir.resolve()
        metadata_path = run_dir / "restricted" / "experiment_metadata.json"
        joined = run_dir / "restricted" / REAGGREGATED_JOINED_NAME
        if not joined.exists():
            joined = run_dir / "restricted" / "joined_evaluation.jsonl"
        aggregate = run_dir / "public" / REAGGREGATED_AGGREGATE_NAME
        if not aggregate.exists():
            aggregate = run_dir / "public" / "aggregate_evaluation.json"
        for path in (metadata_path, joined, aggregate):
            if not path.exists():
                raise FileNotFoundError(path)
        metadata = _load_json(metadata_path)
        report = _load_json(aggregate)
        if "routing_policy_analysis" not in report.get("metrics", {}):
            raise ValueError(f"{run_dir} must be evaluation-only reaggregated first")
        cohort_sha256, cohort_count = _cohort_hash(joined)
        runs.append({
            "run_dir": run_dir,
            "metadata": metadata,
            "report": report,
            "cohort_sha256": cohort_sha256,
            "cohort_count": cohort_count,
        })

    reference = runs[0]
    for run in runs[1:]:
        mismatches = [
            field for field in FROZEN_METADATA_FIELDS
            if run["metadata"].get(field) != reference["metadata"].get(field)
        ]
        if run["metadata"].get("input_artifact_sha256") != reference["metadata"].get("input_artifact_sha256"):
            mismatches.append("input_artifact_sha256")
        if run["cohort_sha256"] != reference["cohort_sha256"]:
            mismatches.append("evaluation_cohort")
        if (
            run["report"]["experiment"].get("evaluation_version")
            != reference["report"]["experiment"].get("evaluation_version")
        ):
            mismatches.append("evaluation_version")
        if mismatches:
            raise ValueError(
                f"{run['run_dir']} does not share frozen comparison fields: {mismatches}"
            )
    retrieval_values = [run["metadata"].get("retrieval_k") for run in runs]
    if len(set(retrieval_values)) != len(retrieval_values):
        raise ValueError("compared runs must have distinct retrieval k values")
    return sorted(runs, key=lambda item: item["metadata"]["retrieval_k"])


def _metric_rate(policy: dict[str, Any], name: str) -> float | None:
    return policy[name]["rate"]


def comparison_rows(runs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for run in runs:
        report = run["report"]
        metrics = report["metrics"]
        policies = metrics["routing_policy_analysis"]
        row: dict[str, Any] = {
            "run_dir": str(run["run_dir"]),
            "retrieval_k": run["metadata"]["retrieval_k"],
            "eligible_patients": metrics["eligible_patients"],
            "generated_total_accuracy": metrics["generated_total"]["exact_accuracy"]["rate"],
            "selective_ai_coverage": metrics["selective_ai"]["coverage"]["rate"],
            "selective_ai_accuracy": metrics["selective_ai"]["exact_accuracy"]["rate"],
            "algorithm_only_accuracy": _metric_rate(
                policies["algorithm_only"], "reviewer_accuracy_among_accepted"
            ),
            "cota_only_accuracy": _metric_rate(
                policies["cota_only"], "reviewer_accuracy_among_accepted"
            ),
        }
        for name in (
            "algorithm_cota_agreement", "algorithm_ai_agreement",
            "cota_ai_agreement", "three_way_agreement",
        ):
            policy = policies[name]
            prefix = name.removesuffix("_agreement")
            row[f"{prefix}_accepted"] = policy["accepted_count"]
            row[f"{prefix}_coverage"] = policy["coverage"]["rate"]
            row[f"{prefix}_correct"] = policy["correct_accepted_count"]
            row[f"{prefix}_incorrect"] = policy["incorrect_accepted_count"]
            row[f"{prefix}_accuracy"] = policy["reviewer_accuracy_among_accepted"]["rate"]
        rows.append(row)
    return rows


def build_comparison(runs: list[dict[str, Any]]) -> dict[str, Any]:
    reference = runs[0]
    return {
        "schema_version": "1.0.0",
        "report_scope": "aggregate_only",
        "comparison": "order-only blind AI routing policies by retrieval k",
        "primary_policy": "three_way_agreement",
        "secondary_policies_are_exploratory": True,
        "frozen_configuration": {
            **{field: reference["metadata"].get(field) for field in FROZEN_METADATA_FIELDS},
            "input_artifact_sha256": reference["metadata"].get("input_artifact_sha256"),
            "evaluation_cohort_sha256": reference["cohort_sha256"],
            "evaluation_cohort_count": reference["cohort_count"],
            "evaluation_version": reference["report"]["experiment"].get("evaluation_version"),
        },
        "rows": comparison_rows(runs),
        "policy_analyses": [
            {
                "retrieval_k": run["metadata"]["retrieval_k"],
                "routing_policy_analysis": run["report"]["metrics"]["routing_policy_analysis"],
                "conditional_agreement_analysis": run["report"]["metrics"]["conditional_agreement_analysis"],
            }
            for run in runs
        ],
        "small_subset_warning": (
            "High precision from very small accepted subsets must not be overinterpreted."
        ),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_comparison(runs: list[dict[str, Any]], output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    result = build_comparison(runs)
    json_path = output_dir / "policy_comparison.json"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    policy_csv = output_dir / "policy_comparison.csv"
    _write_csv(policy_csv, result["rows"])
    pattern_rows = []
    for run in runs:
        patterns = run["report"]["metrics"]["vote_pattern_analysis"]["patterns"]
        for pattern, values in patterns.items():
            pattern_rows.append({
                "retrieval_k": run["metadata"]["retrieval_k"],
                "vote_pattern": pattern,
                "patient_count": values["patient_count"],
                "proportion": values["proportion_of_eligible_patients"]["rate"],
                "algorithm_accuracy": values["reviewer_accuracy"]["algorithm"]["rate"],
                "cota_accuracy": values["reviewer_accuracy"]["cota"]["rate"],
                "usable_ai_accuracy": values["reviewer_accuracy"]["ai_when_usable"]["rate"],
                "algorithm_uniquely_correct": values["winner_attribution"]["algorithm_uniquely_correct"],
                "cota_uniquely_correct": values["winner_attribution"]["cota_uniquely_correct"],
                "ai_uniquely_correct": values["winner_attribution"]["ai_uniquely_correct"],
                "none_correct": values["winner_attribution"]["none_correct"],
            })
    pattern_csv = output_dir / "vote_pattern_comparison.csv"
    _write_csv(pattern_csv, pattern_rows)
    return {"json": json_path, "policy_csv": policy_csv, "vote_pattern_csv": pattern_csv}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    runs = load_runs(args.run_dir)
    outputs = write_comparison(runs, args.output)
    for row in comparison_rows(runs):
        print(
            f"k={row['retrieval_k']}: generated={row['generated_total_accuracy']}, "
            f"usable coverage={row['selective_ai_coverage']}, "
            f"three-way accepted={row['three_way_accepted']}"
        )
    print(f"Comparison outputs: {outputs['json'].parent}")


if __name__ == "__main__":
    main()
