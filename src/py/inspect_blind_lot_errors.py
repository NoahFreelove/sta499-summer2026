#!/usr/bin/env python3
"""Write restricted diagnostics from saved blind LOT artifacts only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from blind_lot.evaluation import POLICY_SPECS, normalize_joined_record
from lot_data import read_jsonl


def _source_dirs(run_dir: Path) -> list[Path]:
    direct = run_dir / "restricted" / "blind_predictions.jsonl"
    if direct.exists():
        return [run_dir]
    source_path = run_dir / "restricted" / "pooled_source_runs.json"
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    sources = json.loads(source_path.read_text(encoding="utf-8"))["sources"]
    benchmarks = next((parent for parent in run_dir.resolve().parents if parent.name == "benchmarks"), None)
    if benchmarks is None:
        raise ValueError("cannot locate artifacts/benchmarks ancestor for pooled run")
    return [benchmarks / "runs" / source["run_id"] for source in sources]


def _index(path: Path) -> dict[str, dict[str, Any]]:
    return {row["case_key"]: row for row in read_jsonl(path)}


def _details(run_dir: Path) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for source in _source_dirs(run_dir):
        restricted = source / "restricted"
        predictions = _index(restricted / "blind_predictions.jsonl")
        retrieval = _index(restricted / "retrieval_debug.jsonl")
        for key in predictions:
            output[key] = {"prediction": predictions[key], "retrieval": retrieval.get(key, {})}
    return output


def _trajectories(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return {row["case_key"]: row.get("trajectory") for row in read_jsonl(path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--compare-run-dir", type=Path)
    parser.add_argument("--policy", default="three_way_agreement", choices=sorted(POLICY_SPECS))
    parser.add_argument("--incorrect-only", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--blind-input", type=Path, default=Path("artifacts/blind/cota_adjudicated.jsonl"))
    args = parser.parse_args()

    joined_name = "joined_evaluation.reaggregated.jsonl"
    rows = [normalize_joined_record(row) for row in read_jsonl(args.run_dir / "restricted" / joined_name)]
    acceptance, prediction = POLICY_SPECS[args.policy]
    selected = [row for row in rows if row[acceptance]]
    if args.incorrect_only:
        selected = [row for row in selected if row[prediction] != row["reviewer_lot"]]
    primary = _details(args.run_dir)
    compare_rows = {}
    compare_details = {}
    if args.compare_run_dir:
        compare_rows = {row["case_key"]: normalize_joined_record(row) for row in read_jsonl(
            args.compare_run_dir / "restricted" / joined_name
        )}
        compare_details = _details(args.compare_run_dir)
    trajectories = _trajectories(args.blind_input)
    output_rows = []
    for row in selected:
        key = row["case_key"]
        detail = primary.get(key, {})
        prediction_row = detail.get("prediction", {})
        response = prediction_row.get("response", {})
        retrieval = detail.get("retrieval", {})
        diagnostic = {
            "case_key": key, "fold": row.get("fold"), "reviewer_lot": row["reviewer_lot"],
            "algorithm_lot": row["algorithm_lot"], "cota_lot": row["cota_lot"],
            "ai_lot": row.get("ai_lot"), "abstained": row["abstained"],
            "usable_ai_vote": row["usable_ai_vote"], "routing_reason": row["routing_reason"],
            "policy": args.policy, "policy_accepted": row[acceptance],
            "policy_correct": row[prediction] == row["reviewer_lot"] if row[acceptance] else None,
            "blind_trajectory": trajectories.get(key),
            "model_transition_decisions": response.get("transitions", []),
            "warnings": response.get("warnings", []),
            "retrieved_example_ids": prediction_row.get("retrieved_case_ids", []),
            "retrieval_hits": [{
                "rank": hit.get("retrieval_rank"), "retrieved_case_key": hit.get("retrieved_case_key"),
                "overall_similarity_score": hit.get("overall_similarity_score"),
                "component_similarity_scores": hit.get("component_similarity_scores"),
            } for hit in retrieval.get("hits", [])],
        }
        if key in compare_rows:
            other = compare_rows[key]
            other_detail = compare_details.get(key, {}).get("prediction", {})
            other_response = other_detail.get("response", {})
            diagnostic["comparison"] = {
                "ai_lot": other.get("ai_lot"), "abstained": other["abstained"],
                "usable_ai_vote": other["usable_ai_vote"], "routing_reason": other["routing_reason"],
                "policy_accepted": other[acceptance],
                "policy_correct": other[prediction] == other["reviewer_lot"] if other[acceptance] else None,
                "transition_decisions": other_response.get("transitions", []),
                "warnings": other_response.get("warnings", []),
                "retrieved_example_ids": other_detail.get("retrieved_case_ids", []),
            }
        output_rows.append(diagnostic)
    path = args.output / "restricted" / "diagnostic_cases.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in sorted(output_rows, key=lambda item: item["case_key"]):
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    print(f"Restricted diagnostic cases: {len(output_rows)} ({path})")


if __name__ == "__main__":
    main()
