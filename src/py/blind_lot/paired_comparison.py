"""Paired patient-level comparison of two completed blind LOT retrieval runs."""

from __future__ import annotations

import csv
import json
import math
import random
from pathlib import Path
from typing import Any, Callable

from lot_data import assert_aggregate_only_report, read_jsonl

from .evaluation import BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED, normalize_joined_record
from .reaggregation import REAGGREGATED_JOINED_NAME


COMPATIBLE_FIELDS = (
    "provider", "model", "reasoning_effort", "temperature", "prompt_version",
    "knowledge_version", "input_artifact_sha256", "random_seed", "bootstrap_replicates",
    "evaluation_version", "folds", "retrieval_training_folds", "run_purpose",
    "repeat_index",
)
POLICIES = {
    "algorithm_ai": ("accepted_by_algorithm_ai", "algorithm_lot"),
    "cota_ai": ("accepted_by_cota_ai", "cota_lot"),
    "three_way": ("accepted_by_three_way", "algorithm_lot"),
}


def _load_run(run_dir: Path, expected_k: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    metadata_path = run_dir / "restricted" / "experiment_metadata.json"
    joined_path = run_dir / "restricted" / REAGGREGATED_JOINED_NAME
    if not joined_path.exists():
        joined_path = run_dir / "restricted" / "joined_evaluation.jsonl"
    for path in (metadata_path, joined_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("retrieval_k") != expected_k:
        raise ValueError(f"expected retrieval_k={expected_k} in {run_dir}")
    rows = [normalize_joined_record(row) for row in read_jsonl(joined_path)]
    keys = [row.get("case_key") for row in rows]
    if None in keys or len(keys) != len(set(keys)):
        raise ValueError(f"missing or duplicate case key in {run_dir}")
    return metadata, rows


def load_aligned_runs(
    k3_run: Path, k5_run: Path,
) -> tuple[dict[str, Any], list[tuple[dict[str, Any], dict[str, Any]]]]:
    metadata3, rows3 = _load_run(k3_run, 3)
    metadata5, rows5 = _load_run(k5_run, 5)
    mismatches = [field for field in COMPATIBLE_FIELDS if metadata3.get(field) != metadata5.get(field)]
    if mismatches:
        raise ValueError(f"pooled runs have incompatible frozen metadata: {mismatches}")
    by3 = {row["case_key"]: row for row in rows3}
    by5 = {row["case_key"]: row for row in rows5}
    if set(by3) != set(by5):
        raise ValueError("k=3 and k=5 cohorts do not contain identical case keys")
    aligned = [(by3[key], by5[key]) for key in sorted(by3)]
    for left, right in aligned:
        if left.get("fold") != right.get("fold") or left["reviewer_lot"] != right["reviewer_lot"]:
            raise ValueError(f"paired truth/fold mismatch for {left['case_key']}")
    return metadata3, aligned


def _percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def _rate(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _condition_metric(rows: list[dict[str, Any]], name: str) -> float | None:
    if name == "generated_total_accuracy":
        cohort = [row for row in rows if row["has_generated_ai_total"]]
        return _rate(sum(row["ai_lot"] == row["reviewer_lot"] for row in cohort), len(cohort))
    if name == "usable_ai_coverage":
        return _rate(sum(row["usable_ai_vote"] for row in rows), len(rows))
    if name == "selective_ai_accuracy":
        cohort = [row for row in rows if row["usable_ai_vote"]]
        return _rate(sum(row["ai_lot"] == row["reviewer_lot"] for row in cohort), len(cohort))
    if name == "three_way_false_accept_rate":
        accepted = [row for row in rows if row["accepted_by_three_way"]]
        return _rate(
            sum(row["algorithm_lot"] != row["reviewer_lot"] for row in accepted),
            len(accepted),
        )
    policy, suffix = name.rsplit("_", 1)
    acceptance, prediction = POLICIES[policy]
    accepted = [row for row in rows if row[acceptance]]
    if suffix == "coverage":
        return _rate(len(accepted), len(rows))
    if suffix == "accuracy":
        return _rate(sum(row[prediction] == row["reviewer_lot"] for row in accepted), len(accepted))
    raise KeyError(name)


METRICS = (
    "generated_total_accuracy", "usable_ai_coverage", "selective_ai_accuracy",
    "three_way_coverage", "three_way_accuracy", "three_way_false_accept_rate",
    "algorithm_ai_coverage", "algorithm_ai_accuracy", "cota_ai_coverage", "cota_ai_accuracy",
)


def paired_bootstrap(
    aligned: list[tuple[dict[str, Any], dict[str, Any]]], *, seed: int = BOOTSTRAP_SEED,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, Any]:
    values = {name: [] for name in METRICS}
    if aligned:
        rng = random.Random(seed)
        for _ in range(replicates):
            sample = [aligned[rng.randrange(len(aligned))] for _ in aligned]
            left, right = [pair[0] for pair in sample], [pair[1] for pair in sample]
            for name in METRICS:
                k3, k5 = _condition_metric(left, name), _condition_metric(right, name)
                if k3 is not None and k5 is not None:
                    values[name].append(k5 - k3)
    result = {}
    left = [pair[0] for pair in aligned]
    right = [pair[1] for pair in aligned]
    for name in METRICS:
        point3, point5 = _condition_metric(left, name), _condition_metric(right, name)
        result[name] = {
            "k3": round(point3, 6) if point3 is not None else None,
            "k5": round(point5, 6) if point5 is not None else None,
            "k5_minus_k3": round(point5 - point3, 6) if point3 is not None and point5 is not None else None,
            "percentile_95_ci": {
                "lower": round(_percentile(values[name], 0.025), 6) if values[name] else None,
                "upper": round(_percentile(values[name], 0.975), 6) if values[name] else None,
                "replicates_used": len(values[name]),
            },
        }
    return result


def _both_categories(
    aligned: list[tuple[dict[str, Any], dict[str, Any]]],
    left: Callable[[dict[str, Any]], bool], right: Callable[[dict[str, Any]], bool],
) -> dict[str, int]:
    return {
        "both": sum(left(a) and right(b) for a, b in aligned),
        "k3_only": sum(left(a) and not right(b) for a, b in aligned),
        "k5_only": sum(not left(a) and right(b) for a, b in aligned),
        "neither": sum(not left(a) and not right(b) for a, b in aligned),
    }


def _restricted_rows(aligned: list[tuple[dict[str, Any], dict[str, Any]]]) -> list[dict[str, Any]]:
    output = []
    for left, right in aligned:
        reviewer = left["reviewer_lot"]
        output.append({
            "case_key": left["case_key"], "fold": left.get("fold"),
            "reviewer_lot": reviewer,
            "k3": {"ai_lot": left.get("ai_lot"), "abstained": left["abstained"],
                   "usable": left["usable_ai_vote"], "generated_total_correct": left["generated_total_correct"],
                   "usable_correct": left["non_abstained_prediction_correct"],
                   "three_way_accepted": left["three_way_agreement"],
                   "three_way_correct": left["ai_lot"] == reviewer if left["three_way_agreement"] else None},
            "k5": {"ai_lot": right.get("ai_lot"), "abstained": right["abstained"],
                   "usable": right["usable_ai_vote"], "generated_total_correct": right["generated_total_correct"],
                   "usable_correct": right["non_abstained_prediction_correct"],
                   "three_way_accepted": right["three_way_agreement"],
                   "three_way_correct": right["ai_lot"] == reviewer if right["three_way_agreement"] else None},
            "changes": {"newly_usable_under_k5": not left["usable_ai_vote"] and right["usable_ai_vote"],
                        "lost_usability_under_k5": left["usable_ai_vote"] and not right["usable_ai_vote"],
                        "gained_three_way_under_k5": not left["three_way_agreement"] and right["three_way_agreement"],
                        "lost_three_way_under_k5": left["three_way_agreement"] and not right["three_way_agreement"],
                        "k5_corrected_k3_generated_error": left["generated_total_correct"] is False and right["generated_total_correct"] is True,
                        "k5_corrected_k3_usable_error": left["non_abstained_prediction_correct"] is False and right["non_abstained_prediction_correct"] is True},
        })
    return output


def build_comparison(
    metadata: dict[str, Any], aligned: list[tuple[dict[str, Any], dict[str, Any]]],
    *, seed: int = BOOTSTRAP_SEED, replicates: int = BOOTSTRAP_REPLICATES,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    usable = _both_categories(aligned, lambda r: r["usable_ai_vote"], lambda r: r["usable_ai_vote"])
    generated_correct = _both_categories(
        aligned, lambda r: r["generated_total_correct"] is True,
        lambda r: r["generated_total_correct"] is True,
    )
    usable_correct = _both_categories(
        [(a, b) for a, b in aligned if a["usable_ai_vote"] and b["usable_ai_vote"]],
        lambda r: r["non_abstained_prediction_correct"] is True,
        lambda r: r["non_abstained_prediction_correct"] is True,
    )
    accepted = _both_categories(aligned, lambda r: r["three_way_agreement"], lambda r: r["three_way_agreement"])
    false_accepted = _both_categories(
        aligned,
        lambda r: r["three_way_agreement"] and r["ai_lot"] != r["reviewer_lot"],
        lambda r: r["three_way_agreement"] and r["ai_lot"] != r["reviewer_lot"],
    )
    false3 = {a["case_key"] for a, _ in aligned if a["three_way_agreement"] and a["ai_lot"] != a["reviewer_lot"]}
    false5 = {b["case_key"] for _, b in aligned if b["three_way_agreement"] and b["ai_lot"] != b["reviewer_lot"]}
    newly = [(a, b) for a, b in aligned if not a["usable_ai_vote"] and b["usable_ai_vote"]]
    lost = [(a, b) for a, b in aligned if a["usable_ai_vote"] and not b["usable_ai_vote"]]
    k3_usable_correct = sum(
        a["usable_ai_vote"] and a["ai_lot"] == a["reviewer_lot"] for a, _ in aligned
    )
    k5_usable_correct = sum(
        b["usable_ai_vote"] and b["ai_lot"] == b["reviewer_lot"] for _, b in aligned
    )
    k3_usable_incorrect = sum(
        a["usable_ai_vote"] and a["ai_lot"] != a["reviewer_lot"] for a, _ in aligned
    )
    k5_usable_incorrect = sum(
        b["usable_ai_vote"] and b["ai_lot"] != b["reviewer_lot"] for _, b in aligned
    )
    analyses = [
        ("all_folds", lambda fold: True, False),
        ("folds_1_through_4", lambda fold: fold != 0, False),
        ("leave_fold_0_out", lambda fold: fold != 0, True),
        *[(f"leave_fold_{fold}_out", lambda value, excluded=fold: value != excluded, False) for fold in range(1, 5)],
    ]
    bootstrap = {}
    for name, predicate, development_label in analyses:
        cohort = [(a, b) for a, b in aligned if predicate(a.get("fold"))]
        bootstrap[name] = {
            "patient_count": len(cohort),
            "fold_0_label": "development-influenced" if name == "all_folds" else (
                "development-influenced fold excluded" if development_label else None
            ),
            "metrics": paired_bootstrap(cohort, seed=seed, replicates=replicates),
        }
    report = {
        "schema_version": "blind-lot-paired-comparison-1.0.0",
        "report_scope": "aggregate_only",
        "comparison": "paired k=5 minus k=3 retrieval analysis",
        "patient_count": len(aligned),
        "frozen_configuration": {field: metadata.get(field) for field in COMPATIBLE_FIELDS},
        "bootstrap": {"method": "paired patient-level nonparametric percentile bootstrap",
                      "seed": seed, "requested_replicates": replicates, "analyses": bootstrap},
        "paired_categories": {
            "abstention_and_usability": {
                "both_abstain": sum(a["abstained"] and b["abstained"] for a, b in aligned),
                "k3_abstains_k5_usable": sum(a["abstained"] and b["usable_ai_vote"] for a, b in aligned),
                "k3_usable_k5_abstains": sum(a["usable_ai_vote"] and b["abstained"] for a, b in aligned),
                "both_usable_same_total": sum(a["usable_ai_vote"] and b["usable_ai_vote"] and a["ai_lot"] == b["ai_lot"] for a, b in aligned),
                "both_usable_different_totals": sum(a["usable_ai_vote"] and b["usable_ai_vote"] and a["ai_lot"] != b["ai_lot"] for a, b in aligned),
            },
            "generated_total_correctness": {"correct_under_both": generated_correct["both"],
                "correct_only_k3": generated_correct["k3_only"], "correct_only_k5": generated_correct["k5_only"],
                "wrong_under_both": generated_correct["neither"]},
            "usable_prediction_correctness_when_comparable": {"correct_under_both": usable_correct["both"],
                "correct_only_k3": usable_correct["k3_only"], "correct_only_k5": usable_correct["k5_only"],
                "wrong_under_both": usable_correct["neither"]},
            "three_way_acceptance": {"accepted_by_both": accepted["both"], "accepted_only_k3": accepted["k3_only"],
                "accepted_only_k5": accepted["k5_only"], "rejected_by_both": accepted["neither"]},
            "three_way_false_acceptance": {"false_accepted_by_both": false_accepted["both"],
                "false_accepted_only_k3": false_accepted["k3_only"], "false_accepted_only_k5": false_accepted["k5_only"]},
        },
        "direct_answers": {
            "same_false_three_way_patient": len(false3) == len(false5) == 1 and false3 == false5,
            "newly_usable_under_k5_count": len(newly),
            "newly_usable_under_k5_correct_count": sum(b["ai_lot"] == b["reviewer_lot"] for _, b in newly),
            "newly_usable_under_k5_incorrect_count": sum(b["ai_lot"] != b["reviewer_lot"] for _, b in newly),
            "lost_usable_under_k5_count": len(lost),
            "net_additional_k5_usable_votes": usable["k5_only"] - usable["k3_only"],
            "net_additional_k5_usable_correct_votes": k5_usable_correct - k3_usable_correct,
            "net_additional_k5_usable_incorrect_votes": k5_usable_incorrect - k3_usable_incorrect,
            "three_way_gained_under_k5_count": accepted["k5_only"],
            "three_way_lost_under_k5_count": accepted["k3_only"],
            "k5_corrected_any_k3_generated_total_error": generated_correct["k5_only"] > 0,
            "k5_corrected_k3_generated_total_error_count": generated_correct["k5_only"],
            "k5_corrected_any_comparable_k3_usable_vote_error": usable_correct["k5_only"] > 0,
            "k5_corrected_comparable_k3_usable_vote_error_count": usable_correct["k5_only"],
            "abstention_stability": {"both_abstain_count": sum(a["abstained"] and b["abstained"] for a, b in aligned),
                "same_abstention_state_count": sum(a["abstained"] == b["abstained"] for a, b in aligned),
                "same_abstention_state_rate": round(sum(a["abstained"] == b["abstained"] for a, b in aligned) / len(aligned), 6) if aligned else None},
            "case_lists": "Available only in restricted paired_case_comparison.jsonl via change flags.",
        },
        "fold_0_note": "Fold 0 was development-influenced.",
        "small_subset_warning": "Small accepted subsets must not be overinterpreted.",
        "pairwise_policy_note": "Pairwise policies are exploratory.",
    }
    assert_aggregate_only_report(report)
    return report, _restricted_rows(aligned)


def write_comparison(
    k3_run: Path, k5_run: Path, output: Path, *, seed: int = BOOTSTRAP_SEED,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, Path]:
    metadata, aligned = load_aligned_runs(k3_run, k5_run)
    report, restricted_rows = build_comparison(metadata, aligned, seed=seed, replicates=replicates)
    restricted_path = output / "restricted" / "paired_case_comparison.jsonl"
    public_json = output / "public" / "paired_retrieval_comparison.json"
    public_csv = output / "public" / "paired_retrieval_comparison.csv"
    restricted_path.parent.mkdir(parents=True, exist_ok=True)
    with restricted_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in restricted_rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    public_json.parent.mkdir(parents=True, exist_ok=True)
    public_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    csv_rows = []
    for group, values in report["paired_categories"].items():
        csv_rows.extend({"section": group, "metric": name, "value": value} for name, value in values.items())
    csv_rows.extend({"section": "direct_answers", "metric": name, "value": value}
                    for name, value in report["direct_answers"].items() if not isinstance(value, dict))
    with public_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("section", "metric", "value"))
        writer.writeheader()
        writer.writerows(csv_rows)
    return {"restricted": restricted_path, "json": public_json, "csv": public_csv}
