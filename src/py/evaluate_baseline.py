"""Evaluate the current COTA LOT algorithm against reviewer consensus."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from lot_data import read_public_jsonl, write_json
from textbook_algo_cota import lot_algorithm_cota


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts" / "normalized" / "cota_patients.jsonl"
DEFAULT_OUTPUT = ROOT / "artifacts" / "reports" / "baseline_evaluation.json"


def ratio(numerator: int, denominator: int) -> float | None:
    return round(numerator / denominator, 6) if denominator else None


def evaluate(records: list[dict[str, Any]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    excluded = 0
    for record in records:
        truth = record["reviewer_lot"]["consensus"]
        cota_lot = record["cota_lot"]
        if truth is None or cota_lot is None:
            excluded += 1
            continue
        sequence = [frozenset(event["drugs"]) for event in record["trajectory"]]
        prediction, flags = lot_algorithm_cota(sequence)
        error = prediction - truth
        agrees = prediction == cota_lot
        rows.append({
            "patient_key": record["patient_key"],
            "reviewer_consensus_lot": truth,
            "algorithm_lot": prediction,
            "cota_lot": cota_lot,
            "error": error,
            "algorithm_cota_agree": agrees,
            "algorithm_flags": sorted(flags),
        })

    n = len(rows)
    exact = sum(row["error"] == 0 for row in rows)
    within_one = sum(abs(row["error"]) <= 1 for row in rows)
    over = sum(row["error"] > 0 for row in rows)
    under = sum(row["error"] < 0 for row in rows)
    agreed = [row for row in rows if row["algorithm_cota_agree"]]
    agreed_correct = sum(row["error"] == 0 for row in agreed)
    both_agree_wrong = [row for row in agreed if row["error"] != 0]

    return {
        "schema_version": "1.0.0",
        "ground_truth": "reviewer consensus",
        "eligible_patients": n,
        "excluded_patients": excluded,
        "metrics": {
            "exact_accuracy": {"count": exact, "rate": ratio(exact, n)},
            "within_one_accuracy": {"count": within_one, "rate": ratio(within_one, n)},
            "over_count": {"count": over, "rate": ratio(over, n)},
            "under_count": {"count": under, "rate": ratio(under, n)},
            "algorithm_cota_agreement_coverage": {
                "count": len(agreed), "rate": ratio(len(agreed), n),
            },
            "accuracy_when_algorithm_cota_agree": {
                "count": agreed_correct, "denominator": len(agreed),
                "rate": ratio(agreed_correct, len(agreed)),
            },
            "both_agree_but_wrong": {
                "count": len(both_agree_wrong),
                "rate_all_patients": ratio(len(both_agree_wrong), n),
                "rate_agreed_patients": ratio(len(both_agree_wrong), len(agreed)),
            },
        },
        "both_agree_but_wrong_cases": both_agree_wrong,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patients", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = evaluate(read_public_jsonl(args.patients))
    write_json(args.output, results)
    metrics = results["metrics"]
    print(f"Evaluated {results['eligible_patients']} patients against reviewer consensus.")
    print(f"Exact accuracy: {metrics['exact_accuracy']['rate']:.1%}")
    print(f"Within one: {metrics['within_one_accuracy']['rate']:.1%}")
    print(f"Results written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
