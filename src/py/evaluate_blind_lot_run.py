#!/usr/bin/env python3
"""Reaggregate a completed blind LOT run without inference or retrieval."""

from __future__ import annotations

import argparse
from pathlib import Path

from blind_lot.reaggregation import reevaluate_run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    result = reevaluate_run(parse_args().run_dir)
    policies = result["metrics"]["routing_policy_analysis"]
    print("Evaluation-only reaggregation complete; provider calls: 0; retrieval executions: 0.")
    for name in (
        "algorithm_cota_agreement", "algorithm_ai_agreement",
        "cota_ai_agreement", "three_way_agreement",
    ):
        policy = policies[name]
        accuracy = policy["reviewer_accuracy_among_accepted"]["rate"]
        rendered = "N/A" if accuracy is None else f"{accuracy:.1%}"
        print(f"{name}: accepted={policy['accepted_count']}, correct={policy['correct_accepted_count']}, accuracy={rendered}")
    print(f"Aggregate output: {result['aggregate_output']}")


if __name__ == "__main__":
    main()
