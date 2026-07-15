"""Aggregate reviewer-consensus evaluation and patient-level bootstrap intervals."""

from __future__ import annotations

import math
import random
from dataclasses import asdict, dataclass
from typing import Any, Callable


BOOTSTRAP_SEED = 4992026
BOOTSTRAP_REPLICATES = 2000
EVALUATION_VERSION = "blind-lot-evaluation-1.2.0"


def is_valid_lot(value: Any) -> bool:
    """Return whether ``value`` is a normalized, non-negative LOT total."""
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


@dataclass(frozen=True)
class PredictionStatus:
    """Canonical derived status for one AI prediction and its three-vote routing."""

    has_generated_ai_total: bool
    usable_ai_vote: bool
    generated_total_correct: bool | None
    non_abstained_prediction_correct: bool | None
    abstained: bool
    invalid_ai_output: bool
    ai_algorithm_agreement: bool
    ai_cota_agreement: bool
    algorithm_cota_agreement: bool
    three_way_agreement: bool
    routing_decision: str
    routing_reason: str


def derive_prediction_status(
    *,
    ai_lot_count: Any,
    abstained: Any,
    reviewer_lot: int,
    algorithm_lot_count: int,
    cota_lot_count: int,
) -> PredictionStatus:
    """Derive validity, correctness, agreement, and routing exactly once.

    Zero is valid because the model response contract uses it for an empty trajectory.
    Python booleans are explicitly not LOT totals despite being ``int`` subclasses.
    """
    has_generated = is_valid_lot(ai_lot_count)
    invalid = not has_generated or not isinstance(abstained, bool)
    is_abstained = abstained is True
    usable = has_generated and abstained is False
    generated_correct = ai_lot_count == reviewer_lot if has_generated else None
    usable_correct = ai_lot_count == reviewer_lot if usable else None
    ai_algorithm = (
        usable and is_valid_lot(algorithm_lot_count)
        and ai_lot_count == algorithm_lot_count
    )
    ai_cota = (
        usable and is_valid_lot(cota_lot_count)
        and ai_lot_count == cota_lot_count
    )
    algorithm_cota = (
        is_valid_lot(algorithm_lot_count)
        and is_valid_lot(cota_lot_count)
        and algorithm_lot_count == cota_lot_count
    )
    three_way = ai_algorithm and ai_cota

    if invalid:
        decision, reason = "human_review", "invalid_ai_output"
    elif is_abstained:
        decision, reason = "human_review", "ai_abstained"
    elif not three_way:
        decision, reason = "human_review", "three_way_disagreement"
    else:
        decision, reason = "consensus_candidate", "three_way_consensus"

    return PredictionStatus(
        has_generated_ai_total=has_generated,
        usable_ai_vote=usable,
        generated_total_correct=generated_correct,
        non_abstained_prediction_correct=usable_correct,
        abstained=is_abstained,
        invalid_ai_output=invalid,
        ai_algorithm_agreement=ai_algorithm,
        ai_cota_agreement=ai_cota,
        algorithm_cota_agreement=algorithm_cota,
        three_way_agreement=three_way,
        routing_decision=decision,
        routing_reason=reason,
    )


def normalize_joined_record(row: dict[str, Any]) -> dict[str, Any]:
    """Return a joined record with canonical prediction-status fields attached."""
    status = derive_prediction_status(
        ai_lot_count=row["ai_lot"],
        abstained=row["ai_abstained"],
        reviewer_lot=row["reviewer_lot"],
        algorithm_lot_count=row["algorithm_lot"],
        cota_lot_count=row["cota_lot"],
    )
    normalized = {**row, **asdict(status)}
    accepted = {
        "algorithm_only": is_valid_lot(row.get("algorithm_lot")),
        "cota_only": is_valid_lot(row.get("cota_lot")),
        "usable_ai_only": status.usable_ai_vote,
        "algorithm_cota": status.algorithm_cota_agreement,
        "algorithm_ai": status.ai_algorithm_agreement,
        "cota_ai": status.ai_cota_agreement,
        "three_way": status.three_way_agreement,
    }
    prediction_fields = {
        "algorithm_only": "algorithm_lot",
        "cota_only": "cota_lot",
        "usable_ai_only": "ai_lot",
        "algorithm_cota": "algorithm_lot",
        "algorithm_ai": "algorithm_lot",
        "cota_ai": "cota_lot",
        "three_way": "algorithm_lot",
    }
    for name, does_accept in accepted.items():
        normalized[f"accepted_by_{name}"] = does_accept
        correctness_name = (
            f"{name}_correct" if name.endswith("_only")
            else f"{name}_policy_correct"
        )
        normalized[correctness_name] = (
            normalized[prediction_fields[name]] == normalized["reviewer_lot"]
            if does_accept else None
        )
    normalized["vote_pattern"] = classify_vote_pattern(normalized)
    return normalized


def classify_vote_pattern(row: dict[str, Any]) -> str:
    """Classify one patient into an exhaustive, mutually exclusive vote stratum."""
    algorithm = row.get("algorithm_lot")
    cota = row.get("cota_lot")
    ai = row.get("ai_lot")
    algorithm_valid = is_valid_lot(algorithm)
    cota_valid = is_valid_lot(cota)

    if row.get("invalid_ai_output"):
        return "invalid_ai_output"
    if row.get("abstained"):
        if algorithm_valid and cota_valid:
            if algorithm == cota:
                return "algorithm_cota_agree_ai_abstains"
            return "algorithm_cota_differ_ai_abstains"
        return "ai_abstains_missing_non_ai_vote"
    if row.get("usable_ai_vote") and algorithm_valid and cota_valid:
        if algorithm == ai == cota:
            return "all_three_agree"
        if algorithm == ai:
            return "algorithm_ai_agree_cota_differs"
        if cota == ai:
            return "cota_ai_agree_algorithm_differs"
        if algorithm == cota:
            return "algorithm_cota_agree_ai_differs"
        return "all_three_differ"
    return "usable_ai_missing_non_ai_vote"


def _rate(numerator: int, denominator: int) -> float | None:
    return round(numerator / denominator, 6) if denominator else None


def _metric(count: int, denominator: int) -> dict[str, Any]:
    return {"count": count, "denominator": denominator, "rate": _rate(count, denominator)}


def _percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def bootstrap_ci(
    rows: list[dict[str, Any]],
    statistic: Callable[[list[dict[str, Any]]], float | None],
    *,
    seed: int = BOOTSTRAP_SEED,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, Any]:
    """Resample eligible patients, then apply subset filtering in ``statistic``."""
    if not rows:
        return {"lower": None, "upper": None, "replicates_used": 0}
    rng = random.Random(seed)
    values: list[float] = []
    for _ in range(replicates):
        sample = [rows[rng.randrange(len(rows))] for _ in rows]
        value = statistic(sample)
        if value is not None:
            values.append(value)
    if not values:
        return {"lower": None, "upper": None, "replicates_used": 0}
    return {
        "lower": round(_percentile(values, 0.025), 6),
        "upper": round(_percentile(values, 0.975), 6),
        "replicates_used": len(values),
    }


def _accuracy(rows: list[dict[str, Any]], prediction: str) -> float | None:
    return sum(row[prediction] == row["reviewer_lot"] for row in rows) / len(rows) if rows else None


def _coverage(rows: list[dict[str, Any]], predicate: Callable[[dict[str, Any]], bool]) -> float | None:
    return sum(predicate(row) for row in rows) / len(rows) if rows else None


POLICY_SPECS: dict[str, tuple[str, str]] = {
    "algorithm_only": ("accepted_by_algorithm_only", "algorithm_lot"),
    "cota_only": ("accepted_by_cota_only", "cota_lot"),
    "usable_ai_only": ("accepted_by_usable_ai_only", "ai_lot"),
    "algorithm_cota_agreement": ("accepted_by_algorithm_cota", "algorithm_lot"),
    "algorithm_ai_agreement": ("accepted_by_algorithm_ai", "algorithm_lot"),
    "cota_ai_agreement": ("accepted_by_cota_ai", "cota_lot"),
    "three_way_agreement": ("accepted_by_three_way", "algorithm_lot"),
}

POLICY_DEFINITIONS = {
    "algorithm_only": "valid algorithm LOT",
    "cota_only": "valid COTA LOT",
    "usable_ai_only": "valid non-abstained AI LOT",
    "algorithm_cota_agreement": "valid algorithm LOT = valid COTA LOT; AI ignored",
    "algorithm_ai_agreement": "usable AI and valid algorithm LOT = AI LOT; COTA ignored",
    "cota_ai_agreement": "usable AI and valid COTA LOT = AI LOT; algorithm ignored",
    "three_way_agreement": "usable AI and algorithm LOT = AI LOT = COTA LOT",
}

VOTE_PATTERN_NAMES = (
    "all_three_agree",
    "algorithm_ai_agree_cota_differs",
    "cota_ai_agree_algorithm_differs",
    "algorithm_cota_agree_ai_differs",
    "algorithm_cota_agree_ai_abstains",
    "all_three_differ",
    "algorithm_cota_differ_ai_abstains",
    "invalid_ai_output",
    "ai_abstains_missing_non_ai_vote",
    "usable_ai_missing_non_ai_vote",
)


def _policy_metric(
    rows: list[dict[str, Any]],
    *,
    acceptance_field: str,
    prediction_field: str,
    bootstrap_seed: int,
    bootstrap_replicates: int,
) -> dict[str, Any]:
    accepted = [row for row in rows if row[acceptance_field]]
    correct = [row for row in accepted if row[prediction_field] == row["reviewer_lot"]]
    incorrect_count = len(accepted) - len(correct)
    result = {
        "eligible_patient_count": len(rows),
        "accepted_count": len(accepted),
        "coverage": _metric(len(accepted), len(rows)),
        "correct_accepted_count": len(correct),
        "incorrect_accepted_count": incorrect_count,
        "reviewer_accuracy_among_accepted": _metric(len(correct), len(accepted)),
        "false_accept_rate_among_accepted": _metric(incorrect_count, len(accepted)),
        "routed_to_review_count": len(rows) - len(accepted),
        "review_proportion": _metric(len(rows) - len(accepted), len(rows)),
    }
    result["coverage"]["bootstrap_95_ci"] = bootstrap_ci(
        rows,
        lambda sample: _coverage(sample, lambda row: row[acceptance_field]),
        seed=bootstrap_seed,
        replicates=bootstrap_replicates,
    )
    result["reviewer_accuracy_among_accepted"]["bootstrap_95_ci"] = bootstrap_ci(
        rows,
        lambda sample: _accuracy(
            [row for row in sample if row[acceptance_field]], prediction_field
        ),
        seed=bootstrap_seed + 1,
        replicates=bootstrap_replicates,
    )
    result["false_accept_rate_among_accepted"]["bootstrap_95_ci"] = bootstrap_ci(
        rows,
        lambda sample: (
            1.0 - value
            if (value := _accuracy(
                [row for row in sample if row[acceptance_field]], prediction_field
            )) is not None else None
        ),
        seed=bootstrap_seed + 2,
        replicates=bootstrap_replicates,
    )
    return result


def _policy_delta(
    rows: list[dict[str, Any]],
    policy_name: str,
    baseline_name: str,
) -> dict[str, Any]:
    acceptance_field, prediction_field = POLICY_SPECS[policy_name]
    baseline_field, baseline_prediction = POLICY_SPECS[baseline_name]
    accepted = [row for row in rows if row[acceptance_field]]
    baseline = [row for row in rows if row[baseline_field]]
    accepted_by_both = [row for row in rows if row[acceptance_field] and row[baseline_field]]
    policy_only = [row for row in rows if row[acceptance_field] and not row[baseline_field]]
    baseline_only = [row for row in rows if row[baseline_field] and not row[acceptance_field]]
    policy_accuracy = _accuracy(accepted, prediction_field)
    baseline_accuracy = _accuracy(baseline, baseline_prediction)
    return {
        "accepted_by_both_count": len(accepted_by_both),
        "accepted_only_by_policy_count": len(policy_only),
        "accepted_only_by_baseline_count": len(baseline_only),
        "net_accepted_count_difference": len(accepted) - len(baseline),
        "correct_by_both_count": sum(
            row[prediction_field] == row["reviewer_lot"]
            and row[baseline_prediction] == row["reviewer_lot"]
            for row in accepted_by_both
        ),
        "incorrect_by_both_count": sum(
            row[prediction_field] != row["reviewer_lot"]
            and row[baseline_prediction] != row["reviewer_lot"]
            for row in accepted_by_both
        ),
        "correct_only_by_policy_count": sum(
            row[prediction_field] == row["reviewer_lot"] for row in policy_only
        ),
        "incorrect_only_by_policy_count": sum(
            row[prediction_field] != row["reviewer_lot"] for row in policy_only
        ),
        "correct_only_by_baseline_count": sum(
            row[baseline_prediction] == row["reviewer_lot"] for row in baseline_only
        ),
        "incorrect_only_by_baseline_count": sum(
            row[baseline_prediction] != row["reviewer_lot"] for row in baseline_only
        ),
        "additional_accepted_patients": len(policy_only),
        "additional_accepted_patients_deprecated": True,
        "additional_accepted_patients_deprecated_alias_for": "accepted_only_by_policy_count",
        "additional_correct_accepted_patients": sum(
            row[prediction_field] == row["reviewer_lot"] for row in policy_only
        ),
        "additional_incorrect_accepted_patients": sum(
            row[prediction_field] != row["reviewer_lot"] for row in policy_only
        ),
        "coverage_difference": round(
            (len(accepted) - len(baseline)) / len(rows), 6
        ) if rows else None,
        "accuracy_difference": (
            round(policy_accuracy - baseline_accuracy, 6)
            if policy_accuracy is not None and baseline_accuracy is not None else None
        ),
    }


def _accuracy_metric_for_source(
    rows: list[dict[str, Any]], source: str, usable_field: str | None = None
) -> dict[str, Any]:
    eligible = [row for row in rows if usable_field is None or row[usable_field]]
    eligible = [row for row in eligible if is_valid_lot(row.get(source))]
    correct = sum(row[source] == row["reviewer_lot"] for row in eligible)
    return _metric(correct, len(eligible))


def _vote_pattern_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(rows)
    patterns: dict[str, Any] = {}
    for name in VOTE_PATTERN_NAMES:
        cohort = [row for row in rows if row["vote_pattern"] == name]
        source_unique = {"algorithm": 0, "cota": 0, "ai": 0}
        none = 0
        multiple = 0
        majority_correct = dissent_correct = 0
        majority_sources: tuple[str, str] | None = None
        dissent_source: str | None = None
        if name == "algorithm_ai_agree_cota_differs":
            majority_sources, dissent_source = ("algorithm_lot", "ai_lot"), "cota_lot"
        elif name == "cota_ai_agree_algorithm_differs":
            majority_sources, dissent_source = ("cota_lot", "ai_lot"), "algorithm_lot"
        elif name == "algorithm_cota_agree_ai_differs":
            majority_sources, dissent_source = ("algorithm_lot", "cota_lot"), "ai_lot"

        match_counts = {"algorithm": 0, "cota": 0, "ai": 0}
        for row in cohort:
            matches = []
            for label, field in (
                ("algorithm", "algorithm_lot"), ("cota", "cota_lot"), ("ai", "ai_lot")
            ):
                available = is_valid_lot(row.get(field)) and (
                    label != "ai" or row["usable_ai_vote"]
                )
                if available and row[field] == row["reviewer_lot"]:
                    matches.append(label)
                    match_counts[label] += 1
            if not matches:
                none += 1
            elif len(matches) > 1:
                multiple += 1
            else:
                source_unique[matches[0]] += 1
            if majority_sources is not None:
                if row[majority_sources[0]] == row["reviewer_lot"]:
                    majority_correct += 1
                elif row[dissent_source] == row["reviewer_lot"]:
                    dissent_correct += 1

        entry = {
            "patient_count": len(cohort),
            "proportion_of_eligible_patients": _metric(len(cohort), n),
            "reviewer_accuracy": {
                "algorithm": _accuracy_metric_for_source(cohort, "algorithm_lot"),
                "cota": _accuracy_metric_for_source(cohort, "cota_lot"),
                "ai_when_usable": _accuracy_metric_for_source(cohort, "ai_lot", "usable_ai_vote"),
            },
            "reviewer_match_counts": match_counts,
            "none_matched_reviewer_count": none,
            "more_than_one_matched_reviewer_count": multiple,
            "winner_attribution": {
                "algorithm_uniquely_correct": source_unique["algorithm"],
                "cota_uniquely_correct": source_unique["cota"],
                "ai_uniquely_correct": source_unique["ai"],
                "none_correct": none,
            },
        }
        if majority_sources is not None:
            entry["majority_vs_dissent"] = {
                "majority_correct": majority_correct,
                "dissenting_vote_correct": dissent_correct,
                "none_correct": none,
            }
        patterns[name] = entry
    return {
        "patterns": patterns,
        "patient_count_sum": sum(item["patient_count"] for item in patterns.values()),
        "eligible_patient_count": n,
        "mutually_exclusive_and_exhaustive": sum(
            item["patient_count"] for item in patterns.values()
        ) == n,
    }


def _conditional_subgroup(
    cohort: list[dict[str, Any]], total_eligible: int, prediction_field: str
) -> dict[str, Any]:
    correct = sum(row[prediction_field] == row["reviewer_lot"] for row in cohort)
    return {
        "count": len(cohort),
        "coverage": _metric(len(cohort), total_eligible),
        "reviewer_accuracy": _metric(correct, len(cohort)),
        "correct_count": correct,
        "incorrect_count": len(cohort) - correct,
    }


def _conditional_agreement_analysis(rows: list[dict[str, Any]]) -> dict[str, Any]:
    algorithm_ai = [row for row in rows if row["accepted_by_algorithm_ai"]]
    algorithm_cota = [row for row in rows if row["accepted_by_algorithm_cota"]]
    cota_ai = [row for row in rows if row["accepted_by_cota_ai"]]
    n = len(rows)
    return {
        "algorithm_ai_agreement_partitioned_by_cota": {
            "base_count": len(algorithm_ai),
            "cota_agrees": _conditional_subgroup(
                [row for row in algorithm_ai if row["cota_lot"] == row["algorithm_lot"]],
                n, "algorithm_lot",
            ),
            "cota_disagrees": _conditional_subgroup(
                [row for row in algorithm_ai if row["cota_lot"] != row["algorithm_lot"]],
                n, "algorithm_lot",
            ),
        },
        "algorithm_cota_agreement_partitioned_by_ai": {
            "base_count": len(algorithm_cota),
            "ai_agrees": _conditional_subgroup(
                [row for row in algorithm_cota if row["usable_ai_vote"] and row["ai_lot"] == row["algorithm_lot"]],
                n, "algorithm_lot",
            ),
            "ai_disagrees": _conditional_subgroup(
                [row for row in algorithm_cota if row["usable_ai_vote"] and row["ai_lot"] != row["algorithm_lot"]],
                n, "algorithm_lot",
            ),
            "ai_abstains": _conditional_subgroup(
                [row for row in algorithm_cota if row["abstained"] and not row["invalid_ai_output"]],
                n, "algorithm_lot",
            ),
            "ai_invalid": _conditional_subgroup(
                [row for row in algorithm_cota if row["invalid_ai_output"]],
                n, "algorithm_lot",
            ),
        },
        "cota_ai_agreement_partitioned_by_algorithm": {
            "base_count": len(cota_ai),
            "algorithm_agrees": _conditional_subgroup(
                [row for row in cota_ai if row["algorithm_lot"] == row["cota_lot"]],
                n, "cota_lot",
            ),
            "algorithm_disagrees": _conditional_subgroup(
                [row for row in cota_ai if row["algorithm_lot"] != row["cota_lot"]],
                n, "cota_lot",
            ),
        },
    }


def routing_policy_analysis(
    rows: list[dict[str, Any]], *, bootstrap_seed: int, bootstrap_replicates: int
) -> dict[str, Any]:
    policies: dict[str, Any] = {}
    for offset, (name, (acceptance, prediction)) in enumerate(POLICY_SPECS.items()):
        policies[name] = _policy_metric(
            rows,
            acceptance_field=acceptance,
            prediction_field=prediction,
            bootstrap_seed=bootstrap_seed + 100 + offset * 3,
            bootstrap_replicates=bootstrap_replicates,
        )
        policies[name]["definition"] = POLICY_DEFINITIONS[name]
    for name in POLICY_SPECS:
        if name == "three_way_agreement":
            continue
        comparisons = {
            "versus_three_way_agreement": _policy_delta(
                rows, name, "three_way_agreement"
            )
        }
        if name != "algorithm_cota_agreement":
            comparisons["versus_algorithm_cota_agreement"] = _policy_delta(
                rows, name, "algorithm_cota_agreement"
            )
        policies[name]["comparisons"] = comparisons
    return {
        "primary_policy": "three_way_agreement",
        "primary_policy_definition": "usable AI and algorithm LOT = AI LOT = COTA LOT",
        "secondary_exploratory_policies": [
            name for name in POLICY_SPECS if name != "three_way_agreement"
        ],
        "exploratory_warning": (
            "Pairwise and single-source analyses are exploratory after fold-0 inspection; "
            "they do not redefine the primary research endpoint."
        ),
        **policies,
    }


def _prediction_metrics(
    cohort: list[dict[str, Any]],
    *,
    coverage_denominator: int,
) -> dict[str, Any]:
    count = len(cohort)
    exact = sum(row["ai_lot"] == row["reviewer_lot"] for row in cohort)
    within = sum(abs(row["ai_lot"] - row["reviewer_lot"]) <= 1 for row in cohort)
    over = sum(row["ai_lot"] > row["reviewer_lot"] for row in cohort)
    under = sum(row["ai_lot"] < row["reviewer_lot"] for row in cohort)
    mae = sum(abs(row["ai_lot"] - row["reviewer_lot"]) for row in cohort) / count if count else None
    return {
        "coverage": _metric(count, coverage_denominator),
        "exact_accuracy": _metric(exact, count),
        "within_one_accuracy": _metric(within, count),
        "mean_absolute_error": round(mae, 6) if mae is not None else None,
        "over_count": _metric(over, count),
        "under_count": _metric(under, count),
    }


def evaluate_joined(
    rows: list[dict[str, Any]],
    *,
    bootstrap_seed: int = BOOTSTRAP_SEED,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
) -> dict[str, Any]:
    """Evaluate eligible adjudicated cases using canonical AI vote semantics."""
    rows = [normalize_joined_record(row) for row in rows]
    n = len(rows)
    generated = [row for row in rows if row["has_generated_ai_total"]]
    usable = [row for row in rows if row["usable_ai_vote"]]
    abstained_generated = [row for row in rows if row["abstained"] and row["has_generated_ai_total"]]

    generated_metrics = _prediction_metrics(generated, coverage_denominator=n)
    selective_metrics = _prediction_metrics(usable, coverage_denominator=n)
    generated_metrics["coverage"]["bootstrap_95_ci"] = bootstrap_ci(
        rows,
        lambda sample: _coverage(sample, lambda row: row["has_generated_ai_total"]),
        seed=bootstrap_seed + 39,
        replicates=bootstrap_replicates,
    )
    generated_metrics["exact_accuracy"]["bootstrap_95_ci"] = bootstrap_ci(
        rows,
        lambda sample: _accuracy([row for row in sample if row["has_generated_ai_total"]], "ai_lot"),
        seed=bootstrap_seed + 40,
        replicates=bootstrap_replicates,
    )
    selective_metrics["coverage"]["bootstrap_95_ci"] = bootstrap_ci(
        rows,
        lambda sample: _coverage(sample, lambda row: row["usable_ai_vote"]),
        seed=bootstrap_seed + 41,
        replicates=bootstrap_replicates,
    )
    selective_metrics["exact_accuracy"]["bootstrap_95_ci"] = bootstrap_ci(
        rows,
        lambda sample: _accuracy([row for row in sample if row["usable_ai_vote"]], "ai_lot"),
        seed=bootstrap_seed + 42,
        replicates=bootstrap_replicates,
    )

    pair_specs = {
        "ai_cota": ("ai_cota_agreement", "ai_lot"),
        "ai_algorithm": ("ai_algorithm_agreement", "ai_lot"),
        "algorithm_cota": ("algorithm_cota_agreement", "algorithm_lot"),
    }
    pairwise: dict[str, Any] = {}
    for offset, (name, (agreement_field, prediction_field)) in enumerate(pair_specs.items()):
        agreed = [row for row in rows if row[agreement_field]]
        correct = sum(row[prediction_field] == row["reviewer_lot"] for row in agreed)
        pairwise[name] = {
            "agreement_coverage": _metric(len(agreed), n),
            "agreement_coverage_denominator_description": "all eligible patients",
            "reviewer_accuracy_when_agree": _metric(correct, len(agreed)),
        }
        pairwise[name]["agreement_coverage"]["bootstrap_95_ci"] = bootstrap_ci(
            rows,
            lambda sample, field=agreement_field: _coverage(sample, lambda row: row[field]),
            seed=bootstrap_seed + offset,
            replicates=bootstrap_replicates,
        )
        pairwise[name]["reviewer_accuracy_when_agree"]["bootstrap_95_ci"] = bootstrap_ci(
            rows,
            lambda sample, field=agreement_field, prediction=prediction_field: _accuracy(
                [row for row in sample if row[field]], prediction
            ),
            seed=bootstrap_seed + 10 + offset,
            replicates=bootstrap_replicates,
        )

    three = [row for row in rows if row["three_way_agreement"]]
    three_correct = sum(row["ai_lot"] == row["reviewer_lot"] for row in three)
    three_wrong = len(three) - three_correct
    three_way = {
        "agreement": _metric(len(three), n),
        "reviewer_accuracy": _metric(three_correct, len(three)),
        "all_three_agree_but_wrong_count": three_wrong,
        "human_review": _metric(n - len(three), n),
    }
    three_way["agreement"]["bootstrap_95_ci"] = bootstrap_ci(
        rows, lambda sample: _coverage(sample, lambda row: row["three_way_agreement"]),
        seed=bootstrap_seed + 20, replicates=bootstrap_replicates,
    )
    three_way["reviewer_accuracy"]["bootstrap_95_ci"] = bootstrap_ci(
        rows, lambda sample: _accuracy([row for row in sample if row["three_way_agreement"]], "ai_lot"),
        seed=bootstrap_seed + 21, replicates=bootstrap_replicates,
    )

    two = [row for row in rows if row["algorithm_cota_agreement"]]
    two_correct = [row for row in two if row["algorithm_lot"] == row["reviewer_lot"]]
    two_wrong = [row for row in two if row["algorithm_lot"] != row["reviewer_lot"]]
    retained = [row for row in two if row["usable_ai_vote"] and row["ai_lot"] == row["algorithm_lot"]]

    def rejected_by_reason(cohort: list[dict[str, Any]], reason: str) -> int:
        return sum(not row["three_way_agreement"] and row["routing_reason"] == reason for row in cohort)

    silent_abstention = rejected_by_reason(two_wrong, "ai_abstained")
    silent_disagreement = rejected_by_reason(two_wrong, "three_way_disagreement")
    silent_invalid = rejected_by_reason(two_wrong, "invalid_ai_output")
    correct_abstention = rejected_by_reason(two_correct, "ai_abstained")
    correct_disagreement = rejected_by_reason(two_correct, "three_way_disagreement")
    correct_invalid = rejected_by_reason(two_correct, "invalid_ai_output")
    silent_prevented = silent_abstention + silent_disagreement + silent_invalid
    correct_rejected = correct_abstention + correct_disagreement + correct_invalid
    total_rejected = len(two) - len(retained)
    two_accuracy = len(two_correct) / len(two) if two else None
    three_accuracy = three_correct / len(three) if three else None
    precision_gain = three_accuracy - two_accuracy if three_accuracy is not None and two_accuracy is not None else None

    def accuracy_difference(sample: list[dict[str, Any]]) -> float | None:
        sample_two = [row for row in sample if row["algorithm_cota_agreement"]]
        sample_three = [row for row in sample if row["three_way_agreement"]]
        left = _accuracy(sample_three, "ai_lot")
        right = _accuracy(sample_two, "algorithm_lot")
        return left - right if left is not None and right is not None else None

    incremental = {
        "two_way_algorithm_cota": {
            "coverage": _metric(len(two), n),
            "reviewer_accuracy": _metric(len(two_correct), len(two)),
        },
        "three_way_algorithm_ai_cota": {
            "coverage": _metric(len(three), n),
            "reviewer_accuracy": _metric(three_correct, len(three)),
        },
        "two_way_silent_failures_rejected_by_ai": silent_prevented,
        "silent_failures_rejected_by_abstention": silent_abstention,
        "silent_failures_rejected_by_disagreement": silent_disagreement,
        "silent_failures_rejected_by_invalid_ai_output": silent_invalid,
        "previously_correct_two_way_agreements_rejected_by_ai": correct_rejected,
        "correct_two_way_candidates_rejected_by_abstention": correct_abstention,
        "correct_two_way_candidates_rejected_by_disagreement": correct_disagreement,
        "correct_two_way_candidates_rejected_by_invalid_ai_output": correct_invalid,
        "two_way_agreement_cases_retained": len(retained),
        "net_reduction_in_incorrect_auto_accept_candidates": len(two_wrong) - three_wrong,
        "coverage_lost_per_silent_failure_prevented": (
            round(correct_rejected / silent_prevented, 6) if silent_prevented else None
        ),
        "all_candidates_rejected_per_silent_failure_prevented": (
            round(total_rejected / silent_prevented, 6) if silent_prevented else None
        ),
        "precision_gain_rate_difference": round(precision_gain, 6) if precision_gain is not None else None,
        "three_way_silent_failures_remaining": three_wrong,
        "any_three_way_silent_failures_remain": three_wrong > 0,
        "consensus_accuracy_difference_bootstrap_95_ci": bootstrap_ci(
            rows, accuracy_difference,
            seed=bootstrap_seed + 30, replicates=bootstrap_replicates,
        ),
    }

    abstained_exact = sum(row["ai_lot"] == row["reviewer_lot"] for row in abstained_generated)
    abstained_within = sum(abs(row["ai_lot"] - row["reviewer_lot"]) <= 1 for row in abstained_generated)
    abstained_diagnostics = {
        "coverage_or_count": _metric(len(abstained_generated), n),
        "exact_accuracy": _metric(abstained_exact, len(abstained_generated)),
        "within_one_accuracy": _metric(abstained_within, len(abstained_generated)),
    }
    abstention = _metric(sum(row["abstained"] for row in rows), n)
    invalid_output = _metric(sum(row["invalid_ai_output"] for row in rows), n)
    policy_analysis = routing_policy_analysis(
        rows,
        bootstrap_seed=bootstrap_seed,
        bootstrap_replicates=bootstrap_replicates,
    )
    vote_patterns = _vote_pattern_metrics(rows)
    conditional_analysis = _conditional_agreement_analysis(rows)

    return {
        "eligible_patients": n,
        "generated_total": generated_metrics,
        "selective_ai": selective_metrics,
        "abstention": abstention,
        "invalid_ai_output": invalid_output,
        "abstained_generated_total": abstained_diagnostics,
        "ai_alone": {
            "deprecated": True,
            "semantic_label": "generated total diagnostics, irrespective of abstention",
            "deprecated_alias_for": "generated_total",
            **generated_metrics,
            "abstention": abstention,
        },
        "pairwise_agreement": pairwise,
        "three_way_consensus": three_way,
        "routing_policy_analysis": policy_analysis,
        "vote_pattern_analysis": vote_patterns,
        "conditional_agreement_analysis": conditional_analysis,
        "incremental_effect_of_ai": incremental,
        "comparison_baseline": {
            "description": "Existing algorithm+COTA agreement subset over all 136 adjudicated patients",
            "reviewer_accuracy": 0.828947,
            "coverage": 0.558824,
            "silent_failure_count": 13,
        },
        "uncertainty": {
            "method": "Patient-level nonparametric bootstrap of the eligible cohort with replacement, followed by subset filtering; percentile 95% intervals. Undefined empty-subset replicates are omitted.",
            "seed": bootstrap_seed,
            "requested_replicates": bootstrap_replicates,
            "small_subset_warning": "Intervals can be unstable for small agreement subsets; do not overstate point estimates.",
        },
    }
