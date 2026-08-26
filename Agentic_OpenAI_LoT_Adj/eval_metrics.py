from __future__ import annotations

from typing import Any


def mapping_from_lines(lines: list[dict[str, Any]]) -> dict[int, int]:
    return {
        int(line["vendor_lot"]): int(line["adjudicated_lot"])
        for line in lines
    }


def mapping_from_result(result: Any) -> dict[int, int]:
    return {
        int(line.vendor_lot): int(line.adjudicated_lot)
        for line in result.corrected_lines
    }


def transition_flags(mapping: dict[int, int]) -> dict[int, str]:
    """
    Convert absolute LoT numbers to the actual classification decision at each
    vendor line: start, split, or merge.

    This correctly ignores downstream carryover numbering differences.
    """
    flags: dict[int, str] = {}
    previous: int | None = None

    for vendor_lot in sorted(mapping):
        corrected = mapping[vendor_lot]

        if previous is None:
            flags[vendor_lot] = "start"
        elif corrected == previous:
            flags[vendor_lot] = "merge"
        elif corrected > previous:
            flags[vendor_lot] = "split"
        else:
            flags[vendor_lot] = "invalid"

        previous = corrected

    return flags


def merge_flags(mapping: dict[int, int]) -> dict[int, bool]:
    return {
        vendor_lot: decision == "merge"
        for vendor_lot, decision in transition_flags(mapping).items()
    }


def score_mapping(
    predicted: dict[int, int] | None,
    truth: dict[int, int],
) -> dict[str, Any]:
    """
    Score split/merge decisions rather than literal LoT-number equality.

    `exact_match` means every transition decision for the patient was correct.
    `correct_lines` includes the first LoT as a correct start decision.
    `correct_transitions` excludes the first LoT and is the pure split/merge
    classification score.
    """
    truth_decisions = transition_flags(truth)
    truth_lots = sorted(truth_decisions)
    transition_lots = truth_lots[1:]

    if predicted is None:
        true_merge_count = sum(
            truth_decisions.get(lot) == "merge" for lot in transition_lots
        )
        return {
            "exact_match": False,
            "correct_lines": 0,
            "total_lines": len(truth_lots),
            "line_accuracy": 0.0,
            "correct_transitions": 0,
            "total_transitions": len(transition_lots),
            "transition_accuracy": 0.0,
            "merge_tp": 0,
            "merge_fp": 0,
            "merge_fn": true_merge_count,
        }

    predicted_decisions = transition_flags(predicted)
    all_lots = sorted(set(truth_decisions) | set(predicted_decisions))

    correct_lines = sum(
        predicted_decisions.get(lot) == truth_decisions.get(lot)
        for lot in all_lots
    )

    comparison_transition_lots = [
        lot for lot in all_lots
        if truth_decisions.get(lot) != "start"
        and predicted_decisions.get(lot) != "start"
    ]

    correct_transitions = sum(
        predicted_decisions.get(lot) == truth_decisions.get(lot)
        for lot in comparison_transition_lots
    )

    exact = (
        set(predicted_decisions) == set(truth_decisions)
        and all(
            predicted_decisions[lot] == truth_decisions[lot]
            for lot in truth_decisions
        )
    )

    true_merge = merge_flags(truth)
    pred_merge = merge_flags(predicted)
    merge_eval_lots = sorted(set(true_merge) | set(pred_merge))[1:]

    merge_tp = sum(
        pred_merge.get(lot, False) and true_merge.get(lot, False)
        for lot in merge_eval_lots
    )
    merge_fp = sum(
        pred_merge.get(lot, False) and not true_merge.get(lot, False)
        for lot in merge_eval_lots
    )
    merge_fn = sum(
        not pred_merge.get(lot, False) and true_merge.get(lot, False)
        for lot in merge_eval_lots
    )

    return {
        "exact_match": exact,
        "correct_lines": correct_lines,
        "total_lines": len(all_lots),
        "line_accuracy": correct_lines / len(all_lots) if all_lots else 1.0,
        "correct_transitions": correct_transitions,
        "total_transitions": len(comparison_transition_lots),
        "transition_accuracy": (
            correct_transitions / len(comparison_transition_lots)
            if comparison_transition_lots
            else 1.0
        ),
        "merge_tp": merge_tp,
        "merge_fp": merge_fp,
        "merge_fn": merge_fn,
    }


def safe_divide(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def aggregate_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate both prediction-quality and patient-triage metrics.

    The primary triage question is whether patients whose original COTA
    split/merge sequence was wrong were sent to human review.
    """
    completed = [row for row in rows if not row.get("error")]
    total = len(rows)
    processed = len(completed)

    review_rows = [row for row in completed if row["status"] == "human_review"]
    auto_rows = [row for row in completed if row["status"] != "human_review"]

    # Patient-level triage confusion matrix:
    # TP: COTA error correctly sent to review
    # FN: COTA error incorrectly auto-accepted
    # FP: correct COTA patient unnecessarily reviewed
    # TN: correct COTA patient safely auto-accepted
    triage_tp = sum(
        bool(row["cota_patient_misclassified"])
        and row["status"] == "human_review"
        for row in completed
    )
    triage_fn = sum(
        bool(row["cota_patient_misclassified"])
        and row["status"] != "human_review"
        for row in completed
    )
    triage_fp = sum(
        not bool(row["cota_patient_misclassified"])
        and row["status"] == "human_review"
        for row in completed
    )
    triage_tn = sum(
        not bool(row["cota_patient_misclassified"])
        and row["status"] != "human_review"
        for row in completed
    )

    true_misclassified = triage_tp + triage_fn
    true_correct = triage_fp + triage_tn

    # A reviewed patient catches every erroneous split/merge transition within
    # that patient. This measures error-transition coverage, not prediction.
    total_cota_errors = sum(
        int(row.get("cota_misclassified_transition_count", 0))
        for row in completed
    )
    caught_cota_errors = sum(
        int(row.get("cota_misclassified_transition_count", 0))
        for row in review_rows
    )
    missed_cota_errors = total_cota_errors - caught_cota_errors

    auto_correct = sum(bool(row["selected_exact_match"]) for row in auto_rows)
    worker_exact = sum(bool(row["worker_exact_match"]) for row in completed)
    auditor_exact = sum(bool(row["auditor_exact_match"]) for row in completed)
    selected_exact = sum(bool(row["selected_exact_match"]) for row in completed)

    total_lines = sum(int(row["truth_line_count"]) for row in completed)
    worker_correct_lines = sum(int(row["worker_correct_lines"]) for row in completed)
    auditor_correct_lines = sum(int(row["auditor_correct_lines"]) for row in completed)
    selected_correct_lines = sum(int(row["selected_correct_lines"]) for row in completed)

    total_transitions = sum(
        int(row.get("truth_transition_count", max(int(row["truth_line_count"]) - 1, 0)))
        for row in completed
    )
    worker_correct_transitions = sum(
        int(row.get("worker_correct_transitions", 0)) for row in completed
    )
    auditor_correct_transitions = sum(
        int(row.get("auditor_correct_transitions", 0)) for row in completed
    )
    selected_correct_transitions = sum(
        int(row.get("selected_correct_transitions", 0)) for row in completed
    )

    merge_tp = sum(int(row["selected_merge_tp"]) for row in completed)
    merge_fp = sum(int(row["selected_merge_fp"]) for row in completed)
    merge_fn = sum(int(row["selected_merge_fn"]) for row in completed)
    merge_precision = safe_divide(merge_tp, merge_tp + merge_fp)
    merge_recall = safe_divide(merge_tp, merge_tp + merge_fn)

    result = {
        "sampled_patients": total,
        "successfully_processed": processed,
        "processing_errors": total - processed,

        # Primary patient-triage metrics.
        "true_cota_misclassified_patients": true_misclassified,
        "true_cota_correct_patients": true_correct,
        "misclassified_patients_caught_for_review": triage_tp,
        "misclassified_patients_missed_and_autoaccepted": triage_fn,
        "correct_patients_unnecessarily_reviewed": triage_fp,
        "correct_patients_safely_autoaccepted": triage_tn,
        "misclassification_catch_rate": safe_divide(triage_tp, true_misclassified),
        "missed_error_rate": safe_divide(triage_fn, true_misclassified),
        "review_precision": safe_divide(triage_tp, triage_tp + triage_fp),
        "unnecessary_review_rate_among_correct_patients": safe_divide(
            triage_fp, true_correct
        ),
        "correct_patient_autoaccept_rate": safe_divide(triage_tn, true_correct),
        "safe_autoaccept_rate": safe_divide(triage_tn, triage_tn + triage_fn),
        "overall_human_review_rate": safe_divide(len(review_rows), processed),
        "overall_autoaccept_rate": safe_divide(len(auto_rows), processed),

        # Transition-error coverage by the patient review queue.
        "true_cota_misclassified_transitions": total_cota_errors,
        "misclassified_transitions_caught_via_patient_review": caught_cota_errors,
        "misclassified_transitions_missed_via_autoaccept": missed_cota_errors,
        "misclassified_transition_catch_rate": safe_divide(
            caught_cota_errors, total_cota_errors
        ),

        # Secondary model-adjudication metrics.
        "worker_exact_patient_accuracy": safe_divide(worker_exact, processed),
        "auditor_exact_patient_accuracy": safe_divide(auditor_exact, processed),
        "selected_exact_patient_accuracy_all_cases": safe_divide(
            selected_exact, processed
        ),
        "worker_lot_classification_accuracy": safe_divide(
            worker_correct_lines, total_lines
        ),
        "auditor_lot_classification_accuracy": safe_divide(
            auditor_correct_lines, total_lines
        ),
        "selected_lot_classification_accuracy_all_cases": safe_divide(
            selected_correct_lines, total_lines
        ),
        "worker_auditor_agreement_rate": safe_divide(
            sum(bool(row["worker_auditor_agree"]) for row in completed), processed
        ),
        "automatic_accept_accuracy_against_doctor_labels": safe_divide(
            auto_correct, len(auto_rows)
        ),
        "automatic_accept_error_rate_against_doctor_labels": safe_divide(
            len(auto_rows) - auto_correct, len(auto_rows)
        ),
        "selected_merge_precision": merge_precision,
        "selected_merge_recall": merge_recall,
        "selected_merge_f1": safe_divide(
            2 * merge_precision * merge_recall,
            merge_precision + merge_recall,
        ),
    }

    if total_transitions:
        result.update(
            {
                "total_split_merge_transitions": total_transitions,
                "worker_split_merge_accuracy": safe_divide(
                    worker_correct_transitions, total_transitions
                ),
                "auditor_split_merge_accuracy": safe_divide(
                    auditor_correct_transitions, total_transitions
                ),
                "selected_split_merge_accuracy_all_cases": safe_divide(
                    selected_correct_transitions, total_transitions
                ),
            }
        )

    return result
