from __future__ import annotations

from typing import Any

from config import HIGH_CONFIDENCE_THRESHOLD
from schemas import AdjudicationResult, AuditCritique, PipelineDecision


def mapping(result: AdjudicationResult) -> dict[int, int]:
    """Return vendor LoT -> adjudicated LoT numbering."""
    return {
        int(line.vendor_lot): int(line.adjudicated_lot)
        for line in result.corrected_lines
    }


def transition_flags(lot_mapping: dict[int, int]) -> dict[int, str]:
    """
    Convert absolute adjudicated LoT numbers into transition decisions.

    For each vendor LoT after the first:
      - "merge" means the corrected LoT stayed the same.
      - "split" means the corrected LoT increased.

    This removes harmless downstream numbering offsets. For example,
    [1, 2, 4, 5] and [1, 2, 3, 4] both make a split at the last two
    transitions, even though the absolute labels differ.
    """
    flags: dict[int, str] = {}
    previous_corrected: int | None = None

    for vendor_lot in sorted(lot_mapping):
        corrected_lot = lot_mapping[vendor_lot]

        if previous_corrected is None:
            flags[vendor_lot] = "start"
        elif corrected_lot == previous_corrected:
            flags[vendor_lot] = "merge"
        elif corrected_lot > previous_corrected:
            flags[vendor_lot] = "split"
        else:
            flags[vendor_lot] = "invalid"

        previous_corrected = corrected_lot

    return flags


def compare_results(
    worker: AdjudicationResult,
    auditor: AdjudicationResult,
) -> dict[str, Any]:
    """
    Compare Worker and Auditor by split/merge decisions, not absolute LoT
    numbers. This prevents downstream carryover numbering from creating false
    disagreements and unnecessary human-review cases.
    """
    differences: list[dict[str, Any]] = []

    if worker.patient_id != auditor.patient_id:
        differences.append(
            {
                "type": "patient_id",
                "worker": worker.patient_id,
                "auditor": auditor.patient_id,
                "major": True,
            }
        )

    worker_map = mapping(worker)
    auditor_map = mapping(auditor)
    worker_transitions = transition_flags(worker_map)
    auditor_transitions = transition_flags(auditor_map)

    all_vendor_lots = sorted(set(worker_transitions) | set(auditor_transitions))

    for vendor_lot in all_vendor_lots:
        worker_value = worker_transitions.get(vendor_lot, "missing")
        auditor_value = auditor_transitions.get(vendor_lot, "missing")

        if worker_value != auditor_value:
            differences.append(
                {
                    "type": "transition_decision",
                    "vendor_lot": vendor_lot,
                    "worker": worker_value,
                    "auditor": auditor_value,
                    "worker_absolute_lot": worker_map.get(vendor_lot),
                    "auditor_absolute_lot": auditor_map.get(vendor_lot),
                    "major": True,
                }
            )

    exact_transition_agreement = not differences

    return {
        # Kept for compatibility with pipeline.py and eval_pipeline.py.
        # It now means exact agreement on all split/merge decisions.
        "exact_mapping_agreement": exact_transition_agreement,
        "exact_transition_agreement": exact_transition_agreement,
        "major_disagreement": any(item["major"] for item in differences),
        "differences": differences,
        "worker_mapping": worker_map,
        "auditor_mapping": auditor_map,
        "worker_transitions": worker_transitions,
        "auditor_transitions": auditor_transitions,
        "absolute_numbering_differs_only": (
            exact_transition_agreement and worker_map != auditor_map
        ),
    }


def _has_review_trigger(result: AdjudicationResult) -> bool:
    """Return True when the result is not safe for automatic acceptance."""
    if result.overall_confidence < HIGH_CONFIDENCE_THRESHOLD:
        return True

    if result.requires_human_review:
        return True

    if not result.corrected_lines:
        return True

    return any(
        line.confidence < HIGH_CONFIDENCE_THRESHOLD
        or line.action == "uncertain"
        for line in result.corrected_lines
    )


def make_decision(
    worker: AdjudicationResult,
    auditor: AdjudicationResult,
    comparison: dict[str, Any],
    critique: AuditCritique | None = None,
) -> PipelineDecision:
    """
    Automatically accept only when Worker and Auditor agree on every
    split/merge transition and both outputs clear confidence and boundary-level
    safeguards. A disagreement always goes to human review.

    The critique argument remains for API compatibility but cannot override a
    split/merge disagreement.
    """
    del critique

    worker_trigger = _has_review_trigger(worker)
    auditor_trigger = _has_review_trigger(auditor)
    transitions_agree = bool(comparison["exact_mapping_agreement"])

    if transitions_agree and not worker_trigger and not auditor_trigger:
        numbering_note = (
            " Absolute LoT numbers differ only because of carryover numbering."
            if comparison.get("absolute_numbering_differs_only")
            else ""
        )

        return PipelineDecision(
            patient_id=worker.patient_id,
            status="accepted_agreement",
            selected_result="worker",
            reason=(
                "Worker and auditor produced the same split/merge decisions "
                "and both cleared confidence and boundary-level review checks."
                + numbering_note
            ),
            worker_confidence=worker.overall_confidence,
            auditor_confidence=auditor.overall_confidence,
            major_disagreement=False,
        )

    if transitions_agree:
        return PipelineDecision(
            patient_id=worker.patient_id,
            status="human_review",
            selected_result="none",
            reason=(
                "The split/merge decisions agree, but at least one agent "
                "reported low confidence, an uncertain line action, or an explicit review flag."
            ),
            worker_confidence=worker.overall_confidence,
            auditor_confidence=auditor.overall_confidence,
            major_disagreement=False,
        )

    return PipelineDecision(
        patient_id=worker.patient_id,
        status="human_review",
        selected_result="none",
        reason="Worker and auditor disagree on at least one split/merge decision.",
        worker_confidence=worker.overall_confidence,
        auditor_confidence=auditor.overall_confidence,
        major_disagreement=True,
    )
