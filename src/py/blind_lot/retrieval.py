"""Deterministic, local retrieval over order-only trajectory features."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable


TRANSPLANT_MARKERS = ("transplant", "sct", "asct", "stem cell", "stem_cell")
CART_MARKERS = ("cart", "car-t", "idecabtagene", "ciltacabtagene", "abecma", "carvykti")
BISPECIFIC_MARKERS = (
    "teclistamab", "elranatamab", "talquetamab", "linvoseltamab",
    "tecvayli", "elrexfio", "talvey",
)
UNKNOWN_MARKERS = ("investigational", "study drug", "trial", "unknown", "unmapped", "protocol")
STEROIDS = {"dexamethasone", "prednisone", "prednisolone", "methylprednisolone", "dex"}


@dataclass(frozen=True)
class RetrievalCandidate:
    case_key: str
    patient_key: str
    fold: int
    exclusion_group: str
    trajectory: list[dict[str, Any]]
    reviewer_consensus_lot: int


@dataclass(frozen=True)
class RetrievalHit:
    candidate: RetrievalCandidate
    score: float
    components: dict[str, float]

    def debug_dict(self) -> dict[str, Any]:
        return {
            "case_key": self.candidate.case_key,
            "fold": self.candidate.fold,
            "score": self.score,
            "component_scores": self.components,
        }


def _sets(trajectory: list[dict[str, Any]]) -> list[frozenset[str]]:
    return [frozenset(str(drug).lower() for drug in event["drugs"]) for event in trajectory]


def transition_types(trajectory: list[dict[str, Any]]) -> list[str]:
    result: list[str] = []
    sets = _sets(trajectory)
    for prior, current in zip(sets, sets[1:]):
        if current == prior:
            kind = "IDENTICAL"
        elif prior < current:
            kind = "ADDITION"
        elif current < prior:
            kind = "REMOVAL"
        elif prior.isdisjoint(current):
            kind = "NONOVERLAP_SWITCH"
        else:
            kind = "OVERLAPPING_CHANGE"
        result.append(kind)
    return result


def _contains(trajectory: list[dict[str, Any]], markers: tuple[str, ...]) -> bool:
    return any(marker in drug.lower() for event in trajectory for drug in event["drugs"] for marker in markers)


def _steroid_only(trajectory: list[dict[str, Any]]) -> bool:
    return any(event["drugs"] and {d.lower() for d in event["drugs"]} <= STEROIDS for event in trajectory)


def _reappearance(trajectory: list[dict[str, Any]]) -> bool:
    sets = _sets(trajectory)
    return any(sets[index] in sets[: index - 1] for index in range(2, len(sets)))


def trajectory_indicators(trajectory: list[dict[str, Any]]) -> dict[str, bool]:
    return {
        "transplant_indicator": _contains(trajectory, TRANSPLANT_MARKERS),
        "cart_indicator": _contains(trajectory, CART_MARKERS),
        "bispecific_indicator": _contains(trajectory, BISPECIFIC_MARKERS),
        "unknown_indicator": _contains(trajectory, UNKNOWN_MARKERS),
        "steroid_only_indicator": _steroid_only(trajectory),
        "reappearance_indicator": _reappearance(trajectory),
    }


def _length_similarity(left: int, right: int) -> float:
    return 1.0 if left == right == 0 else 1.0 - abs(left - right) / max(left, right, 1)


def _sequence_similarity(left: list[str], right: list[str]) -> float:
    if not left and not right:
        return 1.0
    width = max(len(left), len(right))
    return sum(a == b for a, b in zip(left, right)) / width


def _regimen_similarity(left: list[dict[str, Any]], right: list[dict[str, Any]]) -> float:
    left_sets, right_sets = _sets(left), _sets(right)
    if not left_sets and not right_sets:
        return 1.0
    width = max(len(left_sets), len(right_sets))
    total = 0.0
    for index in range(width):
        a = left_sets[index] if index < len(left_sets) else frozenset()
        b = right_sets[index] if index < len(right_sets) else frozenset()
        union = a | b
        total += len(a & b) / len(union) if union else 1.0
    return total / width


def similarity(left: list[dict[str, Any]], right: list[dict[str, Any]]) -> tuple[float, dict[str, float]]:
    left_flags, right_flags = trajectory_indicators(left), trajectory_indicators(right)
    components = {
        "trajectory_length": _length_similarity(len(left), len(right)),
        "transition_type_sequence": _sequence_similarity(transition_types(left), transition_types(right)),
        "regimen_set_jaccard": _regimen_similarity(left, right),
    }
    for name in left_flags:
        components[name] = float(left_flags[name] == right_flags[name])
    score = (
        0.15 * components["trajectory_length"]
        + 0.20 * components["transition_type_sequence"]
        + 0.35 * components["regimen_set_jaccard"]
        + 0.05 * sum(components[name] for name in left_flags)
    )
    return round(score, 12), {key: round(value, 12) for key, value in components.items()}


class LocalRetriever:
    def __init__(self, candidates: Iterable[RetrievalCandidate]) -> None:
        self.candidates = tuple(candidates)

    def retrieve(
        self,
        target_trajectory: list[dict[str, Any]],
        *,
        target_patient_key: str,
        target_fold: int,
        target_exclusion_group: str,
        k: int,
    ) -> list[RetrievalHit]:
        if k not in {0, 3, 5}:
            raise ValueError("retrieval k must be one of 0, 3, or 5")
        eligible = [
            candidate for candidate in self.candidates
            if candidate.patient_key != target_patient_key
            and candidate.fold != target_fold
            and candidate.exclusion_group != target_exclusion_group
        ]
        hits = []
        for candidate in eligible:
            score, components = similarity(target_trajectory, candidate.trajectory)
            hits.append(RetrievalHit(candidate, score, components))
        hits.sort(key=lambda hit: (-hit.score, hit.candidate.case_key))
        return hits[:k]


def assert_retrieval_safe(
    hits: list[RetrievalHit], *, target_patient_key: str, target_fold: int, target_exclusion_group: str
) -> None:
    for hit in hits:
        candidate = hit.candidate
        if candidate.patient_key == target_patient_key:
            raise ValueError("retrieval leaked the target patient")
        if candidate.fold == target_fold:
            raise ValueError("retrieval leaked a held-out fold record")
        if candidate.exclusion_group == target_exclusion_group:
            raise ValueError("retrieval leaked the target exclusion group")
