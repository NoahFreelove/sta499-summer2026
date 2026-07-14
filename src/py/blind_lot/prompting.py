"""Stable-prefix prompt construction and provenance-aware leakage checks."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .models import BlindModelInput
from .retrieval import RetrievalHit


FORBIDDEN_TARGET_FIELDS = {
    "cota_lot",
    "reviewer_lot",
    "reviewer_consensus_lot",
    "source_line_number",
    "source_group_index",
    "source_line_start_date",
    "source_line_end_date",
    "algorithm_lot",
    "algorithm_flags",
    "exclusion_group",
    "patient_key",
    "raw_patient_id",
    "cpid",
}


@dataclass(frozen=True)
class PromptBundle:
    prompt_version: str
    knowledge_version: str
    stable_prefix: str
    patient_prompt: str
    target_payload: dict[str, Any]
    retrieved_example_ids: list[str]
    rendered_retrieval_demonstrations: list[str]
    rendered_retrieval_context: str


def _package_root() -> Path:
    return Path(__file__).resolve().parent


def load_prompt_template() -> tuple[str, str]:
    text = (_package_root() / "prompts" / "order_only_v1.txt").read_text(encoding="utf-8")
    first, _, rest = text.partition("\n")
    return first.split(":", 1)[1].strip(), rest.strip()


def load_knowledge() -> tuple[str, dict[str, Any]]:
    value = json.loads((_package_root() / "knowledge" / "order_only_v1.json").read_text(encoding="utf-8"))
    return value["knowledge_version"], value


def _all_keys(value: Any) -> set[str]:
    if isinstance(value, dict):
        return set(value) | set().union(*(_all_keys(item) for item in value.values()), set())
    if isinstance(value, list):
        return set().union(*(_all_keys(item) for item in value), set())
    return set()


def validate_target_provenance(target: dict[str, Any]) -> BlindModelInput:
    """Validate the source record and explicitly whitelist promptable target fields."""
    model = BlindModelInput.model_validate(target)
    overlap = _all_keys(target) & FORBIDDEN_TARGET_FIELDS
    if overlap:
        raise ValueError(f"forbidden target fields present before prompting: {sorted(overlap)}")
    return model


def build_prompt(target: dict[str, Any], hits: list[RetrievalHit]) -> PromptBundle:
    validated = validate_target_provenance(target)
    prompt_version, template = load_prompt_template()
    knowledge_version, knowledge = load_knowledge()

    # Deliberately reconstruct a minimal payload instead of filtering a restricted record.
    target_payload = {
        "trajectory": [event.model_dump(mode="json") for event in validated.trajectory],
        "number_of_regimen_events": len(validated.trajectory),
    }
    examples = [
        {
            "case_id": hit.candidate.case_key,
            "blind_trajectory": hit.candidate.trajectory,
            "reviewer_consensus_patient_total_lot": hit.candidate.reviewer_consensus_lot,
        }
        for hit in hits
    ]
    rendered_examples = [
        json.dumps(example, sort_keys=True, separators=(",", ":")) for example in examples
    ]
    rendered_retrieval_context = json.dumps(examples, sort_keys=True, separators=(",", ":"))
    stable_prefix = template + "\n\nCOMPACT ORDER-ONLY KNOWLEDGE BASE:\n" + json.dumps(
        knowledge, sort_keys=True, separators=(",", ":")
    )
    patient_value = {
        "benchmark_mode": "order-only blind AI benchmark",
        "target": target_payload,
        "permitted_retrieved_training_examples": examples,
    }
    patient_prompt = json.dumps(patient_value, sort_keys=True, separators=(",", ":"))
    assert_prompt_safe(patient_value)
    return PromptBundle(
        prompt_version=prompt_version,
        knowledge_version=knowledge_version,
        stable_prefix=stable_prefix,
        patient_prompt=patient_prompt,
        target_payload=target_payload,
        retrieved_example_ids=[hit.candidate.case_key for hit in hits],
        rendered_retrieval_demonstrations=rendered_examples,
        rendered_retrieval_context=rendered_retrieval_context,
    )


def assert_prompt_safe(patient_value: dict[str, Any]) -> None:
    """Check target provenance; training totals are allowed only inside examples."""
    if set(patient_value) != {
        "benchmark_mode", "target", "permitted_retrieved_training_examples"
    }:
        raise ValueError("unexpected patient prompt sections")
    target = patient_value["target"]
    if set(target) != {"trajectory", "number_of_regimen_events"}:
        raise ValueError("target prompt payload does not match the explicit allowlist")
    overlap = _all_keys(target) & FORBIDDEN_TARGET_FIELDS
    if overlap:
        raise ValueError(f"target-derived leakage in prompt: {sorted(overlap)}")
    for example in patient_value["permitted_retrieved_training_examples"]:
        if set(example) != {
            "case_id", "blind_trajectory", "reviewer_consensus_patient_total_lot"
        }:
            raise ValueError("retrieved example contains an unapproved field")
