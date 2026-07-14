"""Strict data contracts for blind model inputs and outputs."""

from __future__ import annotations

from datetime import date
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)


class BlindEvent(StrictModel):
    order: int = Field(ge=1)
    drugs: list[str]

    @field_validator("drugs")
    @classmethod
    def unique_drugs(cls, value: list[str]) -> list[str]:
        if len(value) != len(set(value)):
            raise ValueError("event drugs must be unique")
        if any(not drug.strip() for drug in value):
            raise ValueError("drug names must be non-empty")
        return value


class BlindContext(StrictModel):
    diagnosis_date: str | None

    @field_validator("diagnosis_date")
    @classmethod
    def iso_date(cls, value: str | None) -> str | None:
        if value is not None:
            date.fromisoformat(value)
        return value


class BlindModelInput(StrictModel):
    schema_version: Literal["2.0.0"]
    case_key: str = Field(pattern=r"^CASE-[0-9a-f]{20}$")
    trajectory: list[BlindEvent]
    context: BlindContext

    @model_validator(mode="after")
    def sequential_orders(self) -> "BlindModelInput":
        if [event.order for event in self.trajectory] != list(range(1, len(self.trajectory) + 1)):
            raise ValueError("trajectory event orders must be consecutive starting at 1")
        return self


Decision = Literal["SAME_LINE", "NEW_LINE", "INSUFFICIENT_INFORMATION"]
EvidenceStrength = Literal["strong", "moderate", "limited", "conflicting"]


class TransitionDecision(StrictModel):
    transition_index: int = Field(ge=1)
    decision: Decision
    reason_codes: list[str]
    evidence_strength: EvidenceStrength
    retrieved_case_ids: list[str]

    @field_validator("reason_codes")
    @classmethod
    def valid_reason_codes(cls, value: list[str]) -> list[str]:
        if any(not code.strip() for code in value):
            raise ValueError("reason codes must be non-empty")
        return value

    @field_validator("retrieved_case_ids")
    @classmethod
    def valid_case_ids(cls, value: list[str]) -> list[str]:
        import re

        if len(value) != len(set(value)):
            raise ValueError("retrieved case IDs must be unique")
        if any(re.fullmatch(r"CASE-[0-9a-f]{20}", case_id) is None for case_id in value):
            raise ValueError("retrieved case IDs must be pseudonymous CASE identifiers")
        return value


class BlindLOTResponse(StrictModel):
    schema_version: Literal["1.0.0"]
    ai_lot_count: int = Field(ge=0)
    abstained: bool
    transitions: list[TransitionDecision]
    warnings: list[str]

    @model_validator(mode="after")
    def internally_consistent(self) -> "BlindLOTResponse":
        indices = [item.transition_index for item in self.transitions]
        if indices != list(range(1, len(self.transitions) + 1)):
            raise ValueError("transition indices must be consecutive starting at 1")
        if self.transitions:
            expected = 1 + sum(item.decision == "NEW_LINE" for item in self.transitions)
            if self.ai_lot_count != expected:
                raise ValueError("ai_lot_count must equal 1 plus NEW_LINE transitions")
        if any(item.decision == "INSUFFICIENT_INFORMATION" for item in self.transitions):
            if not self.abstained:
                raise ValueError("INSUFFICIENT_INFORMATION requires abstained=true")
        return self


def validate_model_response(payload: str | dict[str, Any], regimen_event_count: int) -> BlindLOTResponse:
    """Reject malformed or logically inconsistent output without repairing it."""
    if isinstance(payload, str):
        response = BlindLOTResponse.model_validate_json(payload)
    else:
        response = BlindLOTResponse.model_validate(payload)
    expected_transitions = max(0, regimen_event_count - 1)
    if len(response.transitions) != expected_transitions:
        raise ValueError(
            f"expected {expected_transitions} transitions for {regimen_event_count} events; "
            f"received {len(response.transitions)}"
        )
    if regimen_event_count == 0:
        if response.ai_lot_count != 0 or not response.abstained:
            raise ValueError("empty trajectories require ai_lot_count=0 and abstained=true")
    elif not response.transitions and response.ai_lot_count != 1:
        raise ValueError("a one-event trajectory requires ai_lot_count=1")
    return response


def response_json_schema() -> dict[str, Any]:
    """JSON schema sent to providers that support constrained output."""
    return BlindLOTResponse.model_json_schema()
