"""Provider-neutral inference interface with mock and OpenAI Responses implementations."""

from __future__ import annotations

import json
import os
import socket
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from typing import Protocol

from pydantic import ValidationError

from .models import BlindLOTResponse, response_json_schema, validate_model_response
from .retrieval import UNKNOWN_MARKERS, STEROIDS


class RetriableProviderError(RuntimeError):
    pass


class Provider(Protocol):
    def complete(self, stable_prefix: str, patient_prompt: str) -> str: ...


@dataclass(frozen=True)
class ProviderConfig:
    model: str
    reasoning_effort: str | None = None
    temperature: float | None = None
    request_timeout: float = 120.0
    retry_count: int = 2


class MockProvider:
    """Deterministic order-only classifier used only for tests and plumbing checks."""

    def __init__(self, config: ProviderConfig) -> None:
        self.config = config

    def complete(self, stable_prefix: str, patient_prompt: str) -> str:
        value = json.loads(patient_prompt)
        trajectory = value["target"]["trajectory"]
        retrieved = [item["case_id"] for item in value["permitted_retrieved_training_examples"]]
        transitions = []
        abstained = not trajectory
        for index, (before, after) in enumerate(zip(trajectory, trajectory[1:]), start=1):
            left, right = set(before["drugs"]), set(after["drugs"])
            text = " ".join(left | right).lower()
            steroid_only = bool(left and left <= STEROIDS) or bool(right and right <= STEROIDS)
            unknown = any(marker in text for marker in UNKNOWN_MARKERS)
            if unknown or steroid_only:
                decision, reason, strength = (
                    "INSUFFICIENT_INFORMATION",
                    "UNMAPPABLE_OR_INTENT_DEPENDENT",
                    "conflicting",
                )
                abstained = True
            elif left == right:
                decision, reason, strength = "SAME_LINE", "IDENTICAL_REGIMEN", "moderate"
            elif right < left:
                decision, reason, strength = "SAME_LINE", "PARTIAL_DRUG_REMOVAL", "limited"
            elif left.isdisjoint(right):
                decision, reason, strength = "NEW_LINE", "FULL_NONOVERLAPPING_SWITCH", "strong"
            elif left < right:
                decision, reason, strength = "NEW_LINE", "DRUG_ADDITION_ORDER_ONLY", "limited"
            else:
                decision, reason, strength = "NEW_LINE", "OVERLAPPING_REGIMEN_CHANGE", "moderate"
            transitions.append({
                "transition_index": index,
                "decision": decision,
                "reason_codes": [reason],
                "evidence_strength": strength,
                "retrieved_case_ids": retrieved,
            })
        count = 0 if not trajectory else 1 + sum(t["decision"] == "NEW_LINE" for t in transitions)
        return json.dumps({
            "schema_version": "1.0.0",
            "ai_lot_count": count,
            "abstained": abstained,
            "transitions": transitions,
            "warnings": ["ORDER_ONLY_TIMING_UNAVAILABLE"],
        }, sort_keys=True)


class OpenAIResponsesProvider:
    """Minimal OpenAI Responses API client using OPENAI_API_KEY."""

    endpoint = "https://api.openai.com/v1/responses"

    def __init__(self, config: ProviderConfig) -> None:
        self.config = config
        self.api_key = os.environ.get("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError("OPENAI_API_KEY is required for --provider openai")

    def complete(self, stable_prefix: str, patient_prompt: str) -> str:
        request_body: dict[str, object] = {
            "model": self.config.model,
            "instructions": stable_prefix,
            "input": patient_prompt,
            "store": False,
            "text": {
                "format": {
                    "type": "json_schema",
                    "name": "blind_lot_response",
                    "strict": True,
                    "schema": response_json_schema(),
                }
            },
        }
        if self.config.reasoning_effort:
            request_body["reasoning"] = {"effort": self.config.reasoning_effort}
        if self.config.temperature is not None:
            request_body["temperature"] = self.config.temperature
        request = urllib.request.Request(
            self.endpoint,
            data=json.dumps(request_body).encode("utf-8"),
            headers={
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=self.config.request_timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            body = error.read().decode("utf-8", errors="replace")
            if error.code == 429 or 500 <= error.code < 600:
                raise RetriableProviderError(f"OpenAI HTTP {error.code}: {body[:500]}") from error
            raise RuntimeError(f"OpenAI HTTP {error.code}: {body[:500]}") from error
        except (urllib.error.URLError, TimeoutError, socket.timeout) as error:
            raise RetriableProviderError(f"OpenAI transport failure: {error}") from error
        return _extract_output_text(payload)


def _extract_output_text(payload: dict[str, object]) -> str:
    if isinstance(payload.get("output_text"), str):
        return str(payload["output_text"])
    texts: list[str] = []
    for output in payload.get("output", []):  # type: ignore[union-attr]
        if not isinstance(output, dict):
            continue
        for content in output.get("content", []):
            if isinstance(content, dict) and content.get("type") == "output_text":
                texts.append(str(content.get("text", "")))
    if not texts:
        raise RetriableProviderError("OpenAI response did not contain output text")
    return "".join(texts)


def complete_with_retries(
    provider: Provider,
    stable_prefix: str,
    patient_prompt: str,
    *,
    regimen_event_count: int,
    retry_count: int,
    permitted_retrieved_case_ids: set[str] | None = None,
) -> tuple[BlindLOTResponse, int]:
    """Retry only transport/rate-limit and invalid JSON/schema failures."""
    attempts = 0
    while True:
        attempts += 1
        try:
            raw = provider.complete(stable_prefix, patient_prompt)
            response = validate_model_response(raw, regimen_event_count)
            if permitted_retrieved_case_ids is not None:
                cited = {
                    case_id
                    for transition in response.transitions
                    for case_id in transition.retrieved_case_ids
                }
                if not cited <= permitted_retrieved_case_ids:
                    raise ValueError("model cited a retrieved case that was not supplied")
            return response, attempts
        except (RetriableProviderError, ValidationError, json.JSONDecodeError, ValueError):
            if attempts > retry_count:
                raise
            time.sleep(min(2 ** (attempts - 1), 4))
