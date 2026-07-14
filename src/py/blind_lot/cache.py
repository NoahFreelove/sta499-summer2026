"""Validated disk cache with configuration-complete keys."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from .models import BlindLOTResponse, validate_model_response


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def cache_key(
    *,
    model: str,
    reasoning_effort: str | None,
    prompt_version: str,
    knowledge_version: str,
    retrieval_k: int,
    retrieved_example_ids: list[str],
    blind_patient_input_hash: str,
    temperature: float | None = None,
) -> str:
    return canonical_hash({
        "model": model,
        "reasoning_effort": reasoning_effort,
        "temperature": temperature,
        "prompt_version": prompt_version,
        "knowledge_version": knowledge_version,
        "retrieval": {"kind": "deterministic_weighted_local_v1", "k": retrieval_k},
        "retrieved_example_ids": retrieved_example_ids,
        "blind_patient_input_hash": blind_patient_input_hash,
    })


class DiskCache:
    def __init__(self, directory: Path) -> None:
        self.directory = directory

    def get(self, key: str, regimen_event_count: int) -> BlindLOTResponse | None:
        path = self.directory / key[:2] / f"{key}.json"
        if not path.exists():
            return None
        return validate_model_response(json.loads(path.read_text(encoding="utf-8")), regimen_event_count)

    def put(self, key: str, response: BlindLOTResponse) -> Path:
        path = self.directory / key[:2] / f"{key}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            path.write_text(
                json.dumps(response.model_dump(mode="json"), sort_keys=True, separators=(",", ":")) + "\n",
                encoding="utf-8",
            )
        return path
