"""Evaluation-only reaggregation of completed blind LOT benchmark runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from lot_data import assert_aggregate_only_report, read_jsonl

from .evaluation import EVALUATION_VERSION, evaluate_joined, normalize_joined_record


REAGGREGATED_JOINED_NAME = "joined_evaluation.reaggregated.jsonl"
REAGGREGATED_AGGREGATE_NAME = "aggregate_evaluation.reaggregated.json"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in sorted(rows, key=lambda item: item["case_key"]):
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def _public_report(metadata: dict[str, Any], metrics: dict[str, Any]) -> dict[str, Any]:
    run_id = metadata["run_id"]
    report = {
        "schema_version": "1.3.0",
        "report_scope": "aggregate_only",
        "benchmark_mode": "order-only blind AI benchmark",
        "ground_truth": "reviewer consensus",
        "cota_role": "independent vote, not a training label",
        "run_id_hash": hashlib.sha256(run_id.encode()).hexdigest(),
        "experiment": {
            "provider": metadata.get("provider"),
            "model": metadata.get("model"),
            "reasoning_effort": metadata.get("reasoning_effort"),
            "prompt_version": metadata.get("prompt_version"),
            "knowledge_version": metadata.get("knowledge_version"),
            "retrieval_k": metadata.get("retrieval_k"),
            "folds": metadata.get("folds"),
            "fold_count_evaluated": len(metadata.get("folds", [])),
            "limited_run": metadata.get("limit") is not None,
            "temperature": metadata.get("temperature"),
            "bootstrap_seed": metadata.get("random_seed"),
            "evaluation_version": EVALUATION_VERSION,
            "metric_semantics": {
                "abstained_outputs_are_votes": False,
                "abstained_numeric_totals_retained_for_diagnostics": True,
                "three_way_requires_non_abstained_ai": True,
                "primary_policy": "all-three agreement with usable AI",
                "pairwise_and_single_source_policies_are_exploratory": True,
            },
        },
        "metrics": metrics,
    }
    assert_aggregate_only_report(report)
    return report


def reevaluate_run(run_dir: Path) -> dict[str, Any]:
    """Reaggregate a completed run without importing providers or changing predictions."""
    run_dir = run_dir.resolve()
    restricted = run_dir / "restricted"
    public = run_dir / "public"
    joined_path = restricted / "joined_evaluation.jsonl"
    metadata_path = restricted / "experiment_metadata.json"
    prediction_path = restricted / "blind_predictions.jsonl"
    retrieval_debug_path = restricted / "retrieval_debug.jsonl"
    for path in (joined_path, metadata_path):
        if not path.exists():
            raise FileNotFoundError(path)

    original_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    source_joined_sha256 = _sha256_file(joined_path)
    prediction_sha256 = _sha256_file(prediction_path) if prediction_path.exists() else None
    retrieval_debug_sha256 = (
        _sha256_file(retrieval_debug_path) if retrieval_debug_path.exists() else None
    )
    rows = [normalize_joined_record(row) for row in read_jsonl(joined_path)]
    seed = int(original_metadata.get("random_seed", 4992026))
    replicates = int(original_metadata.get("bootstrap_replicates", 2000))
    metrics = evaluate_joined(rows, bootstrap_seed=seed, bootstrap_replicates=replicates)

    joined_output = restricted / REAGGREGATED_JOINED_NAME
    aggregate_output = public / REAGGREGATED_AGGREGATE_NAME
    provenance_output = restricted / "evaluation_reaggregation.json"
    _write_jsonl(joined_output, rows)
    _write_json(aggregate_output, _public_report(original_metadata, metrics))
    provenance = {
        "schema_version": "1.0.0",
        "evaluation_version": EVALUATION_VERSION,
        "evaluation_only": True,
        "provider_initialized": False,
        "model_requests_made": 0,
        "retrieval_executed": False,
        "predictions_modified": False,
        "source_experiment_metadata_sha256": _sha256_file(metadata_path),
        "source_joined_evaluation_sha256": source_joined_sha256,
        "blind_predictions_sha256": prediction_sha256,
        "retrieval_debug_sha256": retrieval_debug_sha256,
        "input_artifact_sha256": original_metadata.get("input_artifact_sha256", {}),
        "outputs": {
            "joined_evaluation": REAGGREGATED_JOINED_NAME,
            "aggregate_evaluation": REAGGREGATED_AGGREGATE_NAME,
        },
    }
    _write_json(provenance_output, provenance)
    if prediction_path.exists() and _sha256_file(prediction_path) != prediction_sha256:
        raise RuntimeError("blind prediction artifact changed during evaluation-only reaggregation")
    if retrieval_debug_path.exists() and _sha256_file(retrieval_debug_path) != retrieval_debug_sha256:
        raise RuntimeError("retrieval debug artifact changed during evaluation-only reaggregation")
    if _sha256_file(metadata_path) != provenance["source_experiment_metadata_sha256"]:
        raise RuntimeError("source experiment metadata changed during reaggregation")
    return {
        "run_dir": str(run_dir),
        "joined_output": str(joined_output),
        "aggregate_output": str(aggregate_output),
        "provenance_output": str(provenance_output),
        "metrics": metrics,
    }
