"""Strict evaluation-only pooling of explicitly frozen blind LOT runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from lot_data import assert_aggregate_only_report, read_jsonl

from .evaluation import EVALUATION_VERSION, evaluate_joined, normalize_joined_record
from .reaggregation import REAGGREGATED_AGGREGATE_NAME, REAGGREGATED_JOINED_NAME


ARTIFACTS = {
    "metadata": "experiment_metadata.json",
    "predictions": "blind_predictions.jsonl",
    "retrieval_diagnostics": "retrieval_debug.jsonl",
    "joined_records": "joined_evaluation.jsonl",
}
FROZEN_FIELDS = (
    "provider", "model", "reasoning_effort", "temperature", "prompt_version",
    "knowledge_version", "input_artifact_sha256", "random_seed", "bootstrap_replicates",
    "retrieval_training_folds", "run_purpose", "repeat_index",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in sorted(rows, key=lambda item: item["case_key"]):
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def _validate_manifest(manifest: dict[str, Any]) -> dict[tuple[int, int], str]:
    if manifest.get("schema_version") != "blind-lot-pooling-manifest-1.0.0":
        raise ValueError("unsupported pooling manifest schema_version")
    folds = [int(value) for value in manifest["expected_folds"]]
    if len(folds) != len(set(folds)):
        raise ValueError("expected_folds contains duplicates")
    entries: dict[tuple[int, int], str] = {}
    for source in manifest["source_runs"]:
        pair = (int(source["retrieval_k"]), int(source["fold"]))
        if pair in entries:
            raise ValueError(f"duplicate source run for retrieval_k/fold {pair}")
        entries[pair] = source["run_id"]
    expected = {(k, fold) for k in (3, 5) for fold in folds}
    missing = sorted(expected - set(entries))
    extra = sorted(set(entries) - expected)
    if missing or extra:
        raise ValueError(f"source run pairs do not match manifest; missing={missing}, extra={extra}")
    return entries


def _validate_expected_metadata(metadata: dict[str, Any], manifest: dict[str, Any]) -> None:
    expected = {
        "provider": manifest["provider"],
        "model": manifest["model"],
        "reasoning_effort": manifest["reasoning_effort"],
        "temperature": manifest.get("temperature"),
        "prompt_version": manifest["prompt_version"],
        "knowledge_version": manifest["knowledge_version"],
        "random_seed": manifest["bootstrap_seed"],
        "bootstrap_replicates": manifest["bootstrap_replicates"],
    }
    mismatches = [name for name, value in expected.items() if metadata.get(name) != value]
    if mismatches:
        raise ValueError(f"source metadata does not match manifest: {mismatches}")


def _assert_frozen(reference: dict[str, Any], metadata: dict[str, Any]) -> None:
    mismatches = [field for field in FROZEN_FIELDS if metadata.get(field) != reference.get(field)]
    if mismatches:
        raise ValueError(f"source runs differ on frozen metadata: {mismatches}")


def _public_report(metadata: dict[str, Any], metrics: dict[str, Any]) -> dict[str, Any]:
    report = {
        "schema_version": "1.3.0",
        "report_scope": "aggregate_only",
        "benchmark_mode": "order-only blind AI benchmark",
        "ground_truth": "reviewer consensus",
        "cota_role": "independent vote, not a training label",
        "experiment": {
            "provider": metadata["provider"],
            "model": metadata["model"],
            "reasoning_effort": metadata["reasoning_effort"],
            "temperature": metadata.get("temperature"),
            "prompt_version": metadata["prompt_version"],
            "knowledge_version": metadata["knowledge_version"],
            "retrieval_k": metadata["retrieval_k"],
            "folds": metadata["folds"],
            "bootstrap_seed": metadata["random_seed"],
            "bootstrap_replicates": metadata["bootstrap_replicates"],
            "evaluation_version": EVALUATION_VERSION,
            "evaluation_only": True,
        },
        "metrics": metrics,
    }
    assert_aggregate_only_report(report)
    return report


def _regression_values(metrics: dict[str, Any]) -> dict[str, Any]:
    policies = metrics["routing_policy_analysis"]

    def triple(name: str) -> list[int]:
        return [
            policies[name]["accepted_count"], policies[name]["correct_accepted_count"],
            policies[name]["incorrect_accepted_count"],
        ]
    return {
        "generated_correct": metrics["generated_total"]["exact_accuracy"]["count"],
        "usable": metrics["selective_ai"]["coverage"]["count"],
        "usable_correct": metrics["selective_ai"]["exact_accuracy"]["count"],
        "algorithm_ai": triple("algorithm_ai_agreement"),
        "cota_ai": triple("cota_ai_agreement"),
        "three_way": triple("three_way_agreement"),
    }


def pool_runs(manifest_path: Path, runs_root: Path, output: Path) -> dict[int, Path]:
    """Pool only manifest-listed completed runs, then evaluate saved joined rows."""
    manifest = _json(manifest_path)
    entries = _validate_manifest(manifest)
    expected_count = int(manifest["expected_patient_count"])
    reference: dict[str, Any] | None = None
    all_rows: dict[int, list[dict[str, Any]]] = {3: [], 5: []}
    provenance: dict[int, list[dict[str, Any]]] = {3: [], 5: []}
    before: dict[Path, str] = {}

    for (retrieval_k, fold), run_id in sorted(entries.items()):
        run_dir = runs_root / run_id
        restricted = run_dir / "restricted"
        paths = {name: restricted / filename for name, filename in ARTIFACTS.items()}
        for path in paths.values():
            if not path.is_file():
                raise FileNotFoundError(path)
        hashes = {name: sha256_file(path) for name, path in paths.items()}
        before.update({path: hashes[name] for name, path in paths.items()})
        metadata = _json(paths["metadata"])
        if metadata.get("run_id") != run_id:
            raise ValueError(f"run ID mismatch for {run_id}")
        if metadata.get("limit") is not None:
            raise ValueError(f"limited source run rejected: {run_id}")
        if metadata.get("dry_run") is True or metadata.get("provider") == "dry-run":
            raise ValueError(f"dry source run rejected: {run_id}")
        if metadata.get("retrieval_k") != retrieval_k or metadata.get("folds") != [fold]:
            raise ValueError(f"source run k/fold mismatch: {run_id}")
        _validate_expected_metadata(metadata, manifest)
        source_evaluation = metadata.get("evaluation_version")
        migration_compatible = fold == 0 and source_evaluation == "blind-lot-evaluation-1.1.0"
        if source_evaluation != manifest["evaluation_version"] and not migration_compatible:
            raise ValueError(f"source evaluation version mismatch: {run_id}")
        if reference is None:
            reference = metadata
        else:
            _assert_frozen(reference, metadata)
        rows = read_jsonl(paths["joined_records"])
        predictions = read_jsonl(paths["predictions"])
        diagnostics = read_jsonl(paths["retrieval_diagnostics"])
        if len(rows) != len(predictions) or len(rows) != len(diagnostics):
            raise ValueError(f"source artifact row counts disagree: {run_id}")
        if any(row.get("fold") != fold for row in rows):
            raise ValueError(f"joined row has incorrect fold in {run_id}")
        all_rows[retrieval_k].extend(normalize_joined_record(row) for row in rows)
        provenance[retrieval_k].append({
            "run_id": run_id, "retrieval_k": retrieval_k, "fold": fold,
            "row_count": len(rows), "artifact_sha256": hashes,
            "source_evaluation_version": metadata.get("evaluation_version"),
        })

    if manifest["evaluation_version"] != EVALUATION_VERSION:
        raise ValueError("manifest evaluation_version does not match installed evaluator")
    key_sets: dict[int, set[str]] = {}
    for retrieval_k, rows in all_rows.items():
        keys = [row["case_key"] for row in rows]
        if len(keys) != len(set(keys)):
            raise ValueError(f"duplicate case key in pooled k={retrieval_k} condition")
        if len(set(keys)) != expected_count:
            raise ValueError(f"k={retrieval_k} has {len(set(keys))}, expected {expected_count} patients")
        key_sets[retrieval_k] = set(keys)
    if key_sets[3] != key_sets[5]:
        raise ValueError("k=3 and k=5 pooled cohorts have different case-key sets")

    outputs: dict[int, Path] = {}
    for retrieval_k in (3, 5):
        run_dir = output / f"k{retrieval_k}-all-folds"
        restricted = run_dir / "restricted"
        public = run_dir / "public"
        metadata = {
            **{field: reference.get(field) for field in FROZEN_FIELDS if field != "input_artifact_sha256"},
            "schema_version": "blind-lot-pooled-run-1.0.0",
            "run_id": f"{manifest['pool_id']}-k{retrieval_k}",
            "retrieval_k": retrieval_k,
            "folds": manifest["expected_folds"],
            "input_artifact_sha256": reference["input_artifact_sha256"],
            "evaluation_version": EVALUATION_VERSION,
            "evaluation_only": True,
            "source_run_ids": [item["run_id"] for item in provenance[retrieval_k]],
        }
        rows = sorted(all_rows[retrieval_k], key=lambda row: row["case_key"])
        metrics = evaluate_joined(
            rows, bootstrap_seed=manifest["bootstrap_seed"],
            bootstrap_replicates=manifest["bootstrap_replicates"],
        )
        expected = manifest.get("frozen_regression", {}).get(str(retrieval_k))
        if expected is not None and _regression_values(metrics) != expected:
            raise RuntimeError(
                f"frozen k={retrieval_k} regression mismatch: "
                f"expected {expected}, got {_regression_values(metrics)}"
            )
        baseline = metrics["routing_policy_analysis"]["algorithm_cota_agreement"]
        expected_baseline = manifest.get("frozen_regression", {}).get("algorithm_cota")
        actual_baseline = {
            "accepted": baseline["accepted_count"], "correct": baseline["correct_accepted_count"],
            "incorrect": baseline["incorrect_accepted_count"],
        }
        if expected_baseline is not None and actual_baseline != expected_baseline:
            raise RuntimeError(f"frozen algorithm+COTA regression mismatch: {actual_baseline}")
        _write_json(restricted / "experiment_metadata.json", metadata)
        _write_json(restricted / "pooled_source_runs.json", {
            "schema_version": "1.0.0", "sources": provenance[retrieval_k],
        })
        _write_jsonl(restricted / "joined_evaluation.jsonl", rows)
        _write_jsonl(restricted / REAGGREGATED_JOINED_NAME, rows)
        _write_json(public / REAGGREGATED_AGGREGATE_NAME, _public_report(metadata, metrics))
        _write_json(restricted / "pooling_provenance.json", {
            "schema_version": "1.0.0", "evaluation_only": True,
            "provider_initialized": False, "retrieval_executed": False,
            "model_requests_made": 0, "predictions_modified": False,
            "manifest_sha256": sha256_file(manifest_path),
            "source_artifacts_verified_unchanged": True,
            "source_runs": provenance[retrieval_k],
        })
        outputs[retrieval_k] = run_dir

    changed = [str(path) for path, digest in before.items() if sha256_file(path) != digest]
    if changed:
        raise RuntimeError(f"source artifacts changed during pooling: {changed}")
    return outputs
