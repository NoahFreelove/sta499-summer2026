"""Experiment orchestration for the order-only blind AI benchmark."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from lot_data import (
    assert_aggregate_only_report,
    case_key,
    read_jsonl,
    reconstruct_vendor_line_sequence,
)
from textbook_algo_cota import lot_algorithm_cota

from .cache import DiskCache, cache_key, canonical_hash
from .evaluation import (
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    EVALUATION_VERSION,
    evaluate_joined,
    normalize_joined_record,
)
from .models import BlindModelInput
from .prompting import PromptBundle, build_prompt, load_knowledge, load_prompt_template
from .providers import (
    ProviderConfig,
    complete_with_retries,
    create_provider,
)
from .retrieval import (
    LocalRetriever,
    RetrievalCandidate,
    RetrievalHit,
    assert_retrieval_safe,
)


@dataclass(frozen=True)
class RunConfig:
    root: Path
    output_dir: Path
    folds: tuple[int, ...]
    limit: int | None
    dry_run: bool
    provider: str
    model: str
    reasoning_effort: str | None
    temperature: float | None
    retrieval_k: int
    concurrency: int
    request_timeout: float
    retry_count: int
    use_cache: bool
    bootstrap_seed: int = BOOTSTRAP_SEED
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES
    retrieval_training_folds: tuple[int, ...] | None = None


@dataclass(frozen=True)
class Target:
    patient_key: str
    case_key: str
    fold: int
    exclusion_group: str
    blind: dict[str, Any]
    evaluation: dict[str, Any]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows: list[dict[str, Any]], key: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in sorted(rows, key=lambda item: item[key]):
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def _git_commit(root: Path) -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _load_manifest(path: Path) -> dict[str, dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["patient_key"]: {
                "fold": int(row["fold"]),
                "exclusion_group": row["exclusion_group"],
                "reviewer_consensus_lot": int(row["reviewer_consensus_lot"]),
            }
            for row in csv.DictReader(handle)
        }


def _load_inputs(config: RunConfig) -> tuple[list[Target], list[RetrievalCandidate], dict[str, str]]:
    paths = {
        "blind_adjudicated": config.root / "artifacts/blind/cota_adjudicated.jsonl",
        "restricted_evaluation": config.root / "artifacts/restricted/evaluation/cota_adjudicated.jsonl",
        "cv_manifest": config.root / "artifacts/restricted/cv/cota_adjudicated_5fold_manifest.csv",
        "exclusion_groups": config.root / "artifacts/restricted/overlap/exclusion_groups.jsonl",
        "blind_input_schema": config.root / "schemas/blind_model_input.schema.json",
    }
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)
    blind_rows = read_jsonl(paths["blind_adjudicated"])
    # Pydantic contract is equivalent to and intentionally stricter about event ordering.
    validated_blind = {
        row["case_key"]: BlindModelInput.model_validate(row).model_dump(mode="json")
        for row in blind_rows
    }
    evaluation_rows = read_jsonl(paths["restricted_evaluation"])
    manifest = _load_manifest(paths["cv_manifest"])
    evaluation_by_patient = {row["patient_key"]: row for row in evaluation_rows}

    candidates: list[RetrievalCandidate] = []
    targets: list[Target] = []
    for patient_key_value, manifest_row in sorted(manifest.items()):
        evaluation = evaluation_by_patient[patient_key_value]
        blind_case_key = case_key(patient_key_value)
        blind = validated_blind[blind_case_key]
        if evaluation["reviewer_lot"]["consensus"] != manifest_row["reviewer_consensus_lot"]:
            raise ValueError("manifest reviewer total disagrees with restricted evaluation record")
        candidate = RetrievalCandidate(
            case_key=blind_case_key,
            patient_key=patient_key_value,
            fold=manifest_row["fold"],
            exclusion_group=manifest_row["exclusion_group"],
            trajectory=blind["trajectory"],
            reviewer_consensus_lot=manifest_row["reviewer_consensus_lot"],
        )
        candidates.append(candidate)
        if manifest_row["fold"] in config.folds:
            targets.append(Target(
                patient_key=patient_key_value,
                case_key=blind_case_key,
                fold=manifest_row["fold"],
                exclusion_group=manifest_row["exclusion_group"],
                blind=blind,
                evaluation=evaluation,
            ))
    targets.sort(key=lambda target: (target.fold, target.case_key))
    if config.limit is not None:
        targets = targets[: config.limit]
    hashes = {name: _sha256_file(path) for name, path in paths.items()}
    return targets, candidates, hashes


def _run_id(config: RunConfig, prompt_version: str, knowledge_version: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    digest = canonical_hash({
        "folds": config.folds,
        "limit": config.limit,
        "dry_run": config.dry_run,
        "provider": config.provider,
        "model": config.model,
        "reasoning_effort": config.reasoning_effort,
        "temperature": config.temperature,
        "retrieval_k": config.retrieval_k,
        "retrieval_training_folds": config.retrieval_training_folds,
        "prompt_version": prompt_version,
        "knowledge_version": knowledge_version,
    })[:10]
    return f"{timestamp}-{digest}"


def _metadata(
    config: RunConfig,
    *,
    run_id: str,
    prompt_version: str,
    knowledge_version: str,
    input_hashes: dict[str, str],
) -> dict[str, Any]:
    return {
        "schema_version": "1.2.0",
        "run_id": run_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(config.root),
        "benchmark_mode": "order-only blind AI benchmark",
        "provider": config.provider,
        "model": config.model,
        "reasoning_effort": config.reasoning_effort,
        "temperature": config.temperature,
        "prompt_version": prompt_version,
        "knowledge_version": knowledge_version,
        "retrieval_k": config.retrieval_k,
        "retrieval_training_folds": (
            list(config.retrieval_training_folds)
            if config.retrieval_training_folds is not None
            else None
        ),
        "folds": list(config.folds),
        "limit": config.limit,
        "concurrency": config.concurrency,
        "request_timeout_seconds": config.request_timeout,
        "retry_count": config.retry_count,
        "use_cache": config.use_cache,
        "random_seed": config.bootstrap_seed,
        "bootstrap_replicates": config.bootstrap_replicates,
        "evaluation_version": EVALUATION_VERSION,
        "metric_semantics": {
            "abstained_outputs_are_votes": False,
            "abstained_numeric_totals_retained_for_diagnostics": True,
            "three_way_requires_non_abstained_ai": True,
            "primary_policy": "all-three agreement with usable AI",
            "pairwise_and_single_source_policies_are_exploratory": True,
        },
        "input_artifact_sha256": input_hashes,
    }


def _prepare(
    target: Target, retriever: LocalRetriever, retrieval_k: int
) -> tuple[PromptBundle, list[RetrievalHit]]:
    hits = retriever.retrieve(
        target.blind["trajectory"],
        target_patient_key=target.patient_key,
        target_fold=target.fold,
        target_exclusion_group=target.exclusion_group,
        k=retrieval_k,
    )
    assert_retrieval_safe(
        hits,
        target_patient_key=target.patient_key,
        target_fold=target.fold,
        target_exclusion_group=target.exclusion_group,
    )
    return build_prompt(target.blind, hits), hits


def _public_guard(report: dict[str, Any]) -> None:
    assert_aggregate_only_report(report)
    serialized = json.dumps(report, sort_keys=True)
    for marker in ("CASE-", "COTAOLD-", "COTANEW-", "EXCL-"):
        if marker in serialized:
            raise ValueError(f"public report contains patient-level identifier marker {marker}")


def run(config: RunConfig) -> dict[str, Any]:
    if config.retrieval_k not in {0, 3, 5}:
        raise ValueError("retrieval k must be 0, 3, or 5")
    if config.concurrency < 1 or config.retry_count < 0 or config.request_timeout <= 0:
        raise ValueError("invalid concurrency, retry, or timeout configuration")
    prompt_version, _ = load_prompt_template()
    knowledge_version, _ = load_knowledge()
    targets, candidates, input_hashes = _load_inputs(config)
    if not targets:
        raise ValueError("no eligible targets selected")
    if config.retrieval_training_folds is not None:
        invalid_folds = set(config.retrieval_training_folds) - set(range(5))
        if invalid_folds or not config.retrieval_training_folds:
            raise ValueError("retrieval training folds must be a non-empty subset of 0..4")
        candidates = [
            candidate for candidate in candidates
            if candidate.fold in config.retrieval_training_folds
        ]
    retriever = LocalRetriever(candidates)
    prepared = {target.case_key: _prepare(target, retriever, config.retrieval_k) for target in targets}

    run_id = _run_id(config, prompt_version, knowledge_version)
    run_dir = config.output_dir / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    restricted_dir, public_dir = run_dir / "restricted", run_dir / "public"
    metadata = _metadata(
        config, run_id=run_id, prompt_version=prompt_version,
        knowledge_version=knowledge_version, input_hashes=input_hashes,
    )
    _write_json(restricted_dir / "experiment_metadata.json", metadata)

    retrieval_debug = []
    for target in targets:
        bundle, hits = prepared[target.case_key]
        debug_hits = []
        for rank, (hit, rendered) in enumerate(
            zip(hits, bundle.rendered_retrieval_demonstrations), start=1
        ):
            debug_hits.append({
                "retrieval_rank": rank,
                "retrieved_case_key": hit.candidate.case_key,
                "overall_similarity_score": hit.score,
                "component_similarity_scores": hit.components,
                "blind_trajectory_supplied": hit.candidate.trajectory,
                "reviewer_consensus_lot_supplied": hit.candidate.reviewer_consensus_lot,
                "transition_decisions_or_explanation_supplied": None,
                "exact_rendered_demonstration_text": rendered,
                "rendered_demonstration_sha256": _sha256_text(rendered),
            })
        retrieval_debug.append({
            "case_key": target.case_key,
            "fold": target.fold,
            "retrieval_k": config.retrieval_k,
            "retrieval_training_folds": (
                list(config.retrieval_training_folds)
                if config.retrieval_training_folds is not None
                else None
            ),
            "prompt_version": bundle.prompt_version,
            "knowledge_version": bundle.knowledge_version,
            "hits": debug_hits,
            "exact_rendered_retrieval_context": bundle.rendered_retrieval_context,
            "retrieval_context_sha256": _sha256_text(bundle.rendered_retrieval_context),
            "rendered_model_input_sha256": _sha256_text(bundle.patient_prompt),
        })
    _write_jsonl(restricted_dir / "retrieval_debug.jsonl", retrieval_debug, "case_key")

    if config.dry_run:
        previews = [
            {
                "case_key": target.case_key,
                "fold": target.fold,
                "stable_prefix": prepared[target.case_key][0].stable_prefix,
                "patient_prompt": json.loads(prepared[target.case_key][0].patient_prompt),
            }
            for target in targets[:3]
        ]
        _write_json(restricted_dir / "prompt_previews.json", previews)
        public = {
            "schema_version": "1.0.0",
            "report_scope": "aggregate_only",
            "benchmark_mode": "order-only blind AI benchmark",
            "run_id_hash": hashlib.sha256(run_id.encode()).hexdigest(),
            "dry_run": True,
            "selected_target_count": len(targets),
            "validated_blind_input_count": len(candidates),
            "retrieval_k": config.retrieval_k,
            "retrieval_training_folds": (
                list(config.retrieval_training_folds)
                if config.retrieval_training_folds is not None
                else None
            ),
            "leakage_checks": {
                "blind_schema_valid": True,
                "target_prompt_allowlist_valid": True,
                "held_out_fold_exclusion_valid": True,
                "exclusion_group_exclusion_valid": True,
                "raw_identifier_loaded_or_serialized": False,
            },
        }
        _public_guard(public)
        _write_json(public_dir / "dry_run_report.json", public)
        return {"run_dir": str(run_dir), "public": public}

    provider_config = ProviderConfig(
        model=config.model,
        reasoning_effort=config.reasoning_effort,
        temperature=config.temperature,
        request_timeout=config.request_timeout,
        retry_count=config.retry_count,
    )
    provider = create_provider(config.provider, provider_config)
    cache = DiskCache(config.output_dir / "cache")

    def infer(target: Target) -> dict[str, Any]:
        bundle, hits = prepared[target.case_key]
        key = cache_key(
            provider=config.provider,
            model=config.model,
            reasoning_effort=config.reasoning_effort,
            temperature=config.temperature,
            prompt_version=bundle.prompt_version,
            knowledge_version=bundle.knowledge_version,
            retrieval_k=config.retrieval_k,
            retrieved_example_ids=bundle.retrieved_example_ids,
            blind_patient_input_hash=canonical_hash(target.blind),
        )
        response = cache.get(key, len(target.blind["trajectory"])) if config.use_cache else None
        cache_hit = response is not None
        attempts = 0
        if response is not None:
            cited = {
                case_id_value
                for transition in response.transitions
                for case_id_value in transition.retrieved_case_ids
            }
            if not cited <= set(bundle.retrieved_example_ids):
                raise ValueError("cached response cites a retrieved case that was not supplied")
        if response is None:
            response, attempts = complete_with_retries(
                provider, bundle.stable_prefix, bundle.patient_prompt,
                regimen_event_count=len(target.blind["trajectory"]),
                retry_count=config.retry_count,
                permitted_retrieved_case_ids=set(bundle.retrieved_example_ids),
            )
            cache.put(key, response)
        return {
            "case_key": target.case_key,
            "fold": target.fold,
            "response": response.model_dump(mode="json"),
            "retrieved_case_ids": [hit.candidate.case_key for hit in hits],
            "cache_key": key,
            "cache_hit": cache_hit,
            "provider_attempts": attempts,
        }

    predictions: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=config.concurrency) as executor:
        futures = {executor.submit(infer, target): target.case_key for target in targets}
        for future in as_completed(futures):
            predictions.append(future.result())
    # This file is finalized before loading any target answer into the evaluation join below.
    _write_jsonl(restricted_dir / "blind_predictions.jsonl", predictions, "case_key")

    prediction_by_case = {row["case_key"]: row for row in predictions}
    joined: list[dict[str, Any]] = []
    for target in targets:
        prediction = prediction_by_case[target.case_key]["response"]
        evaluation = target.evaluation
        reviewer = evaluation["reviewer_lot"]["consensus"]
        cota = evaluation["cota_lot"]
        if reviewer is None or cota is None:
            continue
        algorithm_lot, algorithm_flags = lot_algorithm_cota(reconstruct_vendor_line_sequence(evaluation))
        joined.append(normalize_joined_record({
            "case_key": target.case_key,
            "fold": target.fold,
            "reviewer_lot": reviewer,
            "ai_lot": prediction["ai_lot_count"],
            "ai_abstained": prediction["abstained"],
            "algorithm_lot": algorithm_lot,
            "cota_lot": cota,
            "algorithm_flags": sorted(algorithm_flags),
        }))
    _write_jsonl(restricted_dir / "joined_evaluation.jsonl", joined, "case_key")
    metrics = evaluate_joined(
        joined,
        bootstrap_seed=config.bootstrap_seed,
        bootstrap_replicates=config.bootstrap_replicates,
    )
    public = {
        "schema_version": "1.2.0",
        "report_scope": "aggregate_only",
        "benchmark_mode": "order-only blind AI benchmark",
        "ground_truth": "reviewer consensus",
        "cota_role": "independent vote, not a training label",
        "run_id_hash": hashlib.sha256(run_id.encode()).hexdigest(),
        "experiment": {
            "provider": config.provider,
            "model": config.model,
            "reasoning_effort": config.reasoning_effort,
            "prompt_version": prompt_version,
            "knowledge_version": knowledge_version,
            "retrieval_k": config.retrieval_k,
            "retrieval_training_folds": (
                list(config.retrieval_training_folds)
                if config.retrieval_training_folds is not None
                else None
            ),
            "folds": list(config.folds),
            "fold_count_evaluated": len(config.folds),
            "limited_run": config.limit is not None,
            "temperature": config.temperature,
            "bootstrap_seed": config.bootstrap_seed,
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
    _public_guard(public)
    _write_json(public_dir / "aggregate_evaluation.json", public)
    return {"run_dir": str(run_dir), "public": public}
