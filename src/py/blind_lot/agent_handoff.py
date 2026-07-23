"""File handoff between the blind LOT benchmark and interactive coding agents."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

from lot_data import reconstruct_vendor_line_sequence
from textbook_algo_cota import lot_algorithm_cota

from .cache import canonical_hash
from .evaluation import evaluate_joined, normalize_joined_record
from .models import response_json_schema, validate_model_response
from .prompting import load_knowledge, load_prompt_template
from .retrieval import LocalRetriever
from .runner import (
    RunConfig,
    Target,
    _load_inputs,
    _metadata,
    _prepare,
    _public_guard,
    _run_id,
    _write_json,
    _write_jsonl,
)


HANDOFF_SCHEMA_VERSION = "blind-lot-agent-handoff-1.0.0"
SHARD_OUTPUT_SCHEMA_VERSION = "blind-lot-agent-shard-output-1.0.0"
WORKER_BUNDLE_SCHEMA_VERSION = "blind-lot-worker-bundle-1.0.0"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def _partition(
    items: list[dict[str, Any]], shard_count: int
) -> list[list[dict[str, Any]]]:
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    actual_count = min(shard_count, len(items))
    return [items[index::actual_count] for index in range(actual_count)]


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _bundle_skill(orchestrator: str) -> str:
    worker_name = "blind_lot_worker" if orchestrator == "codex" else "blind-lot-worker"
    return f"""---
name: blind-lot-bundle
description: Process a sanitized blind LOT worker bundle with parallel subagents.
---

# Process this blind LOT bundle

Read only `manifest.json` and the packet assigned to each worker. Spawn one
`{worker_name}` per shard, up to the concurrency limit. Workers must write only
their declared `outputs/shard-*.json` file and must include `worker_provenance`
from the packet contract.

Do not search outside this bundle, compare cases, evaluate predictions, or read
another worker's packet. Wait for every shard. Preserve a failed output under
`failed_attempts/` before retrying, increment `worker_provenance.attempt`, and
record the agent/thread ID when available. Return only completion status, shard
attempt counts, and the bundle path.
"""


def _codex_worker(model: str, reasoning_effort: str | None) -> str:
    effort = (
        f"model_reasoning_effort = {json.dumps(reasoning_effort)}\n"
        if reasoning_effort is not None
        else ""
    )
    return f"""name = "blind_lot_worker"
description = "Classify one assigned blind LOT packet without consulting other data."
model = {json.dumps(model)}
{effort}sandbox_mode = "workspace-write"
developer_instructions = \"\"\"
Read only the assigned packet. Do not search outside the worker bundle or read
another packet. Use only stable_prefix and patient_prompt. Write only the
declared output JSON, including truthful worker_provenance. Do not evaluate,
compare, summarize, or alter source inputs.
\"\"\"
"""


def _claude_worker(model: str, reasoning_effort: str | None) -> str:
    effort = (
        f"effort: {reasoning_effort}\n" if reasoning_effort is not None else ""
    )
    return f"""---
name: blind-lot-worker
description: Classify one assigned blind LOT packet without consulting other data.
tools: Read, Write
model: {model}
{effort}---

Read only the assigned packet. Do not search outside the worker bundle or read
another packet. Use only `stable_prefix` and `patient_prompt`. Write only the
declared output JSON, including truthful `worker_provenance`. Do not evaluate,
compare, summarize, or alter source inputs.
"""


def _retrieval_debug_record(
    target: Target, bundle: Any, hits: list[Any], retrieval_k: int
) -> dict[str, Any]:
    debug_hits = []
    for rank, (hit, rendered) in enumerate(
        zip(hits, bundle.rendered_retrieval_demonstrations), start=1
    ):
        debug_hits.append(
            {
                "retrieval_rank": rank,
                "retrieved_case_key": hit.candidate.case_key,
                "overall_similarity_score": hit.score,
                "component_similarity_scores": hit.components,
                "blind_trajectory_supplied": hit.candidate.trajectory,
                "reviewer_consensus_lot_supplied": hit.candidate.reviewer_consensus_lot,
                "transition_decisions_or_explanation_supplied": None,
                "exact_rendered_demonstration_text": rendered,
                "rendered_demonstration_sha256": hashlib.sha256(
                    rendered.encode()
                ).hexdigest(),
            }
        )
    return {
        "case_key": target.case_key,
        "fold": target.fold,
        "retrieval_k": retrieval_k,
        "prompt_version": bundle.prompt_version,
        "knowledge_version": bundle.knowledge_version,
        "hits": debug_hits,
        "exact_rendered_retrieval_context": bundle.rendered_retrieval_context,
        "retrieval_context_sha256": hashlib.sha256(
            bundle.rendered_retrieval_context.encode()
        ).hexdigest(),
        "rendered_model_input_sha256": hashlib.sha256(
            bundle.patient_prompt.encode()
        ).hexdigest(),
    }


def prepare_agent_run(
    config: RunConfig,
    *,
    shard_count: int,
    run_purpose: str = "smoke",
    repeat_index: int = 1,
) -> dict[str, Any]:
    """Prepare blind prompt shards without initializing an API provider."""
    if config.provider not in {"codex-subagents", "claude-code-subagents"}:
        raise ValueError(
            "agent handoff provider must be codex-subagents or claude-code-subagents"
        )
    if config.retrieval_k not in {0, 3, 5}:
        raise ValueError("retrieval k must be 0, 3, or 5")
    if run_purpose not in {"smoke", "screening", "development", "final"}:
        raise ValueError("run purpose must be smoke, screening, development, or final")
    if repeat_index < 1:
        raise ValueError("repeat index must be positive")
    if run_purpose != "smoke" and config.model == "account-default":
        raise ValueError("screening, development, and final runs require an exact model label")
    if run_purpose != "smoke" and config.limit is not None:
        raise ValueError(
            "screening, development, and final benchmark runs cannot use a patient limit"
        )
    expected_target_folds = {
        "screening": (0,),
        "development": (0, 1, 2),
        "final": (3, 4),
    }
    if run_purpose in expected_target_folds and config.folds != expected_target_folds[
        run_purpose
    ]:
        raise ValueError(
            f"{run_purpose} runs require target folds "
            f"{expected_target_folds[run_purpose]}"
        )
    if run_purpose != "smoke" and config.retrieval_training_folds != (0, 1, 2):
        raise ValueError(
            "screening, development, and final runs require retrieval training folds "
            "(0, 1, 2)"
        )

    prompt_version, _ = load_prompt_template()
    knowledge_version, _ = load_knowledge()
    targets, candidates, input_hashes = _load_inputs(config)
    if not targets:
        raise ValueError("no eligible targets selected")
    if config.retrieval_training_folds is not None:
        invalid_folds = set(config.retrieval_training_folds) - set(range(5))
        if invalid_folds or not config.retrieval_training_folds:
            raise ValueError(
                "retrieval training folds must be a non-empty subset of 0..4"
            )
        candidates = [
            candidate
            for candidate in candidates
            if candidate.fold in config.retrieval_training_folds
        ]
    retriever = LocalRetriever(candidates)
    prepared = {
        target.case_key: _prepare(target, retriever, config.retrieval_k)
        for target in targets
    }

    run_id = _run_id(config, prompt_version, knowledge_version)
    run_dir = config.output_dir / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    restricted_dir = run_dir / "restricted"
    handoff_dir = restricted_dir / "agent_handoff"
    packets_dir = handoff_dir / "packets"
    outputs_dir = handoff_dir / "outputs"
    packets_dir.mkdir(parents=True)
    outputs_dir.mkdir()

    metadata = _metadata(
        config,
        run_id=run_id,
        prompt_version=prompt_version,
        knowledge_version=knowledge_version,
        input_hashes=input_hashes,
    )
    metadata.update(
        {
            "inference_mode": "interactive_subagent_handoff",
            "agent_handoff_schema_version": HANDOFF_SCHEMA_VERSION,
            "api_provider_initialized": False,
            "model_api_requests_made_by_python": 0,
            "requested_shard_count": shard_count,
            "run_purpose": run_purpose,
            "repeat_index": repeat_index,
            "exact_model_label_required": run_purpose != "smoke",
        }
    )
    _write_json(restricted_dir / "experiment_metadata.json", metadata)

    retrieval_debug = []
    cases = []
    schema = response_json_schema()
    shared_stable_prefix: str | None = None
    for target in targets:
        bundle, hits = prepared[target.case_key]
        if shared_stable_prefix is None:
            shared_stable_prefix = bundle.stable_prefix
        elif bundle.stable_prefix != shared_stable_prefix:
            raise ValueError("prepared cases do not share one stable prompt prefix")
        retrieval_debug.append(
            _retrieval_debug_record(target, bundle, hits, config.retrieval_k)
        )
        cases.append(
            {
                "case_key": target.case_key,
                "fold": target.fold,
                "regimen_event_count": len(target.blind["trajectory"]),
                "permitted_retrieved_case_ids": bundle.retrieved_example_ids,
                "patient_prompt": bundle.patient_prompt,
                "prompt_sha256": canonical_hash(
                    {
                        "stable_prefix": bundle.stable_prefix,
                        "patient_prompt": bundle.patient_prompt,
                    }
                ),
            }
        )
    _write_jsonl(restricted_dir / "retrieval_debug.jsonl", retrieval_debug, "case_key")

    shard_entries = []
    for index, shard_cases in enumerate(_partition(cases, shard_count), start=1):
        shard_id = f"shard-{index:03d}"
        packet_path = packets_dir / f"{shard_id}.json"
        packet = {
            "schema_version": HANDOFF_SCHEMA_VERSION,
            "run_id": run_id,
            "shard_id": shard_id,
            "stable_prefix": shared_stable_prefix,
            "response_json_schema": schema,
            "instructions": (
                "Classify every case independently. Use only this packet's stable_prefix "
                "and each case's patient_prompt. Return one response per case and no prose."
            ),
            "output_contract": {
                "schema_version": SHARD_OUTPUT_SCHEMA_VERSION,
                "run_id": run_id,
                "shard_id": shard_id,
                "predictions": [
                    {
                        "case_key": "<case_key>",
                        "response": "<response_json_schema object>",
                    }
                ],
                "worker_provenance": {
                    "configured_model": config.model,
                    "configured_reasoning_effort": config.reasoning_effort,
                    "attempt": 1,
                    "agent_id": "<agent/thread id when available, otherwise null>",
                },
            },
            "cases": shard_cases,
        }
        _write_json(packet_path, packet)
        shard_entries.append(
            {
                "shard_id": shard_id,
                "packet": str(packet_path.relative_to(run_dir)),
                "output": str((outputs_dir / f"{shard_id}.json").relative_to(run_dir)),
                "case_count": len(shard_cases),
                "case_keys": [case["case_key"] for case in shard_cases],
                "packet_sha256": _sha256_file(packet_path),
            }
        )

    manifest = {
        "schema_version": HANDOFF_SCHEMA_VERSION,
        "run_id": run_id,
        "provider": config.provider,
        "model": config.model,
        "reasoning_effort": config.reasoning_effort,
        "run_purpose": run_purpose,
        "repeat_index": repeat_index,
        "retrieval_training_folds": (
            list(config.retrieval_training_folds)
            if config.retrieval_training_folds is not None
            else None
        ),
        "case_count": len(cases),
        "shard_count": len(shard_entries),
        "shards": shard_entries,
    }
    manifest_path = handoff_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    return {
        "run_dir": str(run_dir),
        "manifest": str(manifest_path),
        "output_dir": str(outputs_dir),
        "case_count": len(cases),
        "shard_count": len(shard_entries),
    }


def export_worker_bundle(
    run_dir: Path,
    bundle_dir: Path,
    *,
    project_root: Path,
    orchestrator: str,
) -> dict[str, Any]:
    """Create a sanitized worker directory outside the benchmark repository."""
    run_dir = run_dir.resolve()
    project_root = project_root.resolve()
    bundle_dir = bundle_dir.expanduser().resolve()
    if orchestrator not in {"codex", "claude-code"}:
        raise ValueError("orchestrator must be codex or claude-code")
    if _is_within(bundle_dir, project_root):
        raise ValueError("worker bundle must be outside the project root")
    bundle_dir.mkdir(parents=True, exist_ok=False)
    packets_dir = bundle_dir / "packets"
    outputs_dir = bundle_dir / "outputs"
    failed_dir = bundle_dir / "failed_attempts"
    packets_dir.mkdir()
    outputs_dir.mkdir()
    failed_dir.mkdir()

    manifest = _read_json(run_dir / "restricted" / "agent_handoff" / "manifest.json")
    metadata = _read_json(run_dir / "restricted" / "experiment_metadata.json")
    if manifest["run_id"] != metadata["run_id"]:
        raise ValueError("run metadata and handoff manifest disagree")

    worker_shards = []
    for shard in manifest["shards"]:
        source = run_dir / shard["packet"]
        destination = packets_dir / f"{shard['shard_id']}.json"
        if _sha256_file(source) != shard["packet_sha256"]:
            raise ValueError(f"source packet hash mismatch: {source}")
        shutil.copyfile(source, destination)
        if _sha256_file(destination) != shard["packet_sha256"]:
            raise RuntimeError(f"worker packet copy hash mismatch: {destination}")
        worker_shards.append(
            {
                "shard_id": shard["shard_id"],
                "packet": str(destination.relative_to(bundle_dir)),
                "output": str(
                    (outputs_dir / f"{shard['shard_id']}.json").relative_to(bundle_dir)
                ),
                "failed_attempts_dir": str(failed_dir.relative_to(bundle_dir)),
                "case_count": shard["case_count"],
                "case_keys": shard["case_keys"],
                "packet_sha256": shard["packet_sha256"],
            }
        )

    worker_manifest = {
        "schema_version": WORKER_BUNDLE_SCHEMA_VERSION,
        "handoff_schema_version": HANDOFF_SCHEMA_VERSION,
        "run_id": manifest["run_id"],
        "provider": manifest["provider"],
        "orchestrator": orchestrator,
        "configured_model": manifest["model"],
        "configured_reasoning_effort": manifest["reasoning_effort"],
        "run_purpose": manifest["run_purpose"],
        "repeat_index": manifest["repeat_index"],
        "case_count": manifest["case_count"],
        "shard_count": manifest["shard_count"],
        "shards": worker_shards,
    }
    _write_json(bundle_dir / "manifest.json", worker_manifest)

    if orchestrator == "codex":
        skill_dir = bundle_dir / ".agents" / "skills" / "blind-lot-bundle"
        agent_dir = bundle_dir / ".codex" / "agents"
        skill_dir.mkdir(parents=True)
        agent_dir.mkdir(parents=True)
        (skill_dir / "SKILL.md").write_text(
            _bundle_skill(orchestrator), encoding="utf-8"
        )
        (agent_dir / "blind-lot-worker.toml").write_text(
            _codex_worker(manifest["model"], manifest["reasoning_effort"]),
            encoding="utf-8",
        )
    else:
        skill_dir = bundle_dir / ".claude" / "skills" / "blind-lot-bundle"
        agent_dir = bundle_dir / ".claude" / "agents"
        skill_dir.mkdir(parents=True)
        agent_dir.mkdir(parents=True)
        (skill_dir / "SKILL.md").write_text(
            _bundle_skill(orchestrator), encoding="utf-8"
        )
        (agent_dir / "blind-lot-worker.md").write_text(
            _claude_worker(manifest["model"], manifest["reasoning_effort"]),
            encoding="utf-8",
        )

    _write_json(
        run_dir / "restricted" / "agent_handoff" / "worker_bundle.json",
        {
            "schema_version": WORKER_BUNDLE_SCHEMA_VERSION,
            "bundle_dir": str(bundle_dir),
            "bundle_manifest_sha256": _sha256_file(bundle_dir / "manifest.json"),
            "orchestrator": orchestrator,
            "contains_evaluation_data": False,
            "packet_hashes_verified": True,
        },
    )
    return {
        "bundle_dir": str(bundle_dir),
        "manifest": str(bundle_dir / "manifest.json"),
        "shard_count": len(worker_shards),
    }


def import_worker_bundle_outputs(run_dir: Path, bundle_dir: Path) -> dict[str, Any]:
    """Verify and archive worker outputs from a sanitized bundle."""
    run_dir = run_dir.resolve()
    bundle_dir = bundle_dir.expanduser().resolve()
    manifest = _read_json(run_dir / "restricted" / "agent_handoff" / "manifest.json")
    worker_manifest_path = bundle_dir / "manifest.json"
    worker_manifest = _read_json(worker_manifest_path)
    if worker_manifest.get("schema_version") != WORKER_BUNDLE_SCHEMA_VERSION:
        raise ValueError("unsupported worker bundle")
    if worker_manifest.get("run_id") != manifest.get("run_id"):
        raise ValueError("worker bundle does not belong to this run")

    internal_by_id = {shard["shard_id"]: shard for shard in manifest["shards"]}
    output_hashes = {}
    for worker_shard in worker_manifest["shards"]:
        shard_id = worker_shard["shard_id"]
        internal = internal_by_id.get(shard_id)
        if internal is None:
            raise ValueError(f"unexpected worker shard: {shard_id}")
        packet = bundle_dir / worker_shard["packet"]
        if _sha256_file(packet) != internal["packet_sha256"]:
            raise ValueError(f"worker packet changed: {packet}")
        source = bundle_dir / worker_shard["output"]
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = run_dir / internal["output"]
        if destination.exists():
            raise FileExistsError(
                f"refusing to overwrite imported output: {destination}"
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        source_hash = _sha256_file(source)
        if _sha256_file(destination) != source_hash:
            raise RuntimeError(f"worker output copy hash mismatch: {destination}")
        output_hashes[shard_id] = source_hash

    if set(output_hashes) != set(internal_by_id):
        raise ValueError("worker bundle shard coverage is incomplete")
    provenance = {
        "schema_version": WORKER_BUNDLE_SCHEMA_VERSION,
        "bundle_dir": str(bundle_dir),
        "bundle_manifest_sha256": _sha256_file(worker_manifest_path),
        "output_sha256": output_hashes,
        "packet_hashes_verified": True,
        "outputs_archived_before_finalization": True,
    }
    provenance_path = run_dir / "restricted" / "agent_handoff" / "bundle_import.json"
    _write_json(provenance_path, provenance)
    return {"provenance": str(provenance_path), "output_sha256": output_hashes}


def validate_prediction_shards(
    run_dir: Path, manifest: dict[str, Any]
) -> list[dict[str, Any]]:
    """Validate exact shard coverage and model response contracts."""
    if manifest.get("schema_version") != HANDOFF_SCHEMA_VERSION:
        raise ValueError("unsupported agent handoff manifest")
    predictions: list[dict[str, Any]] = []
    seen: set[str] = set()

    for shard in manifest["shards"]:
        packet_path = run_dir / shard["packet"]
        output_path = run_dir / shard["output"]
        if _sha256_file(packet_path) != shard["packet_sha256"]:
            raise ValueError(f"prompt packet changed after preparation: {packet_path}")
        if not output_path.exists():
            raise FileNotFoundError(output_path)
        packet = _read_json(packet_path)
        output = _read_json(output_path)
        for field in ("run_id", "shard_id"):
            if output.get(field) != packet[field]:
                raise ValueError(f"{field} mismatch in {output_path}")
        if output.get("schema_version") != SHARD_OUTPUT_SCHEMA_VERSION:
            raise ValueError(f"unsupported shard output schema in {output_path}")
        provenance = output.get("worker_provenance")
        if manifest.get("run_purpose") != "smoke":
            if not isinstance(provenance, dict):
                raise ValueError(f"worker provenance is required in {output_path}")
            if provenance.get("configured_model") != manifest.get("model"):
                raise ValueError(f"configured model mismatch in {output_path}")
            if provenance.get("configured_reasoning_effort") != manifest.get(
                "reasoning_effort"
            ):
                raise ValueError(
                    f"configured reasoning effort mismatch in {output_path}"
                )
            attempt = provenance.get("attempt")
            if isinstance(attempt, bool) or not isinstance(attempt, int) or attempt < 1:
                raise ValueError(f"invalid worker attempt in {output_path}")

        expected = {case["case_key"]: case for case in packet["cases"]}
        rows = output.get("predictions")
        if not isinstance(rows, list):
            raise ValueError(f"predictions must be an array in {output_path}")
        received_keys = [row.get("case_key") for row in rows if isinstance(row, dict)]
        if len(received_keys) != len(rows) or set(received_keys) != set(expected):
            raise ValueError(f"shard coverage mismatch in {output_path}")
        if len(received_keys) != len(set(received_keys)):
            raise ValueError(f"duplicate case key in {output_path}")

        for row in rows:
            case_key = row["case_key"]
            if case_key in seen:
                raise ValueError(f"duplicate prediction across shards: {case_key}")
            seen.add(case_key)
            case = expected[case_key]
            response = validate_model_response(
                row["response"], case["regimen_event_count"]
            )
            cited = {
                cited_case
                for transition in response.transitions
                for cited_case in transition.retrieved_case_ids
            }
            permitted = set(case["permitted_retrieved_case_ids"])
            if not cited <= permitted:
                raise ValueError(
                    f"{case_key} cites a retrieved case that was not supplied"
                )
            predictions.append(
                {
                    "case_key": case_key,
                    "fold": case["fold"],
                    "response": response.model_dump(mode="json"),
                    "retrieved_case_ids": case["permitted_retrieved_case_ids"],
                    "cache_key": case["prompt_sha256"],
                    "cache_hit": False,
                    "provider_attempts": 1,
                    "agent_shard_id": shard["shard_id"],
                    "agent_output_sha256": _sha256_file(output_path),
                    "worker_provenance": provenance,
                }
            )

    expected_all = {
        case_key for shard in manifest["shards"] for case_key in shard["case_keys"]
    }
    if seen != expected_all or len(seen) != manifest["case_count"]:
        raise ValueError("combined shard coverage does not match the manifest")
    return sorted(predictions, key=lambda row: row["case_key"])


def finalize_agent_run(run_dir: Path, root: Path) -> dict[str, Any]:
    """Validate subagent outputs, freeze predictions, and run the existing evaluator."""
    run_dir = run_dir.resolve()
    restricted_dir = run_dir / "restricted"
    public_dir = run_dir / "public"
    metadata = _read_json(restricted_dir / "experiment_metadata.json")
    manifest = _read_json(restricted_dir / "agent_handoff" / "manifest.json")
    if manifest.get("run_id") != metadata.get("run_id"):
        raise ValueError("handoff manifest does not belong to this run")

    prediction_path = restricted_dir / "blind_predictions.jsonl"
    joined_path = restricted_dir / "joined_evaluation.jsonl"
    aggregate_path = public_dir / "aggregate_evaluation.json"
    for path in (prediction_path, joined_path, aggregate_path):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite finalized artifact: {path}")

    predictions = validate_prediction_shards(run_dir, manifest)
    _write_jsonl(prediction_path, predictions, "case_key")
    frozen_prediction_sha256 = _sha256_file(prediction_path)

    config = RunConfig(
        root=root.resolve(),
        output_dir=run_dir.parent.parent,
        folds=tuple(int(value) for value in metadata["folds"]),
        limit=metadata.get("limit"),
        dry_run=False,
        provider=metadata["provider"],
        model=metadata["model"],
        reasoning_effort=metadata.get("reasoning_effort"),
        temperature=metadata.get("temperature"),
        retrieval_k=int(metadata["retrieval_k"]),
        concurrency=1,
        request_timeout=1,
        retry_count=0,
        use_cache=False,
        bootstrap_seed=int(metadata["random_seed"]),
        bootstrap_replicates=int(metadata["bootstrap_replicates"]),
        retrieval_training_folds=(
            tuple(int(value) for value in metadata["retrieval_training_folds"])
            if metadata.get("retrieval_training_folds") is not None
            else None
        ),
    )
    targets, _, current_hashes = _load_inputs(config)
    if current_hashes != metadata["input_artifact_sha256"]:
        raise ValueError("benchmark inputs changed after prompt preparation")
    target_by_case = {target.case_key: target for target in targets}
    prediction_by_case = {row["case_key"]: row for row in predictions}
    if set(target_by_case) != set(prediction_by_case):
        raise ValueError("final target cohort does not match prediction coverage")

    joined = []
    for case_key in sorted(target_by_case):
        target = target_by_case[case_key]
        prediction = prediction_by_case[case_key]["response"]
        evaluation = target.evaluation
        reviewer = evaluation["reviewer_lot"]["consensus"]
        cota = evaluation["cota_lot"]
        if reviewer is None or cota is None:
            continue
        algorithm_lot, algorithm_flags = lot_algorithm_cota(
            reconstruct_vendor_line_sequence(evaluation)
        )
        joined.append(
            normalize_joined_record(
                {
                    "case_key": case_key,
                    "fold": target.fold,
                    "reviewer_lot": reviewer,
                    "ai_lot": prediction["ai_lot_count"],
                    "ai_abstained": prediction["abstained"],
                    "algorithm_lot": algorithm_lot,
                    "cota_lot": cota,
                    "algorithm_flags": sorted(algorithm_flags),
                }
            )
        )
    _write_jsonl(joined_path, joined, "case_key")
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
        "run_id_hash": hashlib.sha256(metadata["run_id"].encode()).hexdigest(),
        "experiment": {
            "provider": config.provider,
            "model": config.model,
            "reasoning_effort": config.reasoning_effort,
            "prompt_version": metadata["prompt_version"],
            "knowledge_version": metadata["knowledge_version"],
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
            "evaluation_version": metadata["evaluation_version"],
            "inference_mode": metadata["inference_mode"],
            "run_purpose": metadata["run_purpose"],
            "repeat_index": metadata["repeat_index"],
            "metric_semantics": metadata["metric_semantics"],
        },
        "metrics": metrics,
    }
    _public_guard(public)
    _write_json(aggregate_path, public)
    if _sha256_file(prediction_path) != frozen_prediction_sha256:
        raise RuntimeError("blind predictions changed during evaluation")
    _write_json(
        restricted_dir / "agent_handoff" / "finalization.json",
        {
            "schema_version": HANDOFF_SCHEMA_VERSION,
            "predictions_finalized_before_evaluation_join": True,
            "blind_predictions_sha256": frozen_prediction_sha256,
            "model_api_requests_made_by_python": 0,
            "validated_prediction_count": len(predictions),
            "joined_evaluation_count": len(joined),
            "run_purpose": metadata["run_purpose"],
            "repeat_index": metadata["repeat_index"],
            "worker_attempts": {
                shard["shard_id"]: _read_json(run_dir / shard["output"])
                .get("worker_provenance", {})
                .get("attempt")
                for shard in manifest["shards"]
            },
            "worker_agent_ids": {
                shard["shard_id"]: _read_json(run_dir / shard["output"])
                .get("worker_provenance", {})
                .get("agent_id")
                for shard in manifest["shards"]
            },
        },
    )
    return {
        "run_dir": str(run_dir),
        "prediction_output": str(prediction_path),
        "aggregate_output": str(aggregate_path),
        "public": public,
    }
