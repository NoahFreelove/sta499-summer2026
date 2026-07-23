#!/usr/bin/env python3
"""Export leakage-checked blind LOT prompt shards for coding-agent subagents."""

from __future__ import annotations

import argparse
from pathlib import Path

from blind_lot.agent_handoff import export_worker_bundle, prepare_agent_run
from blind_lot.runner import RunConfig
from lot_data import find_project_root
from run_blind_lot_benchmark import parse_folds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--fold", type=int, choices=range(5))
    modes.add_argument("--folds", type=parse_folds)
    modes.add_argument("--all-folds", action="store_true")
    parser.add_argument(
        "--orchestrator", choices=("codex", "claude-code"), default="codex"
    )
    parser.add_argument("--model-label", default="account-default")
    parser.add_argument("--reasoning-effort")
    parser.add_argument(
        "--run-purpose",
        choices=("smoke", "screening", "development", "final"),
        default="smoke",
    )
    parser.add_argument("--repeat-index", type=int, default=1)
    parser.add_argument("--retrieval-k", type=int, choices=(0, 3, 5), default=0)
    parser.add_argument(
        "--retrieval-training-folds",
        type=parse_folds,
        help="restrict few-shot retrieval candidates, for example 0,1,2",
    )
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--worker-bundle-dir",
        type=Path,
        help="sanitized directory outside the repository for worker-only chats",
    )
    parser.add_argument("--project-root", type=Path)
    parser.add_argument("--bootstrap-seed", type=int, default=4992026)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    args = parser.parse_args()
    if args.shards < 1:
        parser.error("--shards must be positive")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.repeat_index < 1:
        parser.error("--repeat-index must be positive")
    if args.run_purpose != "smoke" and args.model_label == "account-default":
        parser.error(
            "screening, development, and final runs require an exact --model-label"
        )
    if args.run_purpose != "smoke" and args.limit is not None:
        parser.error("screening, development, and final runs cannot use --limit")
    if args.run_purpose != "smoke" and args.worker_bundle_dir is None:
        parser.error(
            "screening, development, and final runs require --worker-bundle-dir"
        )
    return args


def main() -> None:
    args = parse_args()
    root = find_project_root(args.project_root)
    output = args.output_dir or Path("artifacts/benchmarks")
    if not output.is_absolute():
        output = (root / output).resolve()
    folds = (
        tuple(range(5))
        if args.all_folds
        else args.folds
        if args.folds is not None
        else (args.fold,)
    )
    expected_folds = {
        "screening": (0,),
        "development": (0, 1, 2),
        "final": (3, 4),
    }
    if args.run_purpose in expected_folds and folds != expected_folds[args.run_purpose]:
        raise SystemExit(
            f"{args.run_purpose} runs require --folds "
            + ",".join(str(value) for value in expected_folds[args.run_purpose])
        )
    if (
        args.run_purpose != "smoke"
        and args.retrieval_training_folds != (0, 1, 2)
    ):
        raise SystemExit(
            "screening, development, and final runs require "
            "--retrieval-training-folds 0,1,2"
        )
    provider = (
        "codex-subagents" if args.orchestrator == "codex" else "claude-code-subagents"
    )
    result = prepare_agent_run(
        RunConfig(
            root=root,
            output_dir=output,
            folds=folds,
            limit=args.limit,
            dry_run=False,
            provider=provider,
            model=args.model_label,
            reasoning_effort=args.reasoning_effort,
            temperature=None,
            retrieval_k=args.retrieval_k,
            concurrency=args.shards,
            request_timeout=1,
            retry_count=0,
            use_cache=False,
            bootstrap_seed=args.bootstrap_seed,
            bootstrap_replicates=args.bootstrap_replicates,
            retrieval_training_folds=args.retrieval_training_folds,
        ),
        shard_count=args.shards,
        run_purpose=args.run_purpose,
        repeat_index=args.repeat_index,
    )
    bundle = None
    if args.worker_bundle_dir is not None:
        bundle = export_worker_bundle(
            Path(result["run_dir"]),
            args.worker_bundle_dir,
            project_root=root,
            orchestrator=args.orchestrator,
        )
    print(f"Prepared {result['case_count']} cases in {result['shard_count']} shards.")
    print("Python API provider calls: 0")
    print(f"Manifest: {result['manifest']}")
    print(f"Worker outputs: {result['output_dir']}")
    if bundle is not None:
        print(f"Sanitized worker bundle: {bundle['bundle_dir']}")
        print(f"Start the worker-only chat in: {bundle['bundle_dir']}")
    print(f"Run directory: {result['run_dir']}")


if __name__ == "__main__":
    main()
