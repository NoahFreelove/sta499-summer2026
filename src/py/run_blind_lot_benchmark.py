#!/usr/bin/env python3
"""Run the leakage-safe order-only blind AI LOT benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path

from dotenv import load_dotenv

from blind_lot.runner import RunConfig, run
from lot_data import find_project_root


def format_metric(label: str, metric: dict[str, object]) -> str:
    """Format a count/denominator/rate metric, including undefined subset rates."""
    rate = metric["rate"]
    rendered_rate = "N/A" if rate is None else f"{rate:.1%}"
    return f"{label}: {rendered_rate} ({metric['count']}/{metric['denominator']})"


def load_project_env(root: Path) -> bool:
    """Load project-local credentials while preserving explicit shell values."""
    return load_dotenv(root / ".env", override=False)


def parse_folds(value: str) -> tuple[int, ...]:
    """Parse a comma-separated, unique subset of the five predefined folds."""
    try:
        folds = tuple(int(item.strip()) for item in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("folds must be comma-separated integers") from error
    if not folds or len(folds) != len(set(folds)) or set(folds) - set(range(5)):
        raise argparse.ArgumentTypeError("folds must be a unique, non-empty subset of 0..4")
    return folds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--fold", type=int, choices=range(5))
    modes.add_argument("--folds", type=parse_folds)
    modes.add_argument("--all-folds", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--provider", choices=("mock", "openai", "claude", "kimi"), default="mock")
    parser.add_argument("--model", default="mock-order-only-v1")
    parser.add_argument("--reasoning-effort")
    parser.add_argument("--temperature", type=float)
    parser.add_argument("--retrieval-k", type=int, choices=(0, 3, 5), default=0)
    parser.add_argument(
        "--retrieval-training-folds", type=parse_folds,
        help="restrict few-shot retrieval candidates, for example 0,1,2",
    )
    parser.add_argument("--concurrency", type=int, default=4)
    parser.add_argument("--request-timeout", type=float, default=120.0)
    parser.add_argument("--retry-count", type=int, default=2)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--use-cache", action="store_true")
    parser.add_argument("--project-root", type=Path)
    parser.add_argument("--bootstrap-seed", type=int, default=4992026)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.provider != "mock" and args.model == "mock-order-only-v1":
        parser.error(f"--provider {args.provider} requires an explicit --model")
    return args


def main() -> None:
    args = parse_args()
    root = find_project_root(args.project_root)
    load_project_env(root)
    output = args.output_dir or Path("artifacts/benchmarks")
    if not output.is_absolute():
        output = (root / output).resolve()
    folds = (
        tuple(range(5)) if args.all_folds
        else args.folds if args.folds is not None
        else (args.fold,)
    )
    result = run(RunConfig(
        root=root,
        output_dir=output,
        folds=folds,
        limit=args.limit,
        dry_run=args.dry_run,
        provider=args.provider,
        model=args.model,
        reasoning_effort=args.reasoning_effort,
        temperature=args.temperature,
        retrieval_k=args.retrieval_k,
        concurrency=args.concurrency,
        request_timeout=args.request_timeout,
        retry_count=args.retry_count,
        use_cache=args.use_cache,
        bootstrap_seed=args.bootstrap_seed,
        bootstrap_replicates=args.bootstrap_replicates,
        retrieval_training_folds=args.retrieval_training_folds,
    ))
    report = result["public"]
    if args.dry_run:
        print(f"Dry run validated {report['selected_target_count']} selected targets with no API calls.")
    else:
        metrics = report["metrics"]
        print(f"Evaluated {metrics['eligible_patients']} patients against reviewer consensus.")
        print(format_metric("Generated-total exact accuracy", metrics["generated_total"]["exact_accuracy"]))
        print(format_metric("Usable AI vote coverage", metrics["selective_ai"]["coverage"]))
        print(format_metric("Non-abstained exact accuracy", metrics["selective_ai"]["exact_accuracy"]))
        print(format_metric("Three-way consensus coverage", metrics["three_way_consensus"]["agreement"]))
    print(f"Run outputs: {result['run_dir']}")


if __name__ == "__main__":
    main()
