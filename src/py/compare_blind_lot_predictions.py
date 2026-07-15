#!/usr/bin/env python3
"""Run paired evaluation-only k=3 versus k=5 comparison."""

from __future__ import annotations

import argparse
from pathlib import Path

from blind_lot.paired_comparison import write_comparison


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k3-run", required=True, type=Path)
    parser.add_argument("--k5-run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=4992026)
    parser.add_argument("--replicates", type=int, default=2000)
    args = parser.parse_args()
    outputs = write_comparison(args.k3_run, args.k5_run, args.output, seed=args.seed, replicates=args.replicates)
    print(f"Paired public comparison: {outputs['json']}")


if __name__ == "__main__":
    main()
