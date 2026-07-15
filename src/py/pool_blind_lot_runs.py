#!/usr/bin/env python3
"""Pool frozen blind LOT fold runs without inference or retrieval."""

from __future__ import annotations

import argparse
from pathlib import Path

from blind_lot.pooling import pool_runs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--runs-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    outputs = pool_runs(args.manifest, args.runs_root, args.output)
    for retrieval_k, path in outputs.items():
        print(f"k={retrieval_k} pooled evaluation: {path}")


if __name__ == "__main__":
    main()
