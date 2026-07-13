"""Build de-identified LOT datasets, quality reports, and CV folds."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from lot_data import (
    duplicate_summary,
    load_cota,
    load_preamble,
    make_folds,
    malformed_summary,
    overlap_summary,
    write_fold_manifest,
    write_json,
    write_jsonl,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OLD = ROOT / "data" / "LoT Adjudication Datasets.xlsx"
DEFAULT_NEW = ROOT / "data" / "Preamble_treatment_cleaned_081425_adjudicated by JK and AR 20250916.xlsx"
DEFAULT_OUTPUT = ROOT / "artifacts"


def build(old_workbook: Path, new_workbook: Path, output_dir: Path, n_folds: int) -> dict:
    cota, cota_stats = load_cota(old_workbook)
    preamble, preamble_stats = load_preamble(new_workbook)

    write_jsonl(output_dir / "normalized" / "cota_patients.jsonl", cota)
    write_jsonl(output_dir / "normalized" / "preamble_patients.jsonl", preamble)

    folds = make_folds(cota, n_folds=n_folds)
    write_fold_manifest(output_dir / "cv" / f"cota_{n_folds}fold_manifest.csv", folds)
    truth_by_fold: dict[str, dict[str, int]] = {}
    for fold in range(n_folds):
        truth_by_fold[str(fold)] = dict(sorted(Counter(
            str(row["reviewer_consensus_lot"])
            for row in folds if row["fold"] == fold
        ).items()))

    consensus_counts = Counter(record.reviewer_consensus_status for record in cota)
    report = {
        "schema_version": "1.0.0",
        "sources": {
            "cota_old_workbook": {
                **cota_stats.public_dict(),
                "reviewer_consensus_status": dict(sorted(consensus_counts.items())),
                "malformed_treatments": malformed_summary(cota),
                "duplicates": duplicate_summary(cota),
            },
            "preamble_new_workbook": {
                **preamble_stats.public_dict(),
                "reviewer_consensus_status": dict(sorted(Counter(
                    record.reviewer_consensus_status for record in preamble
                ).items())),
                "malformed_treatments": malformed_summary(preamble),
                "duplicates": duplicate_summary(preamble),
            },
        },
        "cross_workbook_overlap": overlap_summary(cota, preamble),
        "cross_validation": {
            "algorithm": "Consensus-LOT-stratified SHA-256 ordering with round-robin assignment",
            "fold_count": n_folds,
            "patient_count": len(folds),
            "patients_per_fold": dict(sorted(Counter(row["fold"] for row in folds).items())),
            "reviewer_consensus_lot_distribution_by_fold": truth_by_fold,
        },
    }
    write_json(output_dir / "reports" / "data_quality_report.json", report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-workbook", type=Path, default=DEFAULT_OLD)
    parser.add_argument("--new-workbook", type=Path, default=DEFAULT_NEW)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--folds", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = build(args.old_workbook, args.new_workbook, args.output_dir, args.folds)
    print(f"Prepared {report['cross_validation']['patient_count']} adjudicated COTA patients.")
    print(f"Outputs written under {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
