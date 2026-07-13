"""Prepare private evaluation data, blind inputs, quality reports, and folds."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from lot_data import (
    ADJUDICATED_SOURCE,
    assert_aggregate_only_report,
    assign_grouped_folds,
    build_exclusion_groups,
    duplicate_summary,
    find_project_root,
    load_adjudicated_cota,
    load_new_cota,
    malformed_summary,
    overlap_public_summary,
    schema_difference_summary,
    write_blind_jsonl,
    write_evaluation_jsonl,
    write_fold_manifest,
    write_json,
    write_jsonl_rows,
)


def resolve_path(value: Path | None, root: Path, default: str) -> Path:
    path = value if value is not None else Path(default)
    return path if path.is_absolute() else (root / path).resolve()


def build(
    adjudicated_workbook: Path,
    new_workbook: Path,
    output_dir: Path,
    n_folds: int,
    adjudicated_sheet: str = "Cota",
    new_sheet: str | None = None,
) -> dict:
    old, old_stats = load_adjudicated_cota(adjudicated_workbook, adjudicated_sheet)
    new, new_stats = load_new_cota(new_workbook, new_sheet)
    combined = old + new

    restricted = output_dir / "restricted"
    blind = output_dir / "blind"
    public = output_dir / "public"

    write_evaluation_jsonl(restricted / "evaluation" / "cota_adjudicated.jsonl", old)
    write_evaluation_jsonl(restricted / "evaluation" / "cota_unadjudicated.jsonl", new)
    write_blind_jsonl(blind / "cota_adjudicated.jsonl", old)
    write_blind_jsonl(blind / "cota_unadjudicated.jsonl", new)

    exclusion_rows = build_exclusion_groups(combined)
    folded_rows = assign_grouped_folds(combined, exclusion_rows, n_folds=n_folds)
    write_jsonl_rows(
        restricted / "overlap" / "exclusion_groups.jsonl",
        folded_rows,
        "patient_key",
    )
    write_fold_manifest(
        restricted / "cv" / f"cota_adjudicated_{n_folds}fold_manifest.csv",
        folded_rows,
    )

    old_fold_counts = Counter(
        row["fold"] for row in folded_rows if row["source"] == ADJUDICATED_SOURCE
    )
    consensus_status = Counter(record.reviewer_consensus_status for record in old)
    new_consensus_status = Counter(record.reviewer_consensus_status for record in new)
    report = {
        "schema_version": "2.0.0",
        "report_scope": "aggregate_only",
        "sources": {
            "cota_adjudicated": {
                **old_stats.public_dict(),
                "reviewer_consensus_status": dict(sorted(consensus_status.items())),
                "malformed_records": malformed_summary(old),
                "duplicates": duplicate_summary(old),
            },
            "cota_unadjudicated": {
                **new_stats.public_dict(),
                "reviewer_consensus_status": dict(sorted(new_consensus_status.items())),
                "malformed_records": malformed_summary(new),
                "duplicates": duplicate_summary(new),
            },
        },
        "schema_differences": schema_difference_summary(old_stats, new_stats),
        "cross_workbook_overlap": overlap_public_summary(old, new, exclusion_rows),
        "cross_validation": {
            "algorithm": "connected exclusion groups over raw-ID, exact, and collapsed signatures; deterministic group-aware assignment",
            "fold_count": n_folds,
            "adjudicated_patient_count": len(old),
            "adjudicated_patients_per_fold": {
                str(fold): old_fold_counts[fold] for fold in range(n_folds)
            },
        },
        "blind_input": {
            "boundary_policy": "bracketed regimen groups flattened in source order",
            "vendor_line_numbers_included": False,
            "vendor_line_dates_included": False,
            "reviewer_labels_included": False,
            "cota_lot_included": False,
            "unavoidable_boundary_information": (
                "The workbook orders bracketed regimen groups but provides dates only at "
                "vendor-line granularity. Blind inputs retain group order and omit those dates; "
                "exact within-line group timing cannot be recovered from the source."
            ),
        },
    }
    assert_aggregate_only_report(report)
    write_json(public / "data_quality_report.json", report)

    # Remove v1 outputs that could otherwise be mistaken for the current contract.
    for obsolete in (
        output_dir / "normalized" / "cota_patients.jsonl",
        output_dir / "normalized" / "preamble_patients.jsonl",
        output_dir / "cv" / "cota_5fold_manifest.csv",
        output_dir / "reports" / "data_quality_report.json",
        output_dir / "reports" / "baseline_evaluation.json",
    ):
        obsolete.unlink(missing_ok=True)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path)
    parser.add_argument("--adjudicated-workbook", type=Path)
    parser.add_argument("--new-workbook", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--adjudicated-sheet", default="Cota")
    parser.add_argument("--new-sheet")
    parser.add_argument("--folds", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = find_project_root(args.project_root)
    adjudicated = resolve_path(
        args.adjudicated_workbook, root, "data/LoT Adjudication Datasets.xlsx"
    )
    new = resolve_path(args.new_workbook, root, "data/new_cota_data.xlsx")
    output = resolve_path(args.output_dir, root, "artifacts")
    report = build(
        adjudicated,
        new,
        output,
        args.folds,
        args.adjudicated_sheet,
        args.new_sheet,
    )
    old_count = report["sources"]["cota_adjudicated"]["patient_count"]
    new_count = report["sources"]["cota_unadjudicated"]["patient_count"]
    print(f"Prepared {old_count} adjudicated and {new_count} unadjudicated COTA patients.")
    print(f"Aggregate public and restricted outputs written under {output}")


if __name__ == "__main__":
    main()
