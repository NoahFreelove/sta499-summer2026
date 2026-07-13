from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from openpyxl import Workbook


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "py"))

from evaluate_baseline import evaluate  # noqa: E402
from lot_data import (  # noqa: E402
    ADJUDICATED_SOURCE,
    UNADJUDICATED_SOURCE,
    PatientRecord,
    TreatmentEvent,
    assert_aggregate_only_report,
    assign_grouped_folds,
    build_exclusion_groups,
    load_adjudicated_cota,
    load_new_cota,
    parse_treatment_groups,
)


OLD_HEADERS = [
    "cpid", "diag_dt", "line_of_therapy_c", "line_of_therapy_name",
    "date_start_line_of_therapy", "date_end_line_of_therapy",
    "Alpesh 1 LoT", "Alberto LOT ",
]
NEW_HEADERS = [
    "cpid", "diag_dt", "deathfl", "dthdt_c", "date_of_death_imp",
    "line_of_therapy_c", "line_of_therapy_name", "discontinue_reason",
    "date_start_line_of_therapy", "date_start_line_of_therapy_imp",
    "date_end_line_of_therapy", "date_end_line_of_therapy_imp", "aval", "avaldt",
]


def save_workbook(path: Path, headers: list[str], rows: list[list[object]], title: str) -> None:
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = title
    sheet.append(headers)
    for row in rows:
        sheet.append(row)
    workbook.save(path)


def record(
    raw_id: str,
    sequence: list[set[str]],
    *,
    source: str = ADJUDICATED_SOURCE,
    truth: int | None = 1,
) -> PatientRecord:
    return PatientRecord(
        source=source,
        raw_patient_id=raw_id,
        trajectory=[
            TreatmentEvent(order=index, drugs=frozenset(drugs), source_line_number=index)
            for index, drugs in enumerate(sequence, start=1)
        ],
        diagnosis_date="2020-01-01",
        reviewer_a_lot=truth,
        reviewer_b_lot=truth,
        reviewer_consensus_lot=truth,
        reviewer_consensus_status="agreed" if truth is not None else "unadjudicated",
        cota_lot=len(sequence),
    )


class ParserTests(unittest.TestCase):
    def test_flattens_bracketed_groups_in_order_and_removes_sct(self) -> None:
        groups, flags, removed = parse_treatment_groups(
            "[Bortezomib, DEX], [autologous SCT], [Lenalidomide]"
        )
        self.assertEqual(groups, [
            frozenset({"bortezomib", "dex"}),
            frozenset({"lenalidomide"}),
        ])
        self.assertEqual(flags, [])
        self.assertEqual(removed, 1)

    def test_reports_malformed_treatment(self) -> None:
        groups, flags, _ = parse_treatment_groups("[")
        self.assertEqual(groups, [])
        self.assertEqual(flags, ["empty_drug_set", "empty_regimen_group", "unbalanced_brackets"])

    def test_adjudicated_continuation_rows_are_reconstructed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "old.xlsx"
            save_workbook(path, OLD_HEADERS, [
                ["raw-1", "2020-01-01", 1, "[Drug A, ", "2020-02-01", "2020-03-01", 2, 2],
                [None, None, None, "Drug B], [Drug C]", None, None, None, None],
                ["raw-1", None, 2, "[Drug D]", "2020-04-01", "2020-05-01", 2, 2],
            ], "Cota")
            records, stats = load_adjudicated_cota(path)

        self.assertEqual(stats.continuation_rows, 1)
        self.assertEqual([event.drugs for event in records[0].trajectory], [
            frozenset({"drug a", "drug b"}),
            frozenset({"drug c"}),
            frozenset({"drug d"}),
        ])
        self.assertEqual(records[0].reviewer_consensus_lot, 2)

    def test_new_cota_schema_and_blank_id_new_lines(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "new.xlsx"
            save_workbook(path, NEW_HEADERS, [
                ["new-1", "2021-01-01", "N", None, None, 1, "[Drug A], [Drug B, ", None, "2021-02-01", 0, "2021-03-01", 0, None, None],
                [None, None, None, None, None, None, "Drug C]", None, None, None, None, None, None, None],
                [None, None, None, None, None, 2, "[Drug D]", None, "2021-04-01", 0, "2021-05-01", 0, None, None],
            ], "ActualSheet")
            records, stats = load_new_cota(path)

        self.assertEqual(stats.sheet_names, ["ActualSheet"])
        self.assertEqual(stats.selected_sheet, "ActualSheet")
        self.assertEqual(stats.patient_count, 1)
        self.assertEqual(stats.vendor_line_count, 2)
        self.assertEqual(stats.continuation_rows, 1)
        self.assertEqual(stats.same_patient_new_line_rows, 1)
        self.assertEqual(records[0].reviewer_consensus_lot, None)
        self.assertEqual(records[0].reviewer_consensus_status, "unadjudicated")
        self.assertEqual(records[0].cota_lot, 2)
        self.assertEqual([event.drugs for event in records[0].trajectory], [
            frozenset({"drug a"}),
            frozenset({"drug b", "drug c"}),
            frozenset({"drug d"}),
        ])


class BlindInputTests(unittest.TestCase):
    def test_forbidden_fields_never_appear_anywhere_in_blind_serialization(self) -> None:
        private = record("secret-id", [{"a"}, {"b"}])
        blind = private.blind_dict()
        serialized = json.dumps(blind, sort_keys=True)
        forbidden = {
            "patient_key", "source", "source_line_number", "source_group_index",
            "source_line_start_date", "source_line_end_date", "cota_lot",
            "reviewer_lot", "reviewer_a", "reviewer_b", "consensus",
            "consensus_status", "algorithm_lot", "algorithm_flags",
            "quality_metadata", "line_start_date", "line_end_date",
        }

        def keys(value: object) -> set[str]:
            if isinstance(value, dict):
                return set(value) | set().union(*(keys(item) for item in value.values()))
            if isinstance(value, list):
                return set().union(*(keys(item) for item in value), set())
            return set()

        self.assertFalse(keys(blind) & forbidden)
        for field in forbidden:
            self.assertNotIn(f'"{field}"', serialized)
        self.assertNotIn("secret-id", serialized)
        self.assertEqual(set(blind), {"schema_version", "case_key", "trajectory", "context"})


class LeakageSafeFoldTests(unittest.TestCase):
    def test_exact_and_collapsed_duplicates_remain_in_one_fold(self) -> None:
        records = [
            record("one", [{"a"}, {"b"}], truth=2),
            record("two", [{"a"}, {"b"}], truth=2),
            record("three", [{"a"}, {"a"}, {"b"}], truth=2),
            record("four", [{"c"}], truth=1),
        ]
        groups = build_exclusion_groups(records)
        folded = assign_grouped_folds(records, groups, 3)
        duplicate_folds = {
            row["fold"] for row in folded
            if row["patient_key"] in {record.key for record in records[:3]}
        }
        self.assertEqual(len(duplicate_folds), 1)

    def test_cross_workbook_id_overlap_is_one_exclusion_group(self) -> None:
        old = record("shared-private-id", [{"a"}], truth=1)
        new = record(
            "shared-private-id", [{"different"}], source=UNADJUDICATED_SOURCE, truth=None
        )
        records = [old, new]
        groups = build_exclusion_groups(records)
        self.assertEqual(len({row["exclusion_group"] for row in groups}), 1)
        self.assertIn("cross_workbook_patient_id", groups[0]["group_reasons"])
        self.assertNotIn("shared-private-id", json.dumps(groups))


class PublicReportTests(unittest.TestCase):
    def test_baseline_public_report_is_aggregate_only(self) -> None:
        evaluation_record = record("private", [{"a"}]).evaluation_dict()
        public, private_rows = evaluate([evaluation_record])
        self.assertEqual(len(private_rows), 1)
        self.assertEqual(public["report_scope"], "aggregate_only")
        serialized = json.dumps(public)
        self.assertNotIn("patient_key", serialized)
        self.assertNotIn("CASE-", serialized)
        self.assertNotIn("both_agree_but_wrong_cases", serialized)

    def test_public_report_guard_rejects_patient_level_content(self) -> None:
        assert_aggregate_only_report({
            "report_scope": "aggregate_only",
            "sources": {"new": {"patient_count": 10}},
        })
        with self.assertRaises(ValueError):
            assert_aggregate_only_report({
                "report_scope": "aggregate_only",
                "patient_key": "CASE-should-not-be-public",
            })


if __name__ == "__main__":
    unittest.main()
