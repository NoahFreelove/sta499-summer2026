from __future__ import annotations

import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from openpyxl import Workbook


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "py"))

from lot_data import (  # noqa: E402
    PatientRecord,
    TreatmentEvent,
    load_cota,
    make_folds,
    parse_treatment,
)


class TreatmentParserTests(unittest.TestCase):
    def test_normalizes_sets_and_removes_sct(self) -> None:
        drugs, flags = parse_treatment(
            "[ Bortezomib, DEX ], [autologous SCT]", bracketed=True
        )
        self.assertEqual(drugs, frozenset({"bortezomib", "dex"}))
        self.assertEqual(flags, [])

    def test_reports_unbalanced_and_empty_treatments(self) -> None:
        drugs, flags = parse_treatment("[", bracketed=True)
        self.assertEqual(drugs, frozenset())
        self.assertEqual(flags, ["empty_drug_set", "unbalanced_brackets"])

    def test_reconstructs_continuation_before_parsing(self) -> None:
        workbook = Workbook()
        sheet = workbook.active
        sheet.title = "Cota"
        sheet.append([
            "cpid", "line_of_therapy_c", "line_of_therapy_name",
            "date_start_line_of_therapy", "date_end_line_of_therapy",
            "Alpesh 1 LoT", "Alberto LOT ",
        ])
        sheet.append(["raw-1", 1, "[Drug A, ", "2020-01-01", "2020-02-01", 2, 2])
        sheet.append([None, None, "Drug B]", None, None, None, None])
        sheet.append(["raw-1", 2, "[Drug C]", "2020-03-01", "2020-04-01", 2, 2])

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.xlsx"
            workbook.save(path)
            records, stats = load_cota(path)

        self.assertEqual(len(records), 1)
        self.assertEqual(stats.continuation_rows, 1)
        self.assertEqual(records[0].trajectory[0].drugs, frozenset({"drug a", "drug b"}))
        self.assertEqual(records[0].trajectory[0].continuation_row_count, 1)
        self.assertEqual(records[0].reviewer_consensus_lot, 2)
        self.assertNotIn("raw-1", str(records[0].public_dict()))


class FoldTests(unittest.TestCase):
    @staticmethod
    def records(count: int = 136) -> list[PatientRecord]:
        return [
            PatientRecord(
                source="cota",
                raw_patient_id=f"private-{index}",
                trajectory=[TreatmentEvent(1, 1, frozenset({"drug"}))],
                reviewer_a_lot=(index % 8) + 1,
                reviewer_b_lot=(index % 8) + 1,
                reviewer_consensus_lot=(index % 8) + 1,
                reviewer_consensus_status="agreed",
                cota_lot=(index % 8) + 1,
            )
            for index in range(count)
        ]

    def test_folds_are_deterministic_patient_level_and_balanced(self) -> None:
        records = self.records()
        first = make_folds(records, 5)
        second = make_folds(list(reversed(records)), 5)
        self.assertEqual(first, second)
        self.assertEqual(len({row["patient_key"] for row in first}), 136)
        self.assertEqual(sorted(Counter(row["fold"] for row in first).values()), [27, 27, 27, 27, 28])
        for lot in range(1, 9):
            counts = Counter(
                row["fold"] for row in first if row["reviewer_consensus_lot"] == lot
            )
            self.assertLessEqual(max(counts.values()) - min(counts.values()), 1)


if __name__ == "__main__":
    unittest.main()
