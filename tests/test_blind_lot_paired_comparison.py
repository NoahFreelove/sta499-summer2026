from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "py"))

from blind_lot.paired_comparison import (  # noqa: E402
    build_comparison,
    load_aligned_runs,
    paired_bootstrap,
    write_comparison,
)
from blind_lot.evaluation import normalize_joined_record  # noqa: E402


class PairedComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        fixture = json.loads((ROOT / "tests/fixtures/paired_rows.json").read_text())["rows"]
        self.rows3, self.rows5 = [], []
        for number, item in enumerate(fixture):
            common = {"case_key": f"SYNTH-{number:02d}", "fold": item["fold"],
                      "reviewer_lot": item["reviewer"], "algorithm_lot": item["algorithm"],
                      "cota_lot": item["cota"]}
            self.rows3.append({**common, "ai_lot": item["k3_ai"], "ai_abstained": item["k3_abstain"]})
            self.rows5.append({**common, "ai_lot": item["k5_ai"], "ai_abstained": item["k5_abstain"]})
        self.run3, self.run5 = self.root / "k3", self.root / "k5"
        self._write_run(self.run3, 3, self.rows3)
        self._write_run(self.run5, 5, self.rows5)

    def tearDown(self) -> None:
        self.temp.cleanup()

    @staticmethod
    def _metadata(k: int) -> dict:
        return {"retrieval_k": k, "provider": "openai", "model": "model",
                "reasoning_effort": "medium", "temperature": None, "prompt_version": "p",
                "knowledge_version": "q", "input_artifact_sha256": {"x": "y"},
                "random_seed": 4992026, "bootstrap_replicates": 2000,
                "evaluation_version": "blind-lot-evaluation-1.2.0", "folds": [0, 1, 2, 3, 4]}

    def _write_run(self, path: Path, k: int, rows: list[dict]) -> None:
        restricted = path / "restricted"
        restricted.mkdir(parents=True, exist_ok=True)
        (restricted / "experiment_metadata.json").write_text(json.dumps(self._metadata(k)))
        (restricted / "joined_evaluation.reaggregated.jsonl").write_text(
            "".join(json.dumps(row) + "\n" for row in rows)
        )

    def test_categories_cover_requested_state_changes(self) -> None:
        metadata, aligned = load_aligned_runs(self.run3, self.run5)
        report, restricted = build_comparison(metadata, aligned, replicates=20)
        categories = report["paired_categories"]
        self.assertGreater(categories["abstention_and_usability"]["both_abstain"], 0)
        self.assertGreater(categories["abstention_and_usability"]["k3_abstains_k5_usable"], 0)
        self.assertGreater(categories["abstention_and_usability"]["k3_usable_k5_abstains"], 0)
        self.assertGreater(categories["abstention_and_usability"]["both_usable_same_total"], 0)
        self.assertGreater(categories["abstention_and_usability"]["both_usable_different_totals"], 0)
        self.assertGreater(categories["generated_total_correctness"]["correct_only_k3"], 0)
        self.assertGreater(categories["generated_total_correctness"]["correct_only_k5"], 0)
        self.assertGreater(categories["three_way_acceptance"]["accepted_by_both"], 0)
        self.assertGreater(categories["three_way_acceptance"]["accepted_only_k3"], 0)
        self.assertGreater(categories["three_way_acceptance"]["accepted_only_k5"], 0)
        self.assertGreater(categories["three_way_false_acceptance"]["false_accepted_by_both"], 0)
        self.assertGreater(categories["three_way_false_acceptance"]["false_accepted_only_k3"], 0)
        self.assertGreater(categories["three_way_false_acceptance"]["false_accepted_only_k5"], 0)
        self.assertTrue(all("case_key" in row for row in restricted))

    def test_missing_and_duplicate_case_keys_are_rejected(self) -> None:
        self._write_run(self.run5, 5, self.rows5[:-1])
        with self.assertRaisesRegex(ValueError, "identical case"):
            load_aligned_runs(self.run3, self.run5)
        self._write_run(self.run5, 5, self.rows5 + [self.rows5[0]])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            load_aligned_runs(self.run3, self.run5)

    def test_bootstrap_is_deterministic_and_fold_exclusions_are_correct(self) -> None:
        metadata, aligned = load_aligned_runs(self.run3, self.run5)
        self.assertEqual(paired_bootstrap(aligned, seed=7, replicates=30),
                         paired_bootstrap(aligned, seed=7, replicates=30))
        report, _ = build_comparison(metadata, aligned, seed=7, replicates=10)
        analyses = report["bootstrap"]["analyses"]
        self.assertEqual(analyses["leave_fold_0_out"]["patient_count"],
                         len([row for row in self.rows3 if row["fold"] != 0]))
        self.assertEqual(analyses["leave_fold_3_out"]["patient_count"],
                         len([row for row in self.rows3 if row["fold"] != 3]))

    def test_undefined_accuracy_interval_is_null(self) -> None:
        rows = [(normalize_joined_record({**self.rows3[0], "ai_abstained": True}),
                 normalize_joined_record({**self.rows5[0], "ai_abstained": True}))]
        result = paired_bootstrap(rows, replicates=10)
        interval = result["selective_ai_accuracy"]["percentile_95_ci"]
        self.assertEqual(interval, {"lower": None, "upper": None, "replicates_used": 0})

    def test_public_is_aggregate_only_and_restricted_is_pseudonymous(self) -> None:
        paths = write_comparison(self.run3, self.run5, self.root / "comparison", replicates=10)
        public = paths["json"].read_text()
        self.assertNotIn("SYNTH-", public)
        self.assertNotIn(str(self.root), public)
        restricted = paths["restricted"].read_text()
        self.assertIn("SYNTH-", restricted)
        for marker in ("COTAOLD-", "COTANEW-"):
            self.assertNotIn(marker, restricted)


if __name__ == "__main__":
    unittest.main()
