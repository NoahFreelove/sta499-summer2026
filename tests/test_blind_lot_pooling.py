from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "py"))

from blind_lot.pooling import pool_runs, sha256_file  # noqa: E402


class PoolingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.runs = self.root / "runs"
        self.output = self.root / "pooled"
        self.manifest_path = self.root / "manifest.json"
        self.manifest = {
            "schema_version": "blind-lot-pooling-manifest-1.0.0",
            "pool_id": "synthetic", "expected_folds": [0, 1, 2, 3, 4],
            "expected_patient_count": 5, "provider": "openai", "model": "model",
            "reasoning_effort": "medium", "temperature": None, "prompt_version": "prompt",
            "knowledge_version": "knowledge", "evaluation_version": "blind-lot-evaluation-1.2.0",
            "bootstrap_seed": 4992026, "bootstrap_replicates": 5, "source_runs": [],
        }
        for k in (3, 5):
            for fold in range(5):
                run_id = f"run-k{k}-f{fold}"
                self.manifest["source_runs"].append({"retrieval_k": k, "fold": fold, "run_id": run_id})
                self._write_run(run_id, k, fold)
        self._write_manifest()

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _write_manifest(self) -> None:
        self.manifest_path.write_text(json.dumps(self.manifest), encoding="utf-8")

    def _metadata(self, run_id: str, k: int, fold: int) -> dict:
        return {
            "run_id": run_id, "retrieval_k": k, "folds": [fold], "provider": "openai",
            "model": "model", "reasoning_effort": "medium", "temperature": None,
            "prompt_version": "prompt", "knowledge_version": "knowledge",
            "input_artifact_sha256": {"input": "abc"}, "random_seed": 4992026,
            "bootstrap_replicates": 5, "evaluation_version": "blind-lot-evaluation-1.2.0",
            "limit": None,
        }

    def _write_run(self, run_id: str, k: int, fold: int) -> None:
        restricted = self.runs / run_id / "restricted"
        restricted.mkdir(parents=True)
        row = {"case_key": f"SYNTH-{fold}", "fold": fold, "reviewer_lot": 2,
               "algorithm_lot": 2, "cota_lot": 2, "ai_lot": 2, "ai_abstained": False}
        (restricted / "experiment_metadata.json").write_text(
            json.dumps(self._metadata(run_id, k, fold)), encoding="utf-8"
        )
        for name, value in (("joined_evaluation.jsonl", row),
                            ("blind_predictions.jsonl", {"case_key": row["case_key"]}),
                            ("retrieval_debug.jsonl", {"case_key": row["case_key"]})):
            (restricted / name).write_text(json.dumps(value) + "\n", encoding="utf-8")

    def _metadata_path(self, run_id: str) -> Path:
        return self.runs / run_id / "restricted" / "experiment_metadata.json"

    def _change_metadata(self, run_id: str, field: str, value: object) -> None:
        path = self._metadata_path(run_id)
        metadata = json.loads(path.read_text())
        metadata[field] = value
        path.write_text(json.dumps(metadata), encoding="utf-8")

    def test_valid_disjoint_folds_pool_successfully_and_order_deterministic(self) -> None:
        outputs = pool_runs(self.manifest_path, self.runs, self.output)
        self.assertEqual(set(outputs), {3, 5})
        rows = (outputs[3] / "restricted/joined_evaluation.jsonl").read_text().splitlines()
        self.assertEqual([json.loads(row)["case_key"] for row in rows], [f"SYNTH-{i}" for i in range(5)])

    def test_missing_fold_is_rejected(self) -> None:
        self.manifest["source_runs"].pop()
        self._write_manifest()
        with self.assertRaisesRegex(ValueError, "missing"):
            pool_runs(self.manifest_path, self.runs, self.output)

    def test_duplicate_pair_is_rejected(self) -> None:
        self.manifest["source_runs"].append(dict(self.manifest["source_runs"][0]))
        self._write_manifest()
        with self.assertRaisesRegex(ValueError, "duplicate"):
            pool_runs(self.manifest_path, self.runs, self.output)

    def test_duplicate_case_across_folds_is_rejected(self) -> None:
        path = self.runs / "run-k3-f1/restricted/joined_evaluation.jsonl"
        row = json.loads(path.read_text())
        row["case_key"] = "SYNTH-0"
        path.write_text(json.dumps(row) + "\n")
        with self.assertRaisesRegex(ValueError, "duplicate case"):
            pool_runs(self.manifest_path, self.runs, self.output)

    def test_incorrect_row_fold_is_rejected(self) -> None:
        path = self.runs / "run-k3-f1/restricted/joined_evaluation.jsonl"
        row = json.loads(path.read_text())
        row["fold"] = 0
        path.write_text(json.dumps(row) + "\n")
        with self.assertRaisesRegex(ValueError, "incorrect fold"):
            pool_runs(self.manifest_path, self.runs, self.output)

    def test_limited_and_dry_runs_are_rejected(self) -> None:
        for field, value, message in (("limit", 1, "limited"), ("dry_run", True, "dry")):
            with self.subTest(field=field):
                self._change_metadata("run-k3-f0", field, value)
                with self.assertRaisesRegex(ValueError, message):
                    pool_runs(self.manifest_path, self.runs, self.output)
                self._write_run_fresh("run-k3-f0", 3, 0)

    def _write_run_fresh(self, run_id: str, k: int, fold: int) -> None:
        import shutil
        shutil.rmtree(self.runs / run_id)
        self._write_run(run_id, k, fold)

    def test_frozen_metadata_mismatches_are_rejected(self) -> None:
        cases = (("model", "other"), ("prompt_version", "other"),
                 ("knowledge_version", "other"), ("input_artifact_sha256", {"input": "no"}),
                 ("bootstrap_replicates", 6), ("evaluation_version", "other"))
        for field, value in cases:
            with self.subTest(field=field):
                self._change_metadata("run-k3-f1", field, value)
                with self.assertRaises(ValueError):
                    pool_runs(self.manifest_path, self.runs, self.output)
                self._write_run_fresh("run-k3-f1", 3, 1)

    def test_cohort_mismatch_is_rejected(self) -> None:
        for name in ("joined_evaluation.jsonl", "blind_predictions.jsonl", "retrieval_debug.jsonl"):
            path = self.runs / "run-k5-f4/restricted" / name
            row = json.loads(path.read_text())
            row["case_key"] = "SYNTH-DIFFERENT"
            path.write_text(json.dumps(row) + "\n")
        with self.assertRaisesRegex(ValueError, "different case-key"):
            pool_runs(self.manifest_path, self.runs, self.output)

    def test_sources_unchanged_and_no_api_key_needed(self) -> None:
        paths = list(self.runs.glob("*/restricted/*"))
        before = {path: sha256_file(path) for path in paths}
        with patch.dict(os.environ, {}, clear=True):
            pool_runs(self.manifest_path, self.runs, self.output)
        self.assertEqual(before, {path: sha256_file(path) for path in paths})

    def test_provider_and_retrieval_are_never_executed(self) -> None:
        with patch("blind_lot.providers.OpenAIResponsesProvider.__init__", side_effect=AssertionError), \
             patch("blind_lot.retrieval.LocalRetriever.retrieve", side_effect=AssertionError):
            pool_runs(self.manifest_path, self.runs, self.output)

    def test_repeated_public_outputs_are_byte_identical_and_safe(self) -> None:
        first = pool_runs(self.manifest_path, self.runs, self.output)
        content = (first[3] / "public/aggregate_evaluation.reaggregated.json").read_bytes()
        second = pool_runs(self.manifest_path, self.runs, self.root / "pooled-again")
        other = (second[3] / "public/aggregate_evaluation.reaggregated.json").read_bytes()
        self.assertEqual(content, other)
        text = content.decode()
        for marker in ("CASE-", "COTAOLD-", "COTANEW-", "EXCL-", str(self.root)):
            self.assertNotIn(marker, text)


if __name__ == "__main__":
    unittest.main()
