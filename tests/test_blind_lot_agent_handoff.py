from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "py"))

from blind_lot.agent_handoff import (  # noqa: E402
    HANDOFF_SCHEMA_VERSION,
    SHARD_OUTPUT_SCHEMA_VERSION,
    export_worker_bundle,
    import_worker_bundle_outputs,
    validate_prediction_shards,
)
from compare_blind_lot_model_runs import build_leaderboard  # noqa: E402


def case_id(number: int) -> str:
    return f"CASE-{number:020x}"


def valid_response(*, cited: list[str] | None = None) -> dict:
    return {
        "schema_version": "1.0.0",
        "ai_lot_count": 2,
        "abstained": False,
        "transitions": [
            {
                "transition_index": 1,
                "decision": "NEW_LINE",
                "reason_codes": ["SYNTHETIC"],
                "evidence_strength": "moderate",
                "retrieved_case_ids": cited or [],
            }
        ],
        "warnings": [],
    }


class AgentHandoffValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.run_dir = Path(self.temporary.name)
        packet_dir = self.run_dir / "restricted" / "agent_handoff" / "packets"
        output_dir = self.run_dir / "restricted" / "agent_handoff" / "outputs"
        packet_dir.mkdir(parents=True)
        output_dir.mkdir()
        self.packet_path = packet_dir / "shard-001.json"
        self.output_path = output_dir / "shard-001.json"
        self.case_key = case_id(1)
        self.retrieved_key = case_id(2)
        self.packet = {
            "schema_version": HANDOFF_SCHEMA_VERSION,
            "run_id": "synthetic-run",
            "shard_id": "shard-001",
            "cases": [
                {
                    "case_key": self.case_key,
                    "fold": 0,
                    "regimen_event_count": 2,
                    "permitted_retrieved_case_ids": [self.retrieved_key],
                    "prompt_sha256": "prompt-hash",
                }
            ],
        }
        self.packet_path.write_text(
            json.dumps(self.packet, sort_keys=True), encoding="utf-8"
        )
        packet_hash = hashlib.sha256(self.packet_path.read_bytes()).hexdigest()
        self.manifest = {
            "schema_version": HANDOFF_SCHEMA_VERSION,
            "run_id": "synthetic-run",
            "provider": "codex-subagents",
            "model": "synthetic-model",
            "reasoning_effort": "high",
            "run_purpose": "smoke",
            "repeat_index": 1,
            "case_count": 1,
            "shard_count": 1,
            "shards": [
                {
                    "shard_id": "shard-001",
                    "packet": str(self.packet_path.relative_to(self.run_dir)),
                    "output": str(self.output_path.relative_to(self.run_dir)),
                    "case_count": 1,
                    "case_keys": [self.case_key],
                    "packet_sha256": packet_hash,
                }
            ],
        }
        manifest_path = self.run_dir / "restricted" / "agent_handoff" / "manifest.json"
        manifest_path.write_text(json.dumps(self.manifest), encoding="utf-8")
        (self.run_dir / "restricted" / "experiment_metadata.json").write_text(
            json.dumps({"run_id": "synthetic-run"}), encoding="utf-8"
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def output_value(
        self,
        response: dict,
        case_key: str | None = None,
        *,
        provenance: dict | None = None,
    ) -> dict:
        value = {
            "schema_version": SHARD_OUTPUT_SCHEMA_VERSION,
            "run_id": "synthetic-run",
            "shard_id": "shard-001",
            "predictions": [
                {
                    "case_key": case_key or self.case_key,
                    "response": response,
                }
            ],
        }
        if provenance is not None:
            value["worker_provenance"] = provenance
        return value

    def write_output(
        self,
        response: dict,
        case_key: str | None = None,
        *,
        provenance: dict | None = None,
    ) -> None:
        self.output_path.write_text(
            json.dumps(self.output_value(response, case_key, provenance=provenance)),
            encoding="utf-8",
        )

    def test_valid_shard_is_normalized(self) -> None:
        self.write_output(valid_response(cited=[self.retrieved_key]))
        rows = validate_prediction_shards(self.run_dir, self.manifest)
        self.assertEqual(rows[0]["case_key"], self.case_key)
        self.assertEqual(rows[0]["agent_shard_id"], "shard-001")
        self.assertEqual(rows[0]["response"]["ai_lot_count"], 2)

    def test_missing_case_is_rejected(self) -> None:
        self.write_output(valid_response(), case_key=case_id(99))
        with self.assertRaisesRegex(ValueError, "coverage mismatch"):
            validate_prediction_shards(self.run_dir, self.manifest)

    def test_unpermitted_retrieval_citation_is_rejected(self) -> None:
        self.write_output(valid_response(cited=[case_id(99)]))
        with self.assertRaisesRegex(ValueError, "was not supplied"):
            validate_prediction_shards(self.run_dir, self.manifest)

    def test_changed_packet_is_rejected(self) -> None:
        self.write_output(valid_response())
        self.packet["cases"][0]["fold"] = 4
        self.packet_path.write_text(json.dumps(self.packet), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "changed after preparation"):
            validate_prediction_shards(self.run_dir, self.manifest)

    def test_development_run_requires_matching_worker_provenance(self) -> None:
        self.manifest["run_purpose"] = "development"
        self.write_output(valid_response())
        with self.assertRaisesRegex(ValueError, "worker provenance is required"):
            validate_prediction_shards(self.run_dir, self.manifest)
        self.write_output(
            valid_response(),
            provenance={
                "configured_model": "synthetic-model",
                "configured_reasoning_effort": "high",
                "attempt": 1,
                "agent_id": "agent-123",
            },
        )
        rows = validate_prediction_shards(self.run_dir, self.manifest)
        self.assertEqual(rows[0]["worker_provenance"]["agent_id"], "agent-123")

    def test_external_worker_bundle_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as bundle_parent:
            bundle_dir = Path(bundle_parent) / "worker-bundle"
            exported = export_worker_bundle(
                self.run_dir,
                bundle_dir,
                project_root=self.run_dir,
                orchestrator="codex",
            )
            self.assertTrue(Path(exported["manifest"]).is_file())
            bundle_manifest = json.loads(Path(exported["manifest"]).read_text())
            worker_output = bundle_dir / bundle_manifest["shards"][0]["output"]
            worker_output.write_text(
                json.dumps(
                    self.output_value(
                        valid_response(),
                        provenance={
                            "configured_model": "synthetic-model",
                            "configured_reasoning_effort": "high",
                            "attempt": 1,
                            "agent_id": "agent-123",
                        },
                    )
                ),
                encoding="utf-8",
            )
            imported = import_worker_bundle_outputs(self.run_dir, bundle_dir)
            self.assertTrue(self.output_path.is_file())
            self.assertIn("shard-001", imported["output_sha256"])

    def test_claude_worker_bundle_pins_model_and_effort(self) -> None:
        with tempfile.TemporaryDirectory() as bundle_parent:
            bundle_dir = Path(bundle_parent) / "worker-bundle"
            export_worker_bundle(
                self.run_dir,
                bundle_dir,
                project_root=self.run_dir,
                orchestrator="claude-code",
            )
            worker = (
                bundle_dir / ".claude" / "agents" / "blind-lot-worker.md"
            ).read_text(encoding="utf-8")
            self.assertIn("model: synthetic-model\n", worker)
            self.assertIn("effort: high\n", worker)


class ModelLeaderboardTests(unittest.TestCase):
    def make_run(
        self,
        root: Path,
        name: str,
        *,
        model: str,
        retrieval_k: int,
        repeat: int,
        accuracy: float,
        mae: float,
    ) -> Path:
        run = root / name
        restricted = run / "restricted"
        public = run / "public"
        restricted.mkdir(parents=True)
        public.mkdir()
        metadata = {
            "run_id": name,
            "provider": "codex-subagents",
            "model": model,
            "reasoning_effort": "high",
            "temperature": None,
            "retrieval_k": retrieval_k,
            "folds": [0, 1, 2],
            "retrieval_training_folds": [0, 1, 2],
            "prompt_version": "prompt-v1",
            "knowledge_version": "knowledge-v1",
            "input_artifact_sha256": {"input": "hash"},
            "run_purpose": "development",
            "repeat_index": repeat,
        }
        (restricted / "experiment_metadata.json").write_text(
            json.dumps(metadata), encoding="utf-8"
        )
        (restricted / "joined_evaluation.jsonl").write_text(
            json.dumps({"case_key": case_id(1)}) + "\n", encoding="utf-8"
        )
        report = {
            "report_scope": "aggregate_only",
            "metrics": {
                "generated_total": {
                    "exact_accuracy": {"rate": accuracy},
                    "mean_absolute_error": mae,
                    "within_one_accuracy": {"rate": accuracy},
                },
                "selective_ai": {"coverage": {"rate": 1.0}},
                "three_way_consensus": {
                    "agreement": {"rate": 0.5},
                    "reviewer_accuracy": {"rate": 1.0},
                },
            },
        }
        (public / "aggregate_evaluation.json").write_text(
            json.dumps(report), encoding="utf-8"
        )
        return run

    def test_development_leaderboard_groups_repeats_and_selects_accuracy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runs = [
                self.make_run(
                    root,
                    "a1",
                    model="model-a",
                    retrieval_k=3,
                    repeat=1,
                    accuracy=0.6,
                    mae=1.0,
                ),
                self.make_run(
                    root,
                    "a2",
                    model="model-a",
                    retrieval_k=3,
                    repeat=2,
                    accuracy=0.8,
                    mae=0.8,
                ),
                self.make_run(
                    root,
                    "b1",
                    model="model-b",
                    retrieval_k=0,
                    repeat=1,
                    accuracy=0.5,
                    mae=1.2,
                ),
                self.make_run(
                    root,
                    "b2",
                    model="model-b",
                    retrieval_k=0,
                    repeat=2,
                    accuracy=0.5,
                    mae=1.1,
                ),
            ]
            report = build_leaderboard(runs, phase="development", expected_repeats=2)
            selected = report["selected_condition"]
            self.assertEqual(selected["model"], "model-a")
            self.assertEqual(selected["generated_exact_accuracy_mean"], 0.7)
            self.assertEqual(selected["repeat_count"], 2)


if __name__ == "__main__":
    unittest.main()
