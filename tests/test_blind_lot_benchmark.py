from __future__ import annotations

import json
import hashlib
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pydantic import ValidationError


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "py"))

from blind_lot.cache import DiskCache, cache_key, canonical_hash  # noqa: E402
from blind_lot.evaluation import (  # noqa: E402
    derive_prediction_status,
    evaluate_joined,
    normalize_joined_record,
)
from blind_lot.models import BlindLOTResponse, validate_model_response  # noqa: E402
from blind_lot.prompting import build_prompt  # noqa: E402
from blind_lot.providers import (  # noqa: E402
    OpenAIResponsesProvider,
    ProviderConfig,
    RetriableProviderError,
    complete_with_retries,
)
from blind_lot.retrieval import (  # noqa: E402
    LocalRetriever,
    RetrievalCandidate,
    RetrievalHit,
    assert_retrieval_safe,
)
from blind_lot.runner import _public_guard  # noqa: E402
from run_blind_lot_benchmark import format_metric  # noqa: E402


def case_id(number: int) -> str:
    return f"CASE-{number:020x}"


def patient_id(number: int) -> str:
    return f"COTAOLD-{number:020x}"


def blind(trajectory: list[list[str]]) -> dict:
    return {
        "schema_version": "2.0.0",
        "case_key": case_id(999),
        "trajectory": [
            {"order": index, "drugs": drugs}
            for index, drugs in enumerate(trajectory, start=1)
        ],
        "context": {"diagnosis_date": "2020-01-01"},
    }


def response(*, decision: str = "NEW_LINE", abstained: bool = False) -> dict:
    return {
        "schema_version": "1.0.0",
        "ai_lot_count": 2 if decision == "NEW_LINE" else 1,
        "abstained": abstained,
        "transitions": [{
            "transition_index": 1,
            "decision": decision,
            "reason_codes": ["SYNTHETIC_REASON"],
            "evidence_strength": "moderate",
            "retrieved_case_ids": [],
        }],
        "warnings": [],
    }


def candidate(number: int, fold: int, group: str, trajectory: list[list[str]]) -> RetrievalCandidate:
    return RetrievalCandidate(
        case_key=case_id(number),
        patient_key=patient_id(number),
        fold=fold,
        exclusion_group=group,
        trajectory=[{"order": i, "drugs": drugs} for i, drugs in enumerate(trajectory, 1)],
        reviewer_consensus_lot=2,
    )


class ResponseSchemaTests(unittest.TestCase):
    def test_valid_response(self) -> None:
        result = validate_model_response(response(), 2)
        self.assertEqual(result.ai_lot_count, 2)

    def test_incorrect_transition_count_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            validate_model_response(response(), 3)

    def test_inconsistent_lot_total_is_rejected(self) -> None:
        value = response()
        value["ai_lot_count"] = 4
        with self.assertRaises(ValidationError):
            validate_model_response(value, 2)

    def test_insufficient_information_requires_abstention(self) -> None:
        with self.assertRaises(ValidationError):
            validate_model_response(response(decision="INSUFFICIENT_INFORMATION"), 2)

    def test_empty_trajectory_behavior(self) -> None:
        value = {
            "schema_version": "1.0.0", "ai_lot_count": 0, "abstained": True,
            "transitions": [], "warnings": ["NO_EVENTS"],
        }
        self.assertEqual(validate_model_response(value, 0).ai_lot_count, 0)
        value["abstained"] = False
        with self.assertRaises(ValueError):
            validate_model_response(value, 0)


class RetrievalTests(unittest.TestCase):
    def setUp(self) -> None:
        self.items = [
            candidate(1, 0, "EXCL-A", [["a"], ["b"]]),
            candidate(2, 1, "EXCL-B", [["a"], ["b"]]),
            candidate(3, 2, "EXCL-C", [["a"], ["c"]]),
            candidate(4, 3, "EXCL-D", [["x"], ["y"]]),
            candidate(5, 4, "EXCL-E", [["a"], ["b"], ["c"]]),
            candidate(6, 1, "EXCL-TARGET", [["a"], ["b"]]),
            candidate(7, 2, "EXCL-F", [["q"]]),
        ]
        self.retriever = LocalRetriever(self.items)
        self.kwargs = {
            "target_patient_key": patient_id(1),
            "target_fold": 0,
            "target_exclusion_group": "EXCL-TARGET",
        }
        self.trajectory = [{"order": 1, "drugs": ["a"]}, {"order": 2, "drugs": ["b"]}]

    def test_deterministic_retrieval(self) -> None:
        first = self.retriever.retrieve(self.trajectory, k=5, **self.kwargs)
        second = self.retriever.retrieve(self.trajectory, k=5, **self.kwargs)
        self.assertEqual([hit.debug_dict() for hit in first], [hit.debug_dict() for hit in second])

    def test_fold_and_exclusion_group_are_excluded(self) -> None:
        hits = self.retriever.retrieve(self.trajectory, k=5, **self.kwargs)
        self.assertTrue(all(hit.candidate.fold != 0 for hit in hits))
        self.assertTrue(all(hit.candidate.exclusion_group != "EXCL-TARGET" for hit in hits))
        assert_retrieval_safe(hits, **self.kwargs)

    def test_supported_k_values(self) -> None:
        for k in (0, 3, 5):
            self.assertEqual(len(self.retriever.retrieve(self.trajectory, k=k, **self.kwargs)), k)

    def test_diagnostics_do_not_change_retrieval_or_rendered_model_input(self) -> None:
        hits = self.retriever.retrieve(self.trajectory, k=3, **self.kwargs)
        bundle = build_prompt(blind([["a"], ["b"]]), hits)
        original_ids = [hit.candidate.case_key for hit in hits]
        original_input = bundle.patient_prompt
        expected_context = json.dumps(
            json.loads(original_input)["permitted_retrieved_training_examples"],
            sort_keys=True,
            separators=(",", ":"),
        )
        diagnostics = [
            {
                "rank": rank,
                "case_key": hit.candidate.case_key,
                "rendered": rendered,
                "sha256": hashlib.sha256(rendered.encode()).hexdigest(),
            }
            for rank, (hit, rendered) in enumerate(
                zip(hits, bundle.rendered_retrieval_demonstrations), start=1
            )
        ]
        self.assertEqual([item["case_key"] for item in diagnostics], original_ids)
        self.assertEqual([item["rank"] for item in diagnostics], [1, 2, 3])
        self.assertEqual(bundle.patient_prompt, original_input)
        self.assertEqual(bundle.rendered_retrieval_context, expected_context)
        self.assertEqual(
            hashlib.sha256(bundle.rendered_retrieval_context.encode()).hexdigest(),
            hashlib.sha256(expected_context.encode()).hexdigest(),
        )


class PromptLeakageTests(unittest.TestCase):
    def test_forbidden_target_fields_are_rejected(self) -> None:
        target = blind([["a"], ["b"]])
        target["cota_lot"] = 987
        with self.assertRaises(ValidationError):
            build_prompt(target, [])

    def test_target_answers_and_metadata_absent_from_prompt(self) -> None:
        hit = RetrievalHit(candidate(2, 1, "EXCL-B", [["a"], ["b"]]), 1.0, {})
        bundle = build_prompt(blind([["a"], ["b"]]), [hit])
        value = json.loads(bundle.patient_prompt)
        target = value["target"]
        serialized_target = json.dumps(target)
        for name in (
            "cota_lot", "reviewer_lot", "source_line_number", "source_group_index",
            "source_line_start_date", "source_line_end_date", "algorithm_lot",
            "algorithm_flags", "exclusion_group", "patient_key",
        ):
            self.assertNotIn(name, serialized_target)
        self.assertEqual(
            value["permitted_retrieved_training_examples"][0]["reviewer_consensus_patient_total_lot"],
            2,
        )


class ProviderTests(unittest.TestCase):
    def test_openai_request_separates_stable_prefix_and_disables_storage(self) -> None:
        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return None

            def read(self) -> bytes:
                return json.dumps({
                    "output": [{"content": [{"type": "output_text", "text": json.dumps(response())}]}]
                }).encode()

        captured = {}

        def fake_urlopen(request, timeout):
            captured["body"] = json.loads(request.data)
            captured["timeout"] = timeout
            return FakeResponse()

        with patch.dict(os.environ, {"OPENAI_API_KEY": "synthetic-key"}), patch(
            "blind_lot.providers.urllib.request.urlopen", side_effect=fake_urlopen
        ):
            provider = OpenAIResponsesProvider(ProviderConfig(model="synthetic-model", request_timeout=7))
            provider.complete("stable-prefix", "patient-only")
        self.assertEqual(captured["body"]["instructions"], "stable-prefix")
        self.assertEqual(captured["body"]["input"], "patient-only")
        self.assertFalse(captured["body"]["store"])
        self.assertEqual(captured["body"]["text"]["format"]["type"], "json_schema")
        self.assertEqual(captured["timeout"], 7)

    def test_malformed_output_is_retried(self) -> None:
        class Provider:
            calls = 0

            def complete(self, stable_prefix: str, patient_prompt: str) -> str:
                self.calls += 1
                return "not-json" if self.calls == 1 else json.dumps(response())

        provider = Provider()
        with patch("blind_lot.providers.time.sleep"):
            _, attempts = complete_with_retries(provider, "stable", "patient", regimen_event_count=2, retry_count=1)
        self.assertEqual(attempts, 2)

    def test_transport_failure_is_retried(self) -> None:
        class Provider:
            calls = 0

            def complete(self, stable_prefix: str, patient_prompt: str) -> str:
                self.calls += 1
                if self.calls == 1:
                    raise RetriableProviderError("synthetic transport failure")
                return json.dumps(response())

        provider = Provider()
        with patch("blind_lot.providers.time.sleep"):
            _, attempts = complete_with_retries(provider, "stable", "patient", regimen_event_count=2, retry_count=1)
        self.assertEqual(attempts, 2)


class CacheTests(unittest.TestCase):
    def test_cache_key_stability_and_current_validation(self) -> None:
        kwargs = {
            "model": "model", "reasoning_effort": "low", "temperature": None,
            "prompt_version": "p1", "knowledge_version": "k1", "retrieval_k": 3,
            "retrieved_example_ids": [case_id(2)],
            "blind_patient_input_hash": canonical_hash(blind([["a"], ["b"]])),
        }
        self.assertEqual(cache_key(**kwargs), cache_key(**kwargs))
        with tempfile.TemporaryDirectory() as directory:
            cache = DiskCache(Path(directory))
            model = BlindLOTResponse.model_validate(response())
            key = cache_key(**kwargs)
            cache.put(key, model)
            self.assertEqual(cache.get(key, 2), model)


class AggregateEvaluationTests(unittest.TestCase):
    @staticmethod
    def row(
        reviewer: int, algorithm: int, cota: int, ai: object, abstained: object,
        number: int = 1,
    ) -> dict:
        return {
            "case_key": case_id(number),
            "reviewer_lot": reviewer,
            "ai_lot": ai,
            "ai_abstained": abstained,
            "algorithm_lot": algorithm,
            "cota_lot": cota,
        }

    def test_abstained_but_numerically_correct_routes_to_review(self) -> None:
        normalized = normalize_joined_record(self.row(5, 5, 5, 5, True))
        self.assertTrue(normalized["has_generated_ai_total"])
        self.assertTrue(normalized["generated_total_correct"])
        self.assertFalse(normalized["usable_ai_vote"])
        self.assertIsNone(normalized["non_abstained_prediction_correct"])
        self.assertFalse(normalized["ai_algorithm_agreement"])
        self.assertFalse(normalized["ai_cota_agreement"])
        self.assertFalse(normalized["three_way_agreement"])
        self.assertEqual(normalized["routing_decision"], "human_review")
        self.assertEqual(normalized["routing_reason"], "ai_abstained")
        metrics = evaluate_joined([normalized], bootstrap_replicates=10)
        self.assertEqual(metrics["generated_total"]["exact_accuracy"]["count"], 1)
        self.assertEqual(metrics["selective_ai"]["exact_accuracy"]["denominator"], 0)

    def test_correct_usable_three_way_vote_routes_to_consensus(self) -> None:
        normalized = normalize_joined_record(self.row(4, 4, 4, 4, False))
        self.assertTrue(normalized["usable_ai_vote"])
        self.assertTrue(normalized["non_abstained_prediction_correct"])
        self.assertTrue(normalized["three_way_agreement"])
        self.assertEqual(normalized["routing_decision"], "consensus_candidate")
        self.assertEqual(normalized["routing_reason"], "three_way_consensus")

    def test_incorrect_two_way_candidate_rejection_reasons(self) -> None:
        rows = [
            self.row(4, 5, 5, 4, False, 1),
            self.row(4, 5, 5, 5, True, 2),
        ]
        incremental = evaluate_joined(rows, bootstrap_replicates=10)["incremental_effect_of_ai"]
        self.assertEqual(incremental["silent_failures_rejected_by_disagreement"], 1)
        self.assertEqual(incremental["silent_failures_rejected_by_abstention"], 1)
        self.assertEqual(incremental["two_way_silent_failures_rejected_by_ai"], 2)

    def test_abstention_rejects_correct_two_way_candidate(self) -> None:
        incremental = evaluate_joined(
            [self.row(5, 5, 5, 5, True)], bootstrap_replicates=10
        )["incremental_effect_of_ai"]
        self.assertEqual(incremental["correct_two_way_candidates_rejected_by_abstention"], 1)
        self.assertEqual(incremental["previously_correct_two_way_agreements_rejected_by_ai"], 1)

    def test_zero_usable_votes_has_null_selective_accuracy_and_ci(self) -> None:
        rows = [self.row(5, 5, 5, 5, True, number) for number in range(1, 4)]
        selective = evaluate_joined(rows, bootstrap_replicates=20)["selective_ai"]
        self.assertEqual(selective["coverage"], {
            "count": 0,
            "denominator": 3,
            "rate": 0.0,
            "bootstrap_95_ci": {"lower": 0.0, "upper": 0.0, "replicates_used": 20},
        })
        self.assertEqual(selective["exact_accuracy"]["count"], 0)
        self.assertEqual(selective["exact_accuracy"]["denominator"], 0)
        self.assertIsNone(selective["exact_accuracy"]["rate"])
        self.assertEqual(selective["exact_accuracy"]["bootstrap_95_ci"], {
            "lower": None, "upper": None, "replicates_used": 0,
        })

    def test_invalid_ai_totals_and_boolean_are_not_generated_votes(self) -> None:
        invalid_values = [None, True, False, -1, 1.5, "5"]
        for value in invalid_values:
            with self.subTest(value=value):
                status = derive_prediction_status(
                    ai_lot_count=value,
                    abstained=False,
                    reviewer_lot=5,
                    algorithm_lot_count=5,
                    cota_lot_count=5,
                )
                self.assertFalse(status.has_generated_ai_total)
                self.assertFalse(status.usable_ai_vote)
                self.assertTrue(status.invalid_ai_output)
                self.assertEqual(status.routing_decision, "human_review")
                self.assertEqual(status.routing_reason, "invalid_ai_output")

    def test_saved_five_patient_k3_and_k5_replay(self) -> None:
        expectations = {
            "20260713T231652744088Z-b539b06cfd": (3, 2, 1, 1),
            "20260713T231836398449Z-890c062947": (3, 1, 0, 0),
        }
        for run_id, expected in expectations.items():
            path = ROOT / "artifacts" / "benchmarks" / "runs" / run_id / "restricted" / "joined_evaluation.jsonl"
            rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
            metrics = evaluate_joined(rows, bootstrap_replicates=20)
            observed = (
                metrics["generated_total"]["exact_accuracy"]["count"],
                metrics["selective_ai"]["coverage"]["count"],
                metrics["selective_ai"]["exact_accuracy"]["count"],
                metrics["three_way_consensus"]["agreement"]["count"],
            )
            self.assertEqual(observed, expected)
            self.assertEqual(metrics["eligible_patients"], 5)

    def test_cli_metric_format_handles_null_rate(self) -> None:
        self.assertEqual(
            format_metric("Non-abstained exact accuracy", {"count": 0, "denominator": 0, "rate": None}),
            "Non-abstained exact accuracy: N/A (0/0)",
        )

    def test_two_way_three_way_and_silent_failure_accounting(self) -> None:
        rows = [
            {"case_key": case_id(1), "reviewer_lot": 2, "ai_lot": 2, "ai_abstained": False, "algorithm_lot": 2, "cota_lot": 2},
            {"case_key": case_id(2), "reviewer_lot": 2, "ai_lot": 2, "ai_abstained": False, "algorithm_lot": 3, "cota_lot": 3},
            {"case_key": case_id(3), "reviewer_lot": 4, "ai_lot": 5, "ai_abstained": False, "algorithm_lot": 4, "cota_lot": 4},
            {"case_key": case_id(4), "reviewer_lot": 1, "ai_lot": 2, "ai_abstained": False, "algorithm_lot": 2, "cota_lot": 2},
        ]
        metrics = evaluate_joined(rows, bootstrap_replicates=20)
        incremental = metrics["incremental_effect_of_ai"]
        self.assertEqual(incremental["two_way_silent_failures_rejected_by_ai"], 1)
        self.assertEqual(incremental["previously_correct_two_way_agreements_rejected_by_ai"], 1)
        self.assertEqual(incremental["two_way_agreement_cases_retained"], 2)
        self.assertEqual(incremental["net_reduction_in_incorrect_auto_accept_candidates"], 1)
        self.assertEqual(incremental["three_way_silent_failures_remaining"], 1)
        self.assertEqual(incremental["coverage_lost_per_silent_failure_prevented"], 1.0)
        self.assertEqual(incremental["all_candidates_rejected_per_silent_failure_prevented"], 2.0)

    def test_public_report_guard_rejects_patient_ids(self) -> None:
        _public_guard({"report_scope": "aggregate_only", "count": 3})
        with self.assertRaises(ValueError):
            _public_guard({"report_scope": "aggregate_only", "note": case_id(1)})

    def test_policy_all_three_agree_correctly(self) -> None:
        row = normalize_joined_record(self.row(5, 5, 5, 5, False))
        self.assertEqual(row["vote_pattern"], "all_three_agree")
        self.assertTrue(all(
            row[name] for name in (
                "accepted_by_algorithm_only", "accepted_by_cota_only",
                "accepted_by_usable_ai_only", "accepted_by_algorithm_cota",
                "accepted_by_algorithm_ai", "accepted_by_cota_ai",
                "accepted_by_three_way",
            )
        ))
        self.assertTrue(row["three_way_policy_correct"])

    def test_policy_all_three_agree_incorrectly_is_false_accept(self) -> None:
        metrics = evaluate_joined([self.row(4, 5, 5, 5, False)], bootstrap_replicates=10)
        policy = metrics["routing_policy_analysis"]["three_way_agreement"]
        self.assertEqual(policy["accepted_count"], 1)
        self.assertEqual(policy["incorrect_accepted_count"], 1)
        self.assertEqual(policy["false_accept_rate_among_accepted"]["rate"], 1.0)

    def test_policy_delta_has_unambiguous_overlap_and_net_fields(self) -> None:
        metrics = evaluate_joined([
            self.row(2, 2, 2, 2, False, 1),
            self.row(2, 2, 3, 2, False, 2),
            self.row(2, 2, 2, 3, False, 3),
        ], bootstrap_replicates=5)
        delta = metrics["routing_policy_analysis"]["algorithm_ai_agreement"][
            "comparisons"
        ]["versus_three_way_agreement"]
        self.assertEqual(delta["accepted_by_both_count"], 1)
        self.assertEqual(delta["accepted_only_by_policy_count"], 1)
        self.assertEqual(delta["accepted_only_by_baseline_count"], 0)
        self.assertEqual(delta["net_accepted_count_difference"], 1)
        self.assertTrue(delta["additional_accepted_patients_deprecated"])

    def test_algorithm_ai_agree_cota_dissent_correct(self) -> None:
        metrics = evaluate_joined([self.row(4, 5, 4, 5, False)], bootstrap_replicates=10)
        policies = metrics["routing_policy_analysis"]
        self.assertEqual(policies["algorithm_ai_agreement"]["incorrect_accepted_count"], 1)
        self.assertEqual(policies["cota_only"]["correct_accepted_count"], 1)
        self.assertEqual(policies["three_way_agreement"]["accepted_count"], 0)
        pattern = metrics["vote_pattern_analysis"]["patterns"]["algorithm_ai_agree_cota_differs"]
        self.assertEqual(pattern["majority_vs_dissent"]["dissenting_vote_correct"], 1)

    def test_algorithm_ai_agree_majority_correct(self) -> None:
        metrics = evaluate_joined([self.row(5, 5, 4, 5, False)], bootstrap_replicates=10)
        policies = metrics["routing_policy_analysis"]
        self.assertEqual(policies["algorithm_ai_agreement"]["correct_accepted_count"], 1)
        self.assertEqual(policies["cota_only"]["incorrect_accepted_count"], 1)
        pattern = metrics["vote_pattern_analysis"]["patterns"]["algorithm_ai_agree_cota_differs"]
        self.assertEqual(pattern["majority_vs_dissent"]["majority_correct"], 1)

    def test_algorithm_cota_agree_ai_abstains(self) -> None:
        row = normalize_joined_record(self.row(5, 5, 5, 5, True))
        self.assertTrue(row["accepted_by_algorithm_cota"])
        self.assertFalse(row["accepted_by_algorithm_ai"])
        self.assertFalse(row["accepted_by_cota_ai"])
        self.assertFalse(row["accepted_by_three_way"])
        self.assertEqual(row["vote_pattern"], "algorithm_cota_agree_ai_abstains")

    def test_algorithm_cota_false_accept_ai_dissent_correct(self) -> None:
        metrics = evaluate_joined([self.row(4, 5, 5, 4, False)], bootstrap_replicates=10)
        self.assertEqual(
            metrics["routing_policy_analysis"]["algorithm_cota_agreement"]["incorrect_accepted_count"], 1
        )
        pattern = metrics["vote_pattern_analysis"]["patterns"]["algorithm_cota_agree_ai_differs"]
        self.assertEqual(pattern["winner_attribution"]["ai_uniquely_correct"], 1)
        self.assertEqual(pattern["majority_vs_dissent"]["dissenting_vote_correct"], 1)

    def test_cota_ai_agree_correctly(self) -> None:
        metrics = evaluate_joined([self.row(4, 5, 4, 4, False)], bootstrap_replicates=10)
        policies = metrics["routing_policy_analysis"]
        self.assertEqual(policies["cota_ai_agreement"]["correct_accepted_count"], 1)
        self.assertEqual(policies["algorithm_ai_agreement"]["accepted_count"], 0)
        pattern = metrics["vote_pattern_analysis"]["patterns"]["cota_ai_agree_algorithm_differs"]
        self.assertEqual(pattern["majority_vs_dissent"]["majority_correct"], 1)

    def test_all_three_differ_ai_uniquely_correct(self) -> None:
        metrics = evaluate_joined([self.row(4, 5, 6, 4, False)], bootstrap_replicates=10)
        policies = metrics["routing_policy_analysis"]
        self.assertEqual(policies["algorithm_cota_agreement"]["accepted_count"], 0)
        self.assertEqual(policies["algorithm_ai_agreement"]["accepted_count"], 0)
        self.assertEqual(policies["cota_ai_agreement"]["accepted_count"], 0)
        self.assertEqual(policies["usable_ai_only"]["correct_accepted_count"], 1)
        pattern = metrics["vote_pattern_analysis"]["patterns"]["all_three_differ"]
        self.assertEqual(pattern["winner_attribution"]["ai_uniquely_correct"], 1)

    def test_no_policy_accepts_and_conditional_rates_are_null(self) -> None:
        metrics = evaluate_joined(
            [self.row(4, -1, -1, None, False)], bootstrap_replicates=10
        )
        for policy in metrics["routing_policy_analysis"].values():
            if not isinstance(policy, dict) or "accepted_count" not in policy:
                continue
            self.assertEqual(policy["accepted_count"], 0)
            self.assertIsNone(policy["reviewer_accuracy_among_accepted"]["rate"])
            self.assertIsNone(policy["false_accept_rate_among_accepted"]["rate"])
            self.assertEqual(
                policy["reviewer_accuracy_among_accepted"]["bootstrap_95_ci"]["replicates_used"], 0
            )

    def test_full_k5_policy_regression_and_exhaustive_patterns(self) -> None:
        if os.environ.get("RUN_RESTRICTED_REGRESSION_TESTS") != "1":
            self.skipTest("set RUN_RESTRICTED_REGRESSION_TESTS=1 for local frozen artifacts")
        path = ROOT / "artifacts/benchmarks/runs/20260713T235819111155Z-22d38bbc4e/restricted/joined_evaluation.jsonl"
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
        metrics = evaluate_joined(rows, bootstrap_replicates=20)
        policies = metrics["routing_policy_analysis"]
        self.assertEqual((policies["algorithm_cota_agreement"]["accepted_count"], policies["algorithm_cota_agreement"]["correct_accepted_count"]), (17, 13))
        self.assertEqual((policies["algorithm_ai_agreement"]["accepted_count"], policies["algorithm_ai_agreement"]["correct_accepted_count"]), (5, 3))
        self.assertEqual((policies["cota_ai_agreement"]["accepted_count"], policies["cota_ai_agreement"]["correct_accepted_count"]), (3, 2))
        self.assertEqual((policies["three_way_agreement"]["accepted_count"], policies["three_way_agreement"]["correct_accepted_count"]), (2, 2))
        pattern = metrics["vote_pattern_analysis"]
        self.assertTrue(pattern["mutually_exclusive_and_exhaustive"])
        self.assertEqual(pattern["patient_count_sum"], 27)
        stratum = pattern["patterns"]["algorithm_ai_agree_cota_differs"]
        self.assertEqual(stratum["patient_count"], 3)
        self.assertEqual(stratum["majority_vs_dissent"], {
            "majority_correct": 1, "dissenting_vote_correct": 1, "none_correct": 1,
        })


if __name__ == "__main__":
    unittest.main()
