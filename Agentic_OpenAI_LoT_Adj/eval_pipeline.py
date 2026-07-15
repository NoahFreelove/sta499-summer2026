from __future__ import annotations

import asyncio
import json
import random
import traceback
from pathlib import Path
from typing import Any

import pandas as pd

from config import (
    MAX_CONCURRENT_PATIENTS,
    NUMBER_OF_REFERENCE_EXAMPLES,
    OLD_ADJUDICATED_FILE,
    OLD_COTA_SHEET,
    PATIENT_ID_COLUMN,
    REFERENCE_LABEL_COLUMN,
    REFERENCE_RESIDUAL_COLUMN,
)
from decision import compare_results, make_decision
from eval_metrics import (
    aggregate_metrics,
    mapping_from_result,
    score_mapping,
    transition_flags,
)
from pipeline import run_critique, run_independent_agents, safe_name, write_json
from preprocessing import preprocess_new_patient
from reference_retrieval import (
    ReferenceLibrary,
    build_reference_library_from_dataframe,
    build_patient_reference,
    collapse_historical_rows,
    normalize_columns,
    retrieve_similar_examples_excluding,
)

EVAL_FRACTION = 0.15
RANDOM_SEED = 42
EVAL_OUTPUT_DIR = Path(__file__).resolve().parent / "eval_output"


def load_collapsed_labeled_data() -> pd.DataFrame:
    if not OLD_ADJUDICATED_FILE.exists():
        raise FileNotFoundError(f"Old adjudicated workbook not found: {OLD_ADJUDICATED_FILE}")

    raw = normalize_columns(pd.read_excel(OLD_ADJUDICATED_FILE, sheet_name=OLD_COTA_SHEET))
    collapsed = collapse_historical_rows(raw)
    collapsed[REFERENCE_LABEL_COLUMN] = pd.to_numeric(
        collapsed[REFERENCE_LABEL_COLUMN], errors="coerce"
    )
    collapsed = collapsed.dropna(subset=[PATIENT_ID_COLUMN, REFERENCE_LABEL_COLUMN])
    collapsed[PATIENT_ID_COLUMN] = collapsed[PATIENT_ID_COLUMN].astype(str).str.strip()
    collapsed[REFERENCE_LABEL_COLUMN] = collapsed[REFERENCE_LABEL_COLUMN].astype(int)
    return collapsed


def split_patient_ids(patient_ids: list[str]) -> tuple[list[str], list[str]]:
    ids = sorted(set(patient_ids))
    if len(ids) < 2:
        raise ValueError("At least two labeled patients are required for evaluation.")

    rng = random.Random(RANDOM_SEED)
    rng.shuffle(ids)
    eval_count = max(1, round(len(ids) * EVAL_FRACTION))
    eval_count = min(eval_count, len(ids) - 1)
    eval_ids = sorted(ids[:eval_count])
    train_ids = sorted(ids[eval_count:])
    return train_ids, eval_ids


def prepare_eval_case(patient_df: pd.DataFrame) -> tuple[dict[str, Any], dict[int, int]]:
    patient_id = str(patient_df[PATIENT_ID_COLUMN].iloc[0]).strip()
    reference = build_patient_reference(patient_id, patient_df)
    truth_mapping = {
        int(line["vendor_lot"]): int(line["adjudicated_lot"])
        for line in reference["lines"]
    }

    # The original COTA decision is represented by its own vendor LoT
    # sequence. Compare transition decisions, not absolute numbering, so a
    # downstream carryover shift is not counted as another mistake.
    vendor_mapping = {
        int(line["vendor_lot"]): int(line["vendor_lot"])
        for line in reference["lines"]
    }
    vendor_decisions = transition_flags(vendor_mapping)
    truth_decisions = transition_flags(truth_mapping)
    transition_lots = sorted(truth_decisions)[1:]
    cota_error_lots = [
        lot
        for lot in transition_lots
        if vendor_decisions.get(lot) != truth_decisions.get(lot)
    ]

    hidden = patient_df.drop(
        columns=[REFERENCE_LABEL_COLUMN, REFERENCE_RESIDUAL_COLUMN],
        errors="ignore",
    ).copy()
    patient = preprocess_new_patient(hidden)
    patient["_eval_cota_patient_misclassified"] = bool(cota_error_lots)
    patient["_eval_cota_misclassified_transition_count"] = len(cota_error_lots)
    patient["_eval_cota_misclassified_vendor_lots"] = cota_error_lots
    return patient, truth_mapping


async def evaluate_patient(
    patient_df: pd.DataFrame,
    reference_library: ReferenceLibrary,
    all_eval_ids: set[str],
    semaphore: asyncio.Semaphore,
) -> dict[str, Any]:
    async with semaphore:
        patient_id = str(patient_df[PATIENT_ID_COLUMN].iloc[0]).strip()
        patient_dir = EVAL_OUTPUT_DIR / "patients" / safe_name(patient_id)

        try:
            patient, truth = prepare_eval_case(patient_df)
            cota_patient_misclassified = bool(
                patient.pop("_eval_cota_patient_misclassified")
            )
            cota_misclassified_transition_count = int(
                patient.pop("_eval_cota_misclassified_transition_count")
            )
            cota_misclassified_vendor_lots = list(
                patient.pop("_eval_cota_misclassified_vendor_lots")
            )

            references = retrieve_similar_examples_excluding(
                patient=patient,
                reference_library=reference_library,
                top_k=NUMBER_OF_REFERENCE_EXAMPLES,
                excluded_patient_ids=all_eval_ids,
            )

            leaked_ids = {
                str(item["reference_patient_id"]).strip()
                for item in references
            } & all_eval_ids
            if leaked_ids:
                raise RuntimeError(f"Evaluation leakage detected: {sorted(leaked_ids)}")

            worker, auditor = await run_independent_agents(patient, references)
            comparison = compare_results(worker, auditor)
            critique = None
            if not comparison["exact_mapping_agreement"]:
                critique = await run_critique(patient, references, worker, auditor, comparison)

            decision = make_decision(worker, auditor, comparison, critique)
            selected = None
            if decision.selected_result == "worker":
                selected = worker
            elif decision.selected_result == "auditor":
                selected = auditor

            worker_score = score_mapping(mapping_from_result(worker), truth)
            auditor_score = score_mapping(mapping_from_result(auditor), truth)
            selected_score = score_mapping(
                mapping_from_result(selected) if selected is not None else None,
                truth,
            )

            write_json(patient_dir / "input_labels_hidden.json", patient)
            write_json(patient_dir / "ground_truth.json", truth)
            write_json(patient_dir / "reference_examples.json", references)
            write_json(patient_dir / "worker.json", worker.model_dump(mode="json"))
            write_json(patient_dir / "auditor.json", auditor.model_dump(mode="json"))
            write_json(patient_dir / "comparison.json", comparison)
            if critique is not None:
                write_json(patient_dir / "critique.json", critique.model_dump(mode="json"))
            write_json(patient_dir / "decision.json", decision.model_dump(mode="json"))

            return {
                "patient_id": patient_id,
                "status": decision.status,
                "selected_result": decision.selected_result,
                "cota_patient_misclassified": cota_patient_misclassified,
                "cota_misclassified_transition_count": (
                    cota_misclassified_transition_count
                ),
                "cota_misclassified_vendor_lots": json.dumps(
                    cota_misclassified_vendor_lots
                ),
                "sent_to_human_review": decision.status == "human_review",
                "triage_outcome": (
                    "caught_error_for_review"
                    if cota_patient_misclassified
                    and decision.status == "human_review"
                    else "missed_error_autoaccepted"
                    if cota_patient_misclassified
                    else "unnecessary_review"
                    if decision.status == "human_review"
                    else "safe_autoaccept"
                ),
                "truth_line_count": len(truth),
                "truth_transition_count": worker_score["total_transitions"],
                "worker_exact_match": worker_score["exact_match"],
                "auditor_exact_match": auditor_score["exact_match"],
                "selected_exact_match": selected_score["exact_match"],
                "worker_correct_lines": worker_score["correct_lines"],
                "auditor_correct_lines": auditor_score["correct_lines"],
                "selected_correct_lines": selected_score["correct_lines"],
                "worker_correct_transitions": worker_score["correct_transitions"],
                "auditor_correct_transitions": auditor_score["correct_transitions"],
                "selected_correct_transitions": selected_score["correct_transitions"],
                "worker_auditor_agree": comparison["exact_mapping_agreement"],
                "selected_merge_tp": selected_score["merge_tp"],
                "selected_merge_fp": selected_score["merge_fp"],
                "selected_merge_fn": selected_score["merge_fn"],
                "references_used": len(references),
                "reference_leakage": False,
                "error": "",
            }
        except Exception as error:
            write_json(
                patient_dir / "error.json",
                {
                    "patient_id": patient_id,
                    "error_type": type(error).__name__,
                    "error": str(error),
                    "traceback": traceback.format_exc(),
                },
            )
            return {
                "patient_id": patient_id,
                "status": "processing_error",
                "selected_result": "none",
                "cota_patient_misclassified": False,
                "cota_misclassified_transition_count": 0,
                "cota_misclassified_vendor_lots": "[]",
                "sent_to_human_review": False,
                "triage_outcome": "processing_error",
                "truth_line_count": 0,
                "truth_transition_count": 0,
                "worker_exact_match": False,
                "auditor_exact_match": False,
                "selected_exact_match": False,
                "worker_correct_lines": 0,
                "auditor_correct_lines": 0,
                "selected_correct_lines": 0,
                "worker_correct_transitions": 0,
                "auditor_correct_transitions": 0,
                "selected_correct_transitions": 0,
                "worker_auditor_agree": False,
                "selected_merge_tp": 0,
                "selected_merge_fp": 0,
                "selected_merge_fn": 0,
                "references_used": 0,
                "reference_leakage": False,
                "error": str(error),
            }


async def main() -> None:
    EVAL_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    labeled = load_collapsed_labeled_data()
    train_ids, eval_ids = split_patient_ids(labeled[PATIENT_ID_COLUMN].tolist())

    train_df = labeled[labeled[PATIENT_ID_COLUMN].isin(train_ids)].copy()
    eval_df = labeled[labeled[PATIENT_ID_COLUMN].isin(eval_ids)].copy()
    reference_library = build_reference_library_from_dataframe(train_df)

    split_record = {
        "random_seed": RANDOM_SEED,
        "eval_fraction": EVAL_FRACTION,
        "training_patient_count": len(train_ids),
        "evaluation_patient_count": len(eval_ids),
        "training_patient_ids": train_ids,
        "evaluation_patient_ids": eval_ids,
    }
    write_json(EVAL_OUTPUT_DIR / "split.json", split_record)

    semaphore = asyncio.Semaphore(MAX_CONCURRENT_PATIENTS)
    eval_id_set = set(eval_ids)
    tasks = [
        evaluate_patient(patient_df.copy(), reference_library, eval_id_set, semaphore)
        for _, patient_df in eval_df.groupby(PATIENT_ID_COLUMN, sort=False)
    ]
    rows = await asyncio.gather(*tasks)

    detail_df = pd.DataFrame(rows)
    detail_df.to_csv(EVAL_OUTPUT_DIR / "evaluation_details.csv", index=False)
    metrics = aggregate_metrics(rows)
    write_json(EVAL_OUTPUT_DIR / "evaluation_metrics.json", metrics)

    print(f"Training patients: {len(train_ids)}")
    print(f"Evaluation patients: {len(eval_ids)}")
    print(f"Details: {EVAL_OUTPUT_DIR / 'evaluation_details.csv'}")
    print(f"Metrics: {EVAL_OUTPUT_DIR / 'evaluation_metrics.json'}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    asyncio.run(main())
