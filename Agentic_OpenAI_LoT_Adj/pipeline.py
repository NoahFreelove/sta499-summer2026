from __future__ import annotations

import asyncio
import json
import re
import traceback
from pathlib import Path
from typing import Any

import pandas as pd
from agents import Runner

from agent_definitions import auditor_agent, critique_agent, worker_agent
from config import (
    MAX_CONCURRENT_PATIENTS,
    NEW_COTA_FILE,
    NEW_COTA_SHEET,
    NUMBER_OF_REFERENCE_EXAMPLES,
    OLD_ADJUDICATED_FILE,
    OLD_COTA_SHEET,
    OUTPUT_DIR,
    PATIENT_ID_COLUMN,
)
from decision import compare_results, make_decision
from preprocessing import preprocess_new_patient
from reference_retrieval import ReferenceLibrary, build_reference_library, retrieve_similar_examples
from schemas import AdjudicationResult, AuditCritique


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, default=str)


def safe_name(value: str) -> str:
    return re.sub(r'[<>:"/\\|?*]', "_", value.strip()) or "unknown_patient"


def load_new_data() -> pd.DataFrame:
    if not NEW_COTA_FILE.exists():
        raise FileNotFoundError(f"New COTA workbook not found: {NEW_COTA_FILE}")
    df = pd.read_excel(NEW_COTA_FILE, sheet_name=NEW_COTA_SHEET)
    df.columns = [str(column).strip() for column in df.columns]
    if PATIENT_ID_COLUMN not in df.columns:
        raise ValueError(f"Missing {PATIENT_ID_COLUMN}. Available columns: {list(df.columns)}")
    df["_pipeline_patient_id"] = (
        df[PATIENT_ID_COLUMN].replace(r"^\s*$", pd.NA, regex=True).ffill()
    )
    return df.dropna(subset=["_pipeline_patient_id"])


def payload(patient: dict[str, Any], references: list[dict[str, Any]]) -> str:
    return json.dumps(
        {
            "task": "Correct vendor COTA LoT oversplitting for this patient.",
            "new_patient": patient,
            "historical_labeled_examples": references,
        },
        indent=2,
        ensure_ascii=False,
        default=str,
    )


async def run_independent_agents(
    patient: dict[str, Any], references: list[dict[str, Any]]
) -> tuple[AdjudicationResult, AdjudicationResult]:
    agent_input = payload(patient, references)
    worker_run, auditor_run = await asyncio.gather(
        Runner.run(worker_agent, agent_input, max_turns=3),
        Runner.run(auditor_agent, agent_input, max_turns=3),
    )
    worker, auditor = worker_run.final_output, auditor_run.final_output
    if not isinstance(worker, AdjudicationResult) or not isinstance(auditor, AdjudicationResult):
        raise TypeError("Worker or auditor returned the wrong structured output type.")
    return worker, auditor


async def run_critique(
    patient: dict[str, Any],
    references: list[dict[str, Any]],
    worker: AdjudicationResult,
    auditor: AdjudicationResult,
    comparison: dict[str, Any],
) -> AuditCritique:
    critique_input = json.dumps(
        {
            "new_patient": patient,
            "historical_labeled_examples": references,
            "worker_result": worker.model_dump(mode="json"),
            "auditor_result": auditor.model_dump(mode="json"),
            "deterministic_comparison": comparison,
        },
        indent=2,
        ensure_ascii=False,
        default=str,
    )
    result = await Runner.run(critique_agent, critique_input, max_turns=3)
    if not isinstance(result.final_output, AuditCritique):
        raise TypeError("Critique agent returned the wrong structured output type.")
    return result.final_output


async def process_patient(
    raw_patient_df: pd.DataFrame,
    reference_library: ReferenceLibrary,
    semaphore: asyncio.Semaphore,
) -> dict[str, Any]:
    async with semaphore:
        patient_id = str(raw_patient_df["_pipeline_patient_id"].iloc[0]).strip()
        patient_dir = OUTPUT_DIR / safe_name(patient_id)
        try:
            patient = preprocess_new_patient(
                raw_patient_df.drop(columns=["_pipeline_patient_id"], errors="ignore")
            )
            references = retrieve_similar_examples(
                patient, reference_library, top_k=NUMBER_OF_REFERENCE_EXAMPLES
            )
            write_json(patient_dir / "preprocessed_input.json", patient)
            write_json(patient_dir / "reference_examples.json", references)

            worker, auditor = await run_independent_agents(patient, references)
            comparison = compare_results(worker, auditor)
            write_json(patient_dir / "worker.json", worker.model_dump(mode="json"))
            write_json(patient_dir / "auditor.json", auditor.model_dump(mode="json"))
            write_json(patient_dir / "comparison.json", comparison)

            critique = None
            if not comparison["exact_mapping_agreement"]:
                critique = await run_critique(patient, references, worker, auditor, comparison)
                write_json(patient_dir / "critique.json", critique.model_dump(mode="json"))

            decision = make_decision(worker, auditor, comparison, critique)
            write_json(patient_dir / "decision.json", decision.model_dump(mode="json"))
            selected = worker if decision.selected_result == "worker" else auditor if decision.selected_result == "auditor" else None
            if selected:
                write_json(patient_dir / "selected_result.json", selected.model_dump(mode="json"))

            return {
                "patient_id": patient_id,
                "status": decision.status,
                "selected_result": decision.selected_result,
                "vendor_line_count": len(patient["vendor_lines"]),
                "worker_corrected_line_count": len({line.adjudicated_lot for line in worker.corrected_lines}),
                "auditor_corrected_line_count": len({line.adjudicated_lot for line in auditor.corrected_lines}),
                "worker_confidence": worker.overall_confidence,
                "auditor_confidence": auditor.overall_confidence,
                "exact_mapping_agreement": comparison["exact_mapping_agreement"],
                "major_disagreement": comparison["major_disagreement"],
                "reason": decision.reason,
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
                "reason": "",
                "error": str(error),
            }


async def main() -> None:
    new_data = load_new_data()
    reference_library = build_reference_library(OLD_ADJUDICATED_FILE, OLD_COTA_SHEET)
    groups = [
        group.copy()
        for _, group in new_data.groupby("_pipeline_patient_id", sort=False, dropna=False)
    ]
    print(f"Historical reference patients: {len(reference_library.patients)}")
    print(f"New patients to process: {len(groups)}")

    semaphore = asyncio.Semaphore(MAX_CONCURRENT_PATIENTS)
    summaries = await asyncio.gather(
        *(process_patient(group, reference_library, semaphore) for group in groups)
    )
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUTPUT_DIR / "adjudication_summary.csv", index=False)
    summary[summary["status"] == "human_review"].to_csv(
        OUTPUT_DIR / "human_review_queue.csv", index=False
    )
    print(summary["status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    asyncio.run(main())
