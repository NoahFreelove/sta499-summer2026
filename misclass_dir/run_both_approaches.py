"""
Run BOTH LoT approaches and write their results into one output file.

1. Deterministic approach: algo_tester.py (rule-engine LoT counter).
   Produces Output/lot_counts_with_rules.xlsx.
2. AI approach: the Agentic_OpenAI_LoT_Adj worker/auditor pipeline
   (lives on branch agentic-flow-worker-auditor; fetch it into this branch with
      git checkout origin/agentic-flow-worker-auditor -- Agentic_OpenAI_LoT_Adj
   and set OPENAI_API_KEY to run it).
   Produces Agentic_OpenAI_LoT_Adj/output/adjudication_summary.csv.

Both are combined into Output/lot_counts_combined.xlsx:
  - "LoT_Counts_Long": one row per (patient, approach) - each approach's count
    appears as additional rows in the same file.
  - "LoT_Counts_Wide": one row per patient with a column per approach plus
    differences vs the doctor's adjudicated total.

Usage:
    python run_both_approaches.py                 # run both (AI if available)
    python run_both_approaches.py --skip-ai       # deterministic only
    python run_both_approaches.py --use-cached    # reuse existing outputs
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
REPO_DIR = BASE_DIR.parent

DETERMINISTIC_SCRIPT = BASE_DIR / "algo_tester.py"
DETERMINISTIC_OUTPUT = BASE_DIR / "Output" / "lot_counts_with_rules.xlsx"
CLEANED_INPUT = BASE_DIR / "Output" / "COTA_cleaned.xlsx"

AI_DIR = Path(os.getenv("AGENTIC_LOT_DIR", REPO_DIR / "Agentic_OpenAI_LoT_Adj"))
AI_SCRIPT = AI_DIR / "pipeline.py"
AI_SUMMARY = AI_DIR / "output" / "adjudication_summary.csv"
AI_INPUT = AI_DIR / "data" / "new_cota_data.xlsx"
AI_REFERENCE = AI_DIR / "data" / "LoT Adjudication Datasets.xlsx"
RAW_ADJUDICATION = REPO_DIR / "data" / "LoT Adjudication Datasets.xlsx"

COMBINED_OUTPUT = BASE_DIR / "Output" / "lot_counts_combined.xlsx"


def run_deterministic(use_cached: bool) -> pd.DataFrame:
    """Run algo_tester.py and return per-patient totals for every approach it knows."""
    if not (use_cached and DETERMINISTIC_OUTPUT.exists()):
        print(f"Running deterministic counter: {DETERMINISTIC_SCRIPT.name} ...")
        subprocess.run(
            [sys.executable, str(DETERMINISTIC_SCRIPT)],
            cwd=BASE_DIR, check=True, capture_output=True, text=True,
        )
    else:
        print(f"Reusing cached deterministic output: {DETERMINISTIC_OUTPUT.name}")

    lot = pd.read_excel(DETERMINISTIC_OUTPUT, sheet_name="Detailed_LoT_Assignments")
    per_patient = (
        lot.groupby("cpid")["total_lots_for_patient"].max()
        .rename("deterministic_rules").reset_index()
    )

    # Doctor-adjudicated and vendor totals come from the preprocessing output
    allr = pd.read_excel(CLEANED_INPUT, sheet_name="All_Data")
    allr.columns = allr.columns.str.strip()
    ref = allr.groupby("cpid").agg(
        doctor=("doctor_lot_numeric_for_transition", "max"),
        cota_vendor=("cota_lot_numeric", "max"),
    ).reset_index()

    return per_patient.merge(ref, on="cpid", how="left")


def prepare_ai_inputs() -> None:
    """Create the data files the agentic pipeline expects, if they are missing."""
    AI_INPUT.parent.mkdir(parents=True, exist_ok=True)
    if not AI_REFERENCE.exists() and RAW_ADJUDICATION.exists():
        import shutil
        shutil.copy(RAW_ADJUDICATION, AI_REFERENCE)
        print(f"Copied reference workbook to {AI_REFERENCE}")
    if not AI_INPUT.exists():
        # The pipeline reads unlabeled COTA rows from Sheet1; feed it the raw
        # Cota sheet of the adjudication workbook (labels are simply unused).
        cota = pd.read_excel(RAW_ADJUDICATION, sheet_name="Cota")
        cota.to_excel(AI_INPUT, sheet_name="Sheet1", index=False)
        print(f"Wrote AI input workbook to {AI_INPUT}")


def run_ai(use_cached: bool, skip_ai: bool) -> pd.DataFrame | None:
    """Run the agentic pipeline (or reuse its summary) and return per-patient counts."""
    if skip_ai:
        print("Skipping AI approach (--skip-ai).")
        return None

    if not AI_DIR.exists():
        print(
            f"AI approach not found at {AI_DIR}.\n"
            "  Fetch it with: git checkout origin/agentic-flow-worker-auditor -- Agentic_OpenAI_LoT_Adj"
        )
        return None

    can_run = AI_SCRIPT.exists() and os.getenv("OPENAI_API_KEY")
    if not (use_cached and AI_SUMMARY.exists()):
        if can_run:
            prepare_ai_inputs()
            print(f"Running AI pipeline: {AI_SCRIPT} ...")
            subprocess.run(
                [sys.executable, str(AI_SCRIPT)],
                cwd=AI_DIR, check=True,
            )
        elif AI_SUMMARY.exists():
            print("OPENAI_API_KEY not set - reusing existing adjudication_summary.csv.")
        else:
            print("Cannot run AI approach (no OPENAI_API_KEY) and no cached summary found - skipping.")
            return None
    else:
        print(f"Reusing cached AI output: {AI_SUMMARY}")

    summary = pd.read_csv(AI_SUMMARY)
    summary = summary.rename(columns={"patient_id": "cpid"})

    def selected_count(row):
        if row.get("selected_result") == "worker":
            return row.get("worker_corrected_line_count")
        if row.get("selected_result") == "auditor":
            return row.get("auditor_corrected_line_count")
        return pd.NA

    out = pd.DataFrame({
        "cpid": summary["cpid"],
        "ai_worker": summary.get("worker_corrected_line_count"),
        "ai_auditor": summary.get("auditor_corrected_line_count"),
        "ai_selected": summary.apply(selected_count, axis=1),
        "ai_status": summary.get("status"),
    })
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--skip-ai", action="store_true",
                        help="only run the deterministic approach")
    parser.add_argument("--use-cached", action="store_true",
                        help="reuse existing outputs instead of re-running the pipelines")
    args = parser.parse_args()

    det = run_deterministic(args.use_cached)
    ai = run_ai(args.use_cached, args.skip_ai)

    wide = det if ai is None else det.merge(ai, on="cpid", how="outer")

    count_cols = ["cota_vendor", "doctor", "deterministic_rules"]
    if ai is not None:
        count_cols += ["ai_worker", "ai_auditor", "ai_selected"]

    # Long format: each approach's count is its own row in the same file
    long = wide.melt(
        id_vars="cpid", value_vars=count_cols,
        var_name="approach", value_name="total_lots",
    ).dropna(subset=["total_lots"]).sort_values(["cpid", "approach"])

    # Wide comparison with differences vs the doctor's adjudicated total
    wide = wide.copy()
    wide["deterministic_minus_doctor"] = wide["deterministic_rules"] - wide["doctor"]
    if ai is not None:
        wide["ai_selected_minus_doctor"] = (
            pd.to_numeric(wide["ai_selected"], errors="coerce") - wide["doctor"]
        )

    with pd.ExcelWriter(COMBINED_OUTPUT) as writer:
        long.to_excel(writer, sheet_name="LoT_Counts_Long", index=False)
        wide.to_excel(writer, sheet_name="LoT_Counts_Wide", index=False)

    print(f"\nSaved combined results to {COMBINED_OUTPUT}")
    print(f"  Patients: {wide['cpid'].nunique()}, approaches: {', '.join(count_cols)}")
    for col in count_cols:
        vals = pd.to_numeric(wide[col], errors="coerce")
        err = vals - wide["doctor"]
        print(f"  {col:>20}: mean LoT {vals.mean():.2f}  "
              f"exact vs doctor {(err == 0).mean():.1%}  MAE {err.abs().mean():.2f}")


if __name__ == "__main__":
    main()
