from __future__ import annotations

import re
from typing import Any

import pandas as pd

from config import (
    DISCONTINUE_REASON_COLUMN,
    END_DATE_COLUMN,
    PATIENT_ID_COLUMN,
    REGIMEN_COLUMN,
    START_DATE_COLUMN,
    VENDOR_LOT_COLUMN,
)


def clean_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def format_date(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None
    parsed = pd.to_datetime(value, errors="coerce")
    if pd.isna(parsed):
        text = clean_text(value)
        return text or None
    return parsed.strftime("%Y-%m-%d")


def extract_drugs(regimen: Any) -> list[str]:
    text = clean_text(regimen).lower()
    if not text:
        return []
    text = text.replace("], [", ",").replace("][", ",")
    text = text.replace("[", "").replace("]", "")
    drugs: list[str] = []
    for part in re.split(r"[,;|]+", text):
        drug = re.sub(r"\s+", " ", part.strip())
        if drug and drug not in drugs:
            drugs.append(drug)
    return sorted(drugs)


def collapse_new_patient_rows(patient_df: pd.DataFrame) -> list[dict[str, Any]]:
    lines: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    for source_index, row in patient_df.iterrows():
        vendor_lot = pd.to_numeric(row.get(VENDOR_LOT_COLUMN), errors="coerce")
        regimen_text = clean_text(row.get(REGIMEN_COLUMN))

        if not pd.isna(vendor_lot):
            if current is not None:
                lines.append(current)
            current = {
                "vendor_lot": int(vendor_lot),
                "raw_regimen": regimen_text,
                "start_date": format_date(row.get(START_DATE_COLUMN)),
                "end_date": format_date(row.get(END_DATE_COLUMN)),
                "discontinue_reason": clean_text(row.get(DISCONTINUE_REASON_COLUMN)),
                "source_row_indices": [str(source_index)],
            }
        elif current is not None and regimen_text:
            current["raw_regimen"] = (
                f"{current['raw_regimen']}, {regimen_text}" if current["raw_regimen"] else regimen_text
            )
            current["source_row_indices"].append(str(source_index))

    if current is not None:
        lines.append(current)

    for line in lines:
        line["regimen"] = extract_drugs(line["raw_regimen"])

    return lines


def preprocess_new_patient(patient_df: pd.DataFrame) -> dict[str, Any]:
    required = {
        PATIENT_ID_COLUMN,
        VENDOR_LOT_COLUMN,
        REGIMEN_COLUMN,
        START_DATE_COLUMN,
        END_DATE_COLUMN,
    }
    missing = required - set(patient_df.columns)
    if missing:
        raise ValueError(f"New COTA data is missing columns: {sorted(missing)}")

    patient_ids = (
        patient_df[PATIENT_ID_COLUMN]
        .replace(r"^\s*$", pd.NA, regex=True)
        .dropna()
        .astype(str)
        .str.strip()
        .unique()
    )
    if len(patient_ids) != 1:
        raise ValueError(f"Expected exactly one patient, found: {patient_ids.tolist()}")

    vendor_lines = collapse_new_patient_rows(patient_df)
    if not vendor_lines:
        raise ValueError("No vendor LoT rows were found for the patient.")

    actual_lots = [line["vendor_lot"] for line in vendor_lines]
    if actual_lots != sorted(actual_lots) or len(actual_lots) != len(set(actual_lots)):
        raise ValueError(f"Vendor LoT values must be unique and sorted: {actual_lots}")

    return {
        "patient_id": patient_ids[0],
        "vendor_line_count": len(vendor_lines),
        "vendor_lines": vendor_lines,
    }
