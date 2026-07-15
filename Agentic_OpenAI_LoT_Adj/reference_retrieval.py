from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

from config import (
    DISCONTINUE_REASON_COLUMN,
    END_DATE_COLUMN,
    PATIENT_ID_COLUMN,
    REFERENCE_LABEL_COLUMN,
    REFERENCE_RESIDUAL_COLUMN,
    REGIMEN_COLUMN,
    START_DATE_COLUMN,
    VENDOR_LOT_COLUMN,
)
from preprocessing import clean_text, extract_drugs, format_date


@dataclass
class ReferenceLibrary:
    patients: list[dict[str, Any]]
    documents: list[str]
    vectorizer: TfidfVectorizer
    vectors: Any


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(column).strip() for column in df.columns]
    return df


def collapse_historical_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    for _, row in df.iterrows():
        patient_id = clean_text(row.get(PATIENT_ID_COLUMN))
        vendor_lot = pd.to_numeric(row.get(VENDOR_LOT_COLUMN), errors="coerce")
        regimen = clean_text(row.get(REGIMEN_COLUMN))

        if patient_id and not pd.isna(vendor_lot):
            if current is not None:
                rows.append(current)
            current = row.to_dict()
            current[PATIENT_ID_COLUMN] = patient_id
            current[VENDOR_LOT_COLUMN] = int(vendor_lot)
        elif current is not None and regimen and pd.isna(vendor_lot):
            existing = clean_text(current.get(REGIMEN_COLUMN))
            current[REGIMEN_COLUMN] = f"{existing}, {regimen}" if existing else regimen

    if current is not None:
        rows.append(current)
    return pd.DataFrame(rows)


def determine_action(previous_corrected_lot: int | None, corrected_lot: int) -> str:
    if previous_corrected_lot is None:
        return "keep"
    if corrected_lot == previous_corrected_lot:
        return "merge_with_previous"
    if corrected_lot == previous_corrected_lot + 1:
        return "keep"
    return "historical_label_anomaly"


def build_patient_reference(patient_id: str, patient_df: pd.DataFrame) -> dict[str, Any]:
    lines: list[dict[str, Any]] = []
    previous_corrected: int | None = None

    for _, row in patient_df.sort_values(VENDOR_LOT_COLUMN).iterrows():
        vendor_lot = int(row[VENDOR_LOT_COLUMN])
        corrected_lot = int(row[REFERENCE_LABEL_COLUMN])
        action = determine_action(previous_corrected, corrected_lot)
        lines.append(
            {
                "vendor_lot": vendor_lot,
                "adjudicated_lot": corrected_lot,
                "action": action,
                "regimen": extract_drugs(row.get(REGIMEN_COLUMN)),
                "raw_regimen": clean_text(row.get(REGIMEN_COLUMN)),
                "start_date": format_date(row.get(START_DATE_COLUMN)),
                "end_date": format_date(row.get(END_DATE_COLUMN)),
                "discontinue_reason": clean_text(row.get(DISCONTINUE_REASON_COLUMN)),
                "residual": int(pd.to_numeric(row.get(REFERENCE_RESIDUAL_COLUMN), errors="coerce") or 0),
            }
        )
        previous_corrected = corrected_lot

    merged = [line["vendor_lot"] for line in lines if line["action"] == "merge_with_previous"]
    return {
        "reference_patient_id": patient_id,
        "vendor_line_count": len(lines),
        "adjudicated_line_count": len({line["adjudicated_lot"] for line in lines}),
        "contains_oversplit_correction": bool(merged),
        "merged_vendor_lines": merged,
        "lines": lines,
    }


def transition_tokens(lines: list[dict[str, Any]]) -> list[str]:
    tokens: list[str] = []
    steroids = {"dexamethasone", "prednisone", "prednisolone", "methylprednisolone"}

    for previous, current in zip(lines, lines[1:]):
        before = set(previous["regimen"])
        after = set(current["regimen"])
        added, removed, retained = after - before, before - after, before & after

        tokens.extend(f"added_{drug}" for drug in sorted(added))
        tokens.extend(f"removed_{drug}" for drug in sorted(removed))
        tokens.extend(f"retained_{drug}" for drug in sorted(retained))
        if before == after:
            tokens.append("identical_adjacent_regimen")
        if after < before:
            tokens.append("regimen_reduction")
        if before < after:
            tokens.append("regimen_expansion")
        changed = added | removed
        if changed and changed.issubset(steroids):
            tokens.append("steroid_only_change")
    return tokens


def patient_document(patient: dict[str, Any]) -> str:
    parts = [f"line_count_{patient['vendor_line_count']}"]
    for line in patient["lines"]:
        lot = line["vendor_lot"]
        parts.append(f"vendor_lot_{lot}")
        for drug in line["regimen"]:
            parts.extend([drug, f"lot_{lot}_{drug}"])
    parts.extend(transition_tokens(patient["lines"]))
    return " ".join(parts)


def build_reference_library(workbook_path: str | Path, sheet_name: str = "Cota") -> ReferenceLibrary:
    workbook_path = Path(workbook_path)
    if not workbook_path.exists():
        raise FileNotFoundError(f"Old adjudicated workbook not found: {workbook_path}")

    raw = normalize_columns(pd.read_excel(workbook_path, sheet_name=sheet_name))
    required = {
        PATIENT_ID_COLUMN,
        VENDOR_LOT_COLUMN,
        REGIMEN_COLUMN,
        START_DATE_COLUMN,
        END_DATE_COLUMN,
        REFERENCE_LABEL_COLUMN,
    }
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"Old COTA sheet is missing columns: {sorted(missing)}")

    collapsed = collapse_historical_rows(raw)
    collapsed[REFERENCE_LABEL_COLUMN] = pd.to_numeric(
        collapsed[REFERENCE_LABEL_COLUMN], errors="coerce"
    )
    collapsed = collapsed.dropna(subset=[REFERENCE_LABEL_COLUMN])
    collapsed[REFERENCE_LABEL_COLUMN] = collapsed[REFERENCE_LABEL_COLUMN].astype(int)

    patients = [
        build_patient_reference(str(patient_id), patient_df)
        for patient_id, patient_df in collapsed.groupby(PATIENT_ID_COLUMN, sort=False)
    ]
    if not patients:
        raise ValueError("No labeled historical COTA patients were found.")

    documents = [patient_document(patient) for patient in patients]
    vectorizer = TfidfVectorizer(
        lowercase=True,
        token_pattern=r"(?u)\b[\w_-]+\b",
        ngram_range=(1, 2),
        min_df=1,
    )
    vectors = vectorizer.fit_transform(documents)
    return ReferenceLibrary(patients, documents, vectorizer, vectors)


def convert_new_patient(patient: dict[str, Any]) -> dict[str, Any]:
    return {
        "reference_patient_id": patient["patient_id"],
        "vendor_line_count": len(patient["vendor_lines"]),
        "lines": [
            {
                "vendor_lot": line["vendor_lot"],
                "regimen": line["regimen"],
            }
            for line in patient["vendor_lines"]
        ],
    }


def retrieve_similar_examples(
    patient: dict[str, Any],
    reference_library: ReferenceLibrary,
    top_k: int = 5,
    oversplit_bonus: float = 0.05,
) -> list[dict[str, Any]]:
    query = patient_document(convert_new_patient(patient))
    query_vector = reference_library.vectorizer.transform([query])
    similarities = cosine_similarity(query_vector, reference_library.vectors)[0]

    ranked: list[tuple[float, float, int]] = []
    for index, similarity in enumerate(similarities):
        historical = reference_library.patients[index]
        bonus = oversplit_bonus if historical["contains_oversplit_correction"] else 0.0
        ranked.append((float(similarity + bonus), float(similarity), index))
    ranked.sort(reverse=True)

    results: list[dict[str, Any]] = []
    for _, similarity, index in ranked[: max(1, top_k)]:
        item = dict(reference_library.patients[index])
        item["similarity_score"] = round(similarity, 4)
        results.append(item)
    return results


def build_reference_library_from_dataframe(
    dataframe: pd.DataFrame,
) -> ReferenceLibrary:
    """Build a leakage-safe reference library from an already filtered DataFrame."""
    raw = normalize_columns(dataframe)
    required = {
        PATIENT_ID_COLUMN,
        VENDOR_LOT_COLUMN,
        REGIMEN_COLUMN,
        START_DATE_COLUMN,
        END_DATE_COLUMN,
        REFERENCE_LABEL_COLUMN,
    }
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"Old COTA data is missing columns: {sorted(missing)}")

    collapsed = collapse_historical_rows(raw)
    collapsed[REFERENCE_LABEL_COLUMN] = pd.to_numeric(
        collapsed[REFERENCE_LABEL_COLUMN], errors="coerce"
    )
    collapsed = collapsed.dropna(
        subset=[PATIENT_ID_COLUMN, VENDOR_LOT_COLUMN, REFERENCE_LABEL_COLUMN]
    )
    collapsed[PATIENT_ID_COLUMN] = (
        collapsed[PATIENT_ID_COLUMN].astype(str).str.strip()
    )
    collapsed[VENDOR_LOT_COLUMN] = pd.to_numeric(
        collapsed[VENDOR_LOT_COLUMN], errors="coerce"
    ).astype(int)
    collapsed[REFERENCE_LABEL_COLUMN] = (
        collapsed[REFERENCE_LABEL_COLUMN].astype(int)
    )

    patients = [
        build_patient_reference(str(patient_id), patient_df)
        for patient_id, patient_df in collapsed.groupby(
            PATIENT_ID_COLUMN, sort=False
        )
    ]
    if not patients:
        raise ValueError("No labeled historical COTA patients were found.")

    documents = [patient_document(patient) for patient in patients]
    vectorizer = TfidfVectorizer(
        lowercase=True,
        token_pattern=r"(?u)\b[\w_-]+\b",
        ngram_range=(1, 2),
        min_df=1,
    )
    vectors = vectorizer.fit_transform(documents)
    return ReferenceLibrary(patients, documents, vectorizer, vectors)


def retrieve_similar_examples_excluding(
    patient: dict[str, Any],
    reference_library: ReferenceLibrary,
    top_k: int = 5,
    excluded_patient_ids: set[str] | None = None,
    oversplit_bonus: float = 0.05,
) -> list[dict[str, Any]]:
    """Retrieve similar examples while excluding held-out evaluation patients."""
    excluded = {
        str(patient_id).strip()
        for patient_id in (excluded_patient_ids or set())
    }

    query = patient_document(convert_new_patient(patient))
    query_vector = reference_library.vectorizer.transform([query])
    similarities = cosine_similarity(
        query_vector, reference_library.vectors
    )[0]

    ranked: list[tuple[float, float, int]] = []
    for index, similarity in enumerate(similarities):
        historical = reference_library.patients[index]
        reference_id = str(
            historical["reference_patient_id"]
        ).strip()
        if reference_id in excluded:
            continue

        bonus = (
            oversplit_bonus
            if historical["contains_oversplit_correction"]
            else 0.0
        )
        ranked.append(
            (float(similarity + bonus), float(similarity), index)
        )

    ranked.sort(reverse=True)

    results: list[dict[str, Any]] = []
    for _, similarity, index in ranked[: max(1, top_k)]:
        item = dict(reference_library.patients[index])
        item["similarity_score"] = round(similarity, 4)
        results.append(item)

    return results