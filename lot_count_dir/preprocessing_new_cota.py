"""
Preprocessing for new_cota_data.xlsx (or any similarly-shaped raw COTA
export) -> Output/COTA_cleaned.xlsx, ready for lot_counting_fixed.py
"""

from __future__ import annotations

from pathlib import Path
import re
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR.parent / "data" / "new_cota_data.xlsx"
OUTPUT_PATH = BASE_DIR / "Output" / "COTA_cleaned.xlsx"

# Columns that only ever appear once per PATIENT (not per LOT row) -- forward
# filled at the patient level before any LOT-grouping happens.
PATIENT_LEVEL_COLUMNS = ["cpid", "diag_dt", "deathfl", "dthdt_c"]

# Columns kept in the final output. Everything else (date_of_death_imp,
# date_start_line_of_therapy_imp, date_end_line_of_therapy_imp, aval,
# avaldt) is dropped as unused by lot_counting_fixed.py.
KEPT_RAW_COLUMNS = [
    "cpid", "line_of_therapy_c", "line_of_therapy_name", "discontinue_reason",
    "date_start_line_of_therapy", "date_end_line_of_therapy",
]

STEROIDS = {"dexamethasone", "prednisone", "prednisolone", "methylprednisolone"}
CONDITIONING = {"melphalan", "busulfan", "carmustine"}
PROCEDURE_TOKENS = {"autologous sct", "allogeneic sct", "sct", "stem cell transplant", "transplant"}
CAR_T_TOKENS = {"cart", "car-t", "car t", "chimeric antigen receptor"}
# Specific CAR-T product names seen in practice -- extend as needed.
CAR_T_PRODUCT_TOKENS = {"cilta-cel", "ciltacabtagene", "ide-cel", "idecabtagene", "abecma", "carvykti"}

# Best-effort drug -> family/class map (see caveat in module docstring).
DRUG_FAMILY_MAP = {
    "bortezomib": "proteasome_inhibitor", "carfilzomib": "proteasome_inhibitor", "ixazomib": "proteasome_inhibitor",
    "lenalidomide": "imid", "thalidomide": "imid", "pomalidomide": "imid",
    "daratumumab": "anti_cd38", "isatuximab": "anti_cd38",
    "talquetamab": "bispecific", "teclistamab": "bispecific", "elranatamab": "bispecific",
    "belantamab": "adc",
    "cyclophosphamide": "alkylator", "melphalan": "alkylator",
    "etoposide": "chemo_other", "selinexor": "sine",
    "dexamethasone": "steroid", "prednisone": "steroid", "prednisolone": "steroid", "methylprednisolone": "steroid",
}


def normalize_text(value) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def clean_token(value: str) -> str:
    token = normalize_text(value).lower()
    token = token.strip().strip("[]'").strip('"')
    return re.sub(r"\s+", " ", token.replace("[", "").replace("]", "")).strip()


def parse_regimen_chunks(value) -> list[tuple[str, ...]]:
    text = normalize_text(value).lower()
    if not text:
        return []
    bracket_chunks = re.findall(r"\[([^\]]+)\]", text)
    if not bracket_chunks:
        bracket_chunks = [text.replace("[", "").replace("]", "")]
    chunks = []
    for chunk in bracket_chunks:
        drugs = []
        for raw in chunk.split(","):
            drug = clean_token(raw)
            if drug and drug not in drugs:
                drugs.append(drug)
        if drugs:
            chunks.append(tuple(drugs))
    return chunks


def load_raw(path) -> pd.DataFrame:
    df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def merge_wrapped_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapses physically-wrapped rows into single logical LOT records.

    A new LOT group starts at every row with a non-null line_of_therapy_c.
    Every following row (until the next non-null line_of_therapy_c, or a new
    patient) is a continuation of that group -- its line_of_therapy_name and
    discontinue_reason text (if any) are appended to the group, rather than
    being read as a separate row or silently dropped.
    """
    df = df.copy()

    # Patient-level fields are only populated on each patient's first row.
    for col in PATIENT_LEVEL_COLUMNS:
        if col in df.columns:
            df[col] = df[col].ffill()

    df["_original_row_order"] = df.index

    # Every non-null line_of_therapy_c starts a new group; NaNs inherit the
    # most recent group. (Safe across patient boundaries: every patient's
    # first row always has a non-null line_of_therapy_c -- verified against
    # the real file, no exceptions.)
    df["_lot_group_id"] = df["line_of_therapy_c"].notna().cumsum()

    records = []
    for _, group in df.groupby("_lot_group_id", sort=True):
        anchor = group.iloc[0]

        name_fragments = [normalize_text(v) for v in group["line_of_therapy_name"] if normalize_text(v)]
        reason_fragments = [normalize_text(v) for v in group["discontinue_reason"] if normalize_text(v)]

        record = {
            "cpid": anchor.get("cpid"),
            "line_of_therapy_c": anchor.get("line_of_therapy_c"),
            "line_of_therapy_name": " ".join(name_fragments) if name_fragments else None,
            "discontinue_reason": "".join(reason_fragments) if reason_fragments else None,
            "date_start_line_of_therapy": group["date_start_line_of_therapy"].dropna().iloc[0]
                if group["date_start_line_of_therapy"].notna().any() else None,
            "date_end_line_of_therapy": group["date_end_line_of_therapy"].dropna().iloc[0]
                if group["date_end_line_of_therapy"].notna().any() else None,
            "_original_row_order": int(anchor["_original_row_order"]),
        }
        records.append(record)

    return pd.DataFrame(records)


def compute_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    chunks_col = df["line_of_therapy_name"].apply(parse_regimen_chunks)
    all_tokens = chunks_col.apply(lambda chunks: {t for chunk in chunks for t in chunk})
    tail_tokens = chunks_col.apply(lambda chunks: set(chunks[-1]) if chunks else set())

    df["is_asct"] = all_tokens.apply(lambda toks: bool(toks & PROCEDURE_TOKENS))
    df["is_car_t"] = all_tokens.apply(lambda toks: bool(toks & CAR_T_TOKENS) or any(
        any(prod in t for prod in CAR_T_PRODUCT_TOKENS) for t in toks
    ))

    def active_only(tokens: set[str]) -> set[str]:
        return {t for t in tokens if t not in STEROIDS and t not in CONDITIONING
                and t not in PROCEDURE_TOKENS and t not in CAR_T_TOKENS}

    df["drugs_active"] = all_tokens.apply(lambda toks: ", ".join(sorted(active_only(toks))))
    df["tail_drugs_active"] = tail_tokens.apply(lambda toks: ", ".join(sorted(active_only(toks))))
    df["is_steroid_only"] = all_tokens.apply(lambda toks: bool(toks) and toks.issubset(STEROIDS))

    def family_classes(tokens: set[str]) -> str:
        classes = {DRUG_FAMILY_MAP[t] for t in tokens if t in DRUG_FAMILY_MAP}
        return ", ".join(sorted(classes))

    df["mapped_family_classes"] = all_tokens.apply(family_classes)

    start = pd.to_datetime(df["date_start_line_of_therapy"], errors="coerce")
    end = pd.to_datetime(df["date_end_line_of_therapy"], errors="coerce")
    df["segment_duration_days"] = (end - start).dt.days

    # asct_index: running count of ASCT occurrences seen so far, per patient,
    # forward-filled so the row(s) right after an ASCT can see "which ASCT
    # number just happened" (matches how rule_asct_rule reads it: checks
    # asct_index == 1 for the first-transplant maintenance transition, etc.)
    df = df.sort_values(["cpid", "_original_row_order"])
    df["asct_index"] = (
        df.groupby("cpid")["is_asct"]
        .apply(lambda s: s.cumsum().where(s.cumsum() > 0))
        .reset_index(level=0, drop=True)
    )

    return df


def preprocess(input_path=INPUT_PATH, output_path=OUTPUT_PATH) -> pd.DataFrame:
    raw = load_raw(input_path)
    merged = merge_wrapped_rows(raw)
    final = compute_derived_columns(merged)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path) as writer:
        final.to_excel(writer, sheet_name="All_Data", index=False)

    print(f"Input:  {input_path} ({len(raw)} raw rows)")
    print(f"Output: {output_path} ({len(final)} merged LOT records, "
          f"{final['cpid'].nunique()} patients)")
    print(f"  is_asct: {int(final['is_asct'].sum())} rows")
    print(f"  is_car_t: {int(final['is_car_t'].sum())} rows")
    print(f"  is_steroid_only: {int(final['is_steroid_only'].sum())} rows")
    return final


if __name__ == "__main__":
    preprocess()
