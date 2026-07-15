from pathlib import Path
import pandas as pd
import numpy as np
import re

BASE_DIR = Path(__file__).resolve().parent

input_path = BASE_DIR / "Output" / "COTA_family_reconstructed.xlsx"
output_path = BASE_DIR / "Output" / "COTA_misclassified_rows_UDP.xlsx"

df = pd.read_excel(input_path, sheet_name="COTA_Reconstructed")
df.columns = df.columns.str.strip()

print("Reading file:", input_path)
print("Rows loaded before split-row repair:", len(df))


def is_missing(value):
    if pd.isna(value):
        return True
    value = str(value).strip()
    return value in ["", ".", "N/A", "NA", "nan", "None", "Missing"]


def clean_text(value):
    if is_missing(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def to_number(value):
    if is_missing(value):
        return np.nan
    try:
        return float(str(value).strip())
    except ValueError:
        return np.nan


def unique_join(values, sep=" | "):
    seen = []
    for value in values:
        value = clean_text(value)
        if value and value not in seen:
            seen.append(value)
    return sep.join(seen)


def first_non_missing(values):
    for value in values:
        if not is_missing(value):
            return value
    return np.nan


def merge_regimen_fragments(values):
    """
    Merge split line_of_therapy_name fragments belonging to the same patient + COTA LoT.
    Example:
      [bortezomib, dexamethasone,
      lenalidomide], [dexamethasone,
      pomalidomide], [dexamethasone]
    becomes one full regimen string.
    """
    parts = [clean_text(v) for v in values if clean_text(v)]
    merged = " ".join(parts)
    merged = re.sub(r"\s+", " ", merged).strip()
    merged = merged.replace("[ ", "[").replace(" ]", "]")
    return merged


def parse_regimen_drugs(regimen):
    """
    Parse drugs after merged regimen is reconstructed.
    Handles normal bracketed lists and imperfect leftover fragments.
    """
    text = clean_text(regimen).lower()
    if not text:
        return []

    drugs = []

    # Prefer complete bracket chunks.
    bracket_chunks = re.findall(r"\[([^\]]+)\]", text)

    if bracket_chunks:
        chunks = bracket_chunks
    else:
        # Fallback if brackets are still imperfect.
        text = text.replace("[", "").replace("]", "")
        chunks = [text]

    for chunk in chunks:
        for raw_drug in chunk.split(","):
            drug = raw_drug.strip().strip("'\"")
            drug = drug.replace("[", "").replace("]", "")
            drug = re.sub(r"\s+", " ", drug).strip()

            if drug and drug not in drugs:
                drugs.append(drug)

    return drugs


# ---------------------------------------------------------------------
# STEP 1: Repair split / continuation rows
# ---------------------------------------------------------------------

df["_original_row_order"] = range(len(df))

# Keep original values so we know what was blank before forward-fill.
df["_original_cpid"] = df["cpid"] if "cpid" in df.columns else np.nan
df["_original_line_of_therapy_c"] = df["line_of_therapy_c"] if "line_of_therapy_c" in df.columns else np.nan

# Treat "Missing" as missing.
for col in ["cpid", "line_of_therapy_c", "Alpesh 1 LoT", "Alberto LOT"]:
    if col in df.columns:
        df[col] = df[col].apply(lambda x: np.nan if is_missing(x) else x)

# Forward-fill identifiers/reviewer LoTs so split fragments inherit the parent row values.
fill_down_cols = [
    "cpid",
    "line_of_therapy_c",
    "Alpesh 1 LoT",
    "Alberto LOT",
    "Reviewer 1",
    "Reviewer 2",
]

for col in fill_down_cols:
    if col in df.columns:
        df[col] = df[col].ffill()

# Mark rows that were likely continuation fragments.
df["was_split_continuation_row"] = (
    df["_original_cpid"].apply(is_missing)
    | df["_original_line_of_therapy_c"].apply(is_missing)
)

# Group by patient + COTA LoT. This is the key repair.
group_keys = ["cpid", "line_of_therapy_c"]

rows = []

for _, group in df.groupby(group_keys, dropna=False, sort=False):
    group = group.sort_values("_original_row_order")

    row = group.iloc[0].copy()

    # Merge the regimen fragments.
    row["line_of_therapy_name"] = merge_regimen_fragments(group["line_of_therapy_name"])

    # Re-parse drugs from repaired regimen.
    repaired_drugs = parse_regimen_drugs(row["line_of_therapy_name"])
    row["parsed_drugs_repaired"] = ", ".join(repaired_drugs)

    # Preserve evidence that this row was merged.
    row["merged_from_split_rows"] = bool(group["was_split_continuation_row"].any())
    row["merged_original_row_count"] = len(group)
    row["merged_original_row_orders"] = ", ".join(
        str(int(x)) for x in group["_original_row_order"].tolist()
    )

    # Merge useful text fields.
    for col in ["discontinue_reason", "aval", "avaldt"]:
        if col in group.columns:
            row[col] = unique_join(group[col])

    # For reviewer/date fields, keep first non-missing value.
    for col in [
        "diag_dt",
        "deathfl",
        "dthdt_c",
        "refdt",
        "date_start_line_of_therapy",
        "date_end_line_of_therapy",
        "Reviewer 1",
        "Reviewer 1 LoT",
        "Reviewer 2",
        "Reviewer 2 LoT",
        "Alpesh 1 LoT",
        "Alberto LOT",
        "residual",
    ]:
        if col in group.columns:
            row[col] = first_non_missing(group[col])

    rows.append(row)

df = pd.DataFrame(rows)

print("Rows loaded after split-row repair:", len(df))


# ---------------------------------------------------------------------
# STEP 2: Transition logic after split-row repair
# ---------------------------------------------------------------------

df["cota_lot_numeric"] = df["line_of_therapy_c"].apply(to_number)
df["reviewer1_lot_numeric"] = df["Alpesh 1 LoT"].apply(to_number)
df["reviewer2_lot_numeric"] = df["Alberto LOT"].apply(to_number)

# Use Alpesh first, then Alberto if Alpesh is missing.
df["doctor_lot_numeric_for_transition"] = df["reviewer1_lot_numeric"]
df.loc[
    df["doctor_lot_numeric_for_transition"].isna(),
    "doctor_lot_numeric_for_transition"
] = df["reviewer2_lot_numeric"]

# Sort by patient and original row order.
df = df.sort_values(["cpid", "_original_row_order"]).copy()

df["previous_cota_lot_numeric"] = df.groupby("cpid")["cota_lot_numeric"].shift(1)
df["previous_doctor_lot_numeric"] = df.groupby("cpid")[
    "doctor_lot_numeric_for_transition"
].shift(1)

df["is_first_patient_row"] = df.groupby("cpid").cumcount().eq(0)

df["cota_started_new_lot_vs_previous_row"] = (
    df["cota_lot_numeric"] != df["previous_cota_lot_numeric"]
)

df["doctor_started_new_lot_vs_previous_row"] = (
    df["doctor_lot_numeric_for_transition"] != df["previous_doctor_lot_numeric"]
)

# First row per patient is baseline/fine.
df.loc[df["is_first_patient_row"], "cota_started_new_lot_vs_previous_row"] = True
df.loc[df["is_first_patient_row"], "doctor_started_new_lot_vs_previous_row"] = True

missing_required = (
    df["cota_lot_numeric"].isna()
    | df["doctor_lot_numeric_for_transition"].isna()
)

df["transition_alignment_status"] = "transition_aligned"

df.loc[
    missing_required,
    "transition_alignment_status"
] = "needs_context_missing_or_unparseable"

df.loc[
    (~missing_required)
    & (~df["is_first_patient_row"])
    & (df["cota_started_new_lot_vs_previous_row"])
    & (~df["doctor_started_new_lot_vs_previous_row"]),
    "transition_alignment_status"
] = "candidate_misclassification_cota_over_split"

df.loc[
    (~missing_required)
    & (~df["is_first_patient_row"])
    & (~df["cota_started_new_lot_vs_previous_row"])
    & (df["doctor_started_new_lot_vs_previous_row"]),
    "transition_alignment_status"
] = "candidate_misclassification_cota_under_split"

df.loc[
    df["is_first_patient_row"] & (~missing_required),
    "transition_alignment_status"
] = "baseline_first_lot_row"

df["cota_minus_doctor_lot_shift"] = (
    df["cota_lot_numeric"] - df["doctor_lot_numeric_for_transition"]
)


# ---------------------------------------------------------------------
# STEP 3: Extract output tables
# ---------------------------------------------------------------------

misclassified = df[
    df["transition_alignment_status"].isin([
        "candidate_misclassification_cota_over_split",
        "candidate_misclassification_cota_under_split",
    ])
].copy()

needs_context = df[
    df["transition_alignment_status"].eq("needs_context_missing_or_unparseable")
].copy()

preferred_columns = [
    "cpid",
    "line_of_therapy_c",
    "line_of_therapy_name",
    "parsed_drugs_repaired",
    "Alpesh 1 LoT",
    "Alberto LOT",
    "doctor_lot_numeric_for_transition",
    "previous_cota_lot_numeric",
    "cota_lot_numeric",
    "previous_doctor_lot_numeric",
    "cota_started_new_lot_vs_previous_row",
    "doctor_started_new_lot_vs_previous_row",
    "transition_alignment_status",
    "cota_minus_doctor_lot_shift",
    "reconstructed_family_combination",
    "in_provided_fiona_category_list",
    "parsed_drugs",
    "unknown_or_unmapped_drugs",
    "merged_from_split_rows",
    "merged_original_row_count",
    "merged_original_row_orders",
]

preferred_columns = [col for col in preferred_columns if col in df.columns]

misclassified = misclassified[
    preferred_columns + [col for col in misclassified.columns if col not in preferred_columns]
]

needs_context = needs_context[
    preferred_columns + [col for col in needs_context.columns if col not in preferred_columns]
]

summary = (
    misclassified
    .groupby(
        ["transition_alignment_status", "reconstructed_family_combination"],
        dropna=False
    )
    .size()
    .reset_index(name="misclassified_row_count")
    .sort_values("misclassified_row_count", ascending=False)
)

patient_summary = (
    misclassified
    .groupby("cpid", dropna=False)
    .size()
    .reset_index(name="misclassified_row_count")
    .sort_values("misclassified_row_count", ascending=False)
)

split_repair_audit = df[
    df["merged_from_split_rows"].eq(True)
].copy()

with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
    misclassified.to_excel(writer, sheet_name="Misclassified_Rows", index=False)
    summary.to_excel(writer, sheet_name="Misclassification_Summary", index=False)
    patient_summary.to_excel(writer, sheet_name="Patient_Summary", index=False)
    needs_context.to_excel(writer, sheet_name="Needs_Context", index=False)
    split_repair_audit.to_excel(writer, sheet_name="Split_Repair_Audit", index=False)
    df.to_excel(writer, sheet_name="All_COTA_With_Transitions", index=False)

print("Saved:", output_path)
print("Misclassified rows:", len(misclassified))
print("Needs context rows:", len(needs_context))
print("Split-repaired LoT rows:", len(split_repair_audit))
print(df["transition_alignment_status"].value_counts(dropna=False))