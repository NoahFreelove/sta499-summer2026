"""
Combines 2 pre-processing scripts from Ana for simplified usage.

Updates to be incorporated:
- New pre-processing rules from the document should be included.
- SCRIPT 2 output contains many sheets.
  Is COTA_reconstructed sheet still necessary if we're using the fixed ALL_COTA_WITH_TRANSITIONS sheet?
  Consider removing it.
- for SCRIPT 3, we'll need to change the output... to be in the original file?
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd
import numpy as np
import os

"""
SCRIPT 1: 
  - original file: reconstruct_cota_family_combinations.py
  - purpose: reconstruct COTA family combinations using Fiona's provided family-combination list.

Inputs:
  1) Raw LoT adjudication workbook, sheet: Cota
  2) Fiona/COTA misclassification workbook, sheet: Line_by_Line

Main logic:
  - Read all unique family combinations from Fiona's Line_by_Line sheet.
  - Parse COTA line_of_therapy_name into individual drugs.
  - Map drugs to drug-family classes.
  - Reconstruct the family combination in the same class order used by Fiona.
  - If the reconstructed combination is not in Fiona's provided list, keep the
    reconstructed pattern but append: (Not in provided Fiona's category list)
"""

DEFAULT_CLASS_ORDER = [
    "Alkylating agents",
    "BCL-2 inhibitor",
    "CAR-T",
    "HDAC inhibitor",
    "IMiDs",
    "Monoclonal antibodies",
    "Other chemotherapy",
    "Proteasome inhibitors",
    "Steroids",
    "Transplant",
]

# Edit this dictionary as needed after clinical review.
# Keys must be lowercase because parsed drug names are normalized to lowercase.
DRUG_TO_FAMILY = {
    # IMiDs
    "lenalidomide": "IMiDs",
    "pomalidomide": "IMiDs",
    "thalidomide": "IMiDs",

    # Proteasome inhibitors
    "bortezomib": "Proteasome inhibitors",
    "carfilzomib": "Proteasome inhibitors",
    "ixazomib": "Proteasome inhibitors",

    # Steroids
    "dexamethasone": "Steroids",
    "prednisone": "Steroids",
    "prednisolone": "Steroids",
    "methylprednisolone": "Steroids",

    # Monoclonal antibodies / antibody-based myeloma therapies
    "daratumumab": "Monoclonal antibodies",
    "elotuzumab": "Monoclonal antibodies",
    "isatuximab": "Monoclonal antibodies",
    "belantamab mafodotin": "Monoclonal antibodies",
    "elranatamab": "Monoclonal antibodies",
    "talquetamab": "Monoclonal antibodies",
    "teclistamab": "Monoclonal antibodies",

    # Alkylating agents
    "bendamustine": "Alkylating agents",
    "busulfan": "Alkylating agents",
    "carmustine": "Alkylating agents",
    "cyclophosphamide": "Alkylating agents",
    "melphalan": "Alkylating agents",
    "melphalan flufenamide": "Alkylating agents",

    # Other chemotherapy
    "carboplatin": "Other chemotherapy",
    "cisplatin": "Other chemotherapy",
    "cytarabine": "Other chemotherapy",
    "docetaxel": "Other chemotherapy",
    "doxorubicin": "Other chemotherapy",
    "epirubicin": "Other chemotherapy",
    "etoposide": "Other chemotherapy",
    "fludarabine": "Other chemotherapy",
    "ifosfamide": "Other chemotherapy",
    "pegylated liposomal doxorubicin": "Other chemotherapy",
    "vincristine": "Other chemotherapy",

    # Targeted/special categories that Fiona used
    "venetoclax": "BCL-2 inhibitor",
    "panobinostat": "HDAC inhibitor",
    "cart": "CAR-T",
    "investigational - cart": "CAR-T",

    # Transplant
    "autologous sct": "Transplant",
    "autologous sct1": "Transplant",
    "autologous sct2": "Transplant",
    "stem_cell_boost": "Transplant",

    # Categories not present in Fiona's list may still be useful for review.
    # These will usually produce a pattern marked as Not in Fiona's category list.
    "selinexor": "XPO1 inhibitor",
    "nivolumab": "Checkpoint inhibitors",
    "pembrolizumab": "Checkpoint inhibitors",
    "clarithromycin": "Other non-chemotherapy",
    "investigational - regimen": "Investigational regimen",
}


def normalize_text(value: object) -> str:
    """Normalize strings for matching while preserving readable family labels elsewhere."""
    if pd.isna(value):
        return ""
    value = str(value).strip()
    value = re.sub(r"\s+", " ", value)
    return value


def clean_drug_name(value: object) -> str:
    """
    Clean a parsed drug token before mapping it to a drug family.

    Important for malformed / split Excel values where a token may still contain
    leftover square brackets, for example: ']dexamethasone' or
    '[belantamab mafodotin'.
    """
    drug = normalize_text(value).lower()

    # Remove square brackets anywhere in the token, not only at the ends.
    drug = drug.replace("[", "").replace("]", "")

    # Remove common wrapping punctuation / quotes after bracket cleanup.
    drug = drug.strip().strip("'\".,;:")

    # Normalize whitespace again after cleanup.
    drug = re.sub(r"\s+", " ", drug).strip()
    return drug


def parse_lot_drugs(regimen: object) -> list[str]:
    """
    Parse COTA line_of_therapy_name values such as:
      [bortezomib, dexamethasone], [lenalidomide]

    Returns unique lowercase drug names in first-seen order.
    """
    text = normalize_text(regimen).lower()
    if not text:
        return []

    drugs: list[str] = []
    # COTA regimens appear as bracketed lists. Extract each bracket first.
    bracket_chunks = re.findall(r"\[([^\]]+)\]", text)
    if bracket_chunks:
        chunks = bracket_chunks
    else:
        # Fallback for unexpected formatting.
        chunks = [text]

    for chunk in chunks:
        for raw_drug in chunk.split(","):
            drug = clean_drug_name(raw_drug)
            if drug and drug not in drugs:
                drugs.append(drug)
    return drugs


def build_class_order(fiona_categories: Iterable[str]) -> list[str]:
    """
    Start with the known clinical order, then append any extra family labels
    found in Fiona's categories or in the drug mapping.
    """
    order = list(DEFAULT_CLASS_ORDER)
    seen = set(order)

    for category in fiona_categories:
        for part in normalize_text(category).split(" + "):
            if part and part not in seen:
                order.append(part)
                seen.add(part)

    for family in DRUG_TO_FAMILY.values():
        if family not in seen:
            order.append(family)
            seen.add(family)

    return order


def reconstruct_family(drugs: list[str], class_order: list[str]) -> tuple[str, list[str], list[str]]:
    """
    Map parsed drugs to family classes and return:
      family combination string, unknown drugs, mapped classes
    """
    mapped_classes = []
    unknown_drugs = []

    for drug in drugs:
        cleaned_drug = clean_drug_name(drug)

        # Try mapping with the cleaned value first. This prevents tokens with
        # leftover brackets from being falsely marked as unknown.
        family = DRUG_TO_FAMILY.get(cleaned_drug)

        if family:
            if family not in mapped_classes:
                mapped_classes.append(family)
        else:
            unknown_drugs.append(cleaned_drug)

    # Sort classes using Fiona/clinical order. Unknown family labels fall to the end.
    order_index = {family: i for i, family in enumerate(class_order)}
    mapped_classes = sorted(mapped_classes, key=lambda x: order_index.get(x, 10_000))

    family_combo = " + ".join(mapped_classes) if mapped_classes else "Unmapped regimen"
    return family_combo, unknown_drugs, mapped_classes


def reconstruct_family_combinations(cota_df: pd.DataFrame, fiona_line_df: pd.DataFrame) -> pd.DataFrame:
    """
    Main function to reconstruct family combinations from COTA data.
    """
    if "Family Combination" not in fiona_line_df.columns:
        raise ValueError("Fiona Line_by_Line sheet must contain a 'Family Combination' column.")
    if "line_of_therapy_name" not in cota_df.columns:
        raise ValueError("Raw COTA sheet must contain a 'line_of_therapy_name' column.")

    fiona_categories = sorted(
        {normalize_text(x) for x in fiona_line_df["Family Combination"].dropna() if normalize_text(x)}
    )
    fiona_category_set = set(fiona_categories)
    class_order = build_class_order(fiona_categories)

    parsed_drugs_col = []
    family_col = []
    fiona_match_col = []
    in_fiona_col = []
    unknown_drugs_col = []
    mapped_classes_col = []

    for regimen in cota_df["line_of_therapy_name"]:
        drugs = parse_lot_drugs(regimen)
        family_combo, unknown_drugs, mapped_classes = reconstruct_family(drugs, class_order)
        in_fiona = family_combo in fiona_category_set

        parsed_drugs_col.append(", ".join(drugs))
        mapped_classes_col.append(" + ".join(mapped_classes))
        unknown_drugs_col.append(", ".join(unknown_drugs))
        fiona_match_col.append(family_combo if in_fiona else "")
        in_fiona_col.append(in_fiona)

        if in_fiona:
            family_col.append(family_combo)
        else:
            family_col.append(f"{family_combo}")

    cota_df.insert(len(cota_df.columns), "parsed_drugs", parsed_drugs_col)
    cota_df.insert(len(cota_df.columns), "mapped_family_classes", mapped_classes_col)
    cota_df.insert(len(cota_df.columns), "reconstructed_family_combination", family_col)
    cota_df.insert(len(cota_df.columns), "in_provided_fiona_category_list", in_fiona_col)
    cota_df.insert(len(cota_df.columns), "fiona_category_match", fiona_match_col)
    cota_df.insert(len(cota_df.columns), "unknown_or_unmapped_drugs", unknown_drugs_col)

    return cota_df, fiona_categories, fiona_category_set

"""
SCRIPT 2: 
  - original file: all_misclass_in_one_file.py
  - purpose: repairs split/continuation rows from Script 1 output, 
    comparing COTA LoT numbers against reviewer-assigned numbers to find discrepancies
"""

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


def repair_split_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Repair split/continuation rows by merging fragments.
    """
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

    return pd.DataFrame(rows)


def add_transition_logic(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add transition alignment logic to identify misclassifications.
    """
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

    return df


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    # Input paths
    raw_path = ".." + os.sep + "data" + os.sep + "LoT Adjudication Datasets.xlsx"
    fiona_path = "COTA misclassification.xlsx"

    # Step 1: Read input files
    print("Reading input files...")
    cota = pd.read_excel(raw_path, sheet_name="Cota")
    fiona_line = pd.read_excel(fiona_path, sheet_name="Line_by_Line")

    # Step 2: Reconstruct family combinations (SCRIPT 1)
    print("Reconstructing family combinations...")
    cota, fiona_categories, fiona_category_set = reconstruct_family_combinations(cota, fiona_line)

    print(f"COTA rows processed: {len(cota)}")
    print(f"Fiona categories loaded from Line_by_Line: {len(fiona_categories)}")
    print(f"Rows matching Fiona categories: {sum(cota['in_provided_fiona_category_list'])}")
    print(f"Rows not matching Fiona categories: {len(cota) - sum(cota['in_provided_fiona_category_list'])}")

    # Step 3: Repair split rows and add transition logic (SCRIPT 2)
    print("\nRepairing split rows and adding transition logic...")

    # Ensure columns are properly stripped
    cota.columns = cota.columns.str.strip()

    print(f"Rows before split-row repair: {len(cota)}")

    # Repair split rows
    cota = repair_split_rows(cota)
    print(f"Rows after split-row repair: {len(cota)}")

    # Add transition logic
    cota = add_transition_logic(cota)

    # Step 4: Extract output tables
    output_path = Path("Output/COTA_misclassified_rows_UPD_2.xlsx")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create summary tables from SCRIPT 1
    summary_family = (
        cota.groupby(["reconstructed_family_combination", "in_provided_fiona_category_list"], dropna=False)
        .size()
        .reset_index(name="raw_cota_total_rows")
        .sort_values(["in_provided_fiona_category_list", "raw_cota_total_rows"], ascending=[False, False])
    )

    fiona_counts = Counter(normalize_text(x) for x in fiona_line["Family Combination"].dropna())
    fiona_categories_df = pd.DataFrame(
        [{"provided_fiona_family_combination": cat, "fiona_line_by_line_count": fiona_counts[cat]} for cat in
         fiona_categories]
    )

    unknown_summary = (
        cota.loc[cota["unknown_or_unmapped_drugs"].astype(bool), "unknown_or_unmapped_drugs"]
        .str.split(", ")
        .explode()
        .value_counts()
        .rename_axis("unknown_or_unmapped_drug")
        .reset_index(name="row_count")
    )

    # Extract misclassified rows
    misclassified = cota[
        cota["transition_alignment_status"].isin([
            "candidate_misclassification_cota_over_split",
            "candidate_misclassification_cota_under_split",
        ])
    ].copy()

    needs_context = cota[
        cota["transition_alignment_status"].eq("needs_context_missing_or_unparseable")
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

    preferred_columns = [col for col in preferred_columns if col in cota.columns]

    misclassified = misclassified[
        preferred_columns + [col for col in misclassified.columns if col not in preferred_columns]
        ]

    needs_context = needs_context[
        preferred_columns + [col for col in needs_context.columns if col not in preferred_columns]
        ]

    summary_misclass = (
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

    split_repair_audit = cota[
        cota["merged_from_split_rows"].eq(True)
    ].copy()

    # Write all sheets to final output
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # SCRIPT 1 sheets
        cota.to_excel(writer, sheet_name="COTA_Reconstructed", index=False)
        summary_family.to_excel(writer, sheet_name="Pattern_Summary", index=False)
        fiona_categories_df.to_excel(writer, sheet_name="Fiona_Category_List", index=False)
        unknown_summary.to_excel(writer, sheet_name="Unmapped_Drugs", index=False)

        # SCRIPT 2 sheets
        misclassified.to_excel(writer, sheet_name="Misclassified_Rows", index=False)
        summary_misclass.to_excel(writer, sheet_name="Misclassification_Summary", index=False)
        patient_summary.to_excel(writer, sheet_name="Patient_Summary", index=False)
        needs_context.to_excel(writer, sheet_name="Needs_Context", index=False)
        split_repair_audit.to_excel(writer, sheet_name="Split_Repair_Audit", index=False)
        cota.to_excel(writer, sheet_name="All_COTA_With_Transitions", index=False)

    print(f"\nFinal output saved: {output_path}")
    print(f"Misclassified rows: {len(misclassified)}")
    print(f"Needs context rows: {len(needs_context)}")
    print(f"Split-repaired LoT rows: {len(split_repair_audit)}")
    print("\nTransition alignment status counts:")
    print(cota["transition_alignment_status"].value_counts(dropna=False))


if __name__ == "__main__":
    main()