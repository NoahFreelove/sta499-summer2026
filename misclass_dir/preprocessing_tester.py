"""
Updated and Testing Code from COTA_preprocessing_newrules
- See if I can remove redundant code, and check the accuracy of the newly added code.
- ask if the conditioning drugs should be listed individually or as
conditioning regimens instead (are they all excluded from the active drug set?)
conditioning_regimens = [
        {"melphalan"},  # High-dose melphalan alone
        {"bortezomib", "melphalan"},  # Bor-HDM
        {"busulfan", "melphalan"},  # BUMEL
        {"busulfan", "melphalan", "bortezomib"},  # BuMelVel
        # ... keep your existing BEAM definition if needed
        {"carmustine", "cytarabine", "etoposide", "melphalan"},  # BEAM
    ]

Pre-processing and reconstruction of data for LoT Algorithm
- Combines 2 pre-processing scripts from Ana for simplified usage.
- Includes pre-processing rules from the updated document, such as .
    1. identify_active_drugs() - Identifies active drugs by excluding:
        - Steroids (dexamethasone, prednisone, prednisolone, methylprednisolone)
        - Conditioning drugs (melphalan, busulfan, carmustine???)
    2. identify_steroid_only_segment() - Checks if a segment contains only steroid drugs
    3. is_car_t_segment() - Identifies CAR-T segments
    4. is_asct_segment() - Identifies ASCT segments
    **is_car_t_segment() and is_asct_segment() replaced with check_segment_type
    5. calculate_gap_days() - Calculates gap in days between two dates -- modified
    6. assign_asct_index() - Assigns chronological ASCT index per patient:
        - First ASCT → asct_index = 1
        - Second ASCT → asct_index = 2
        - Independent of LOT numbering
    7. removed: identify_steroid_only_absorption() - Identifies steroid-only first segments that should be absorbed:
        - Only steroids and gap to first non-steroid segment ≤ 7 days
    **apply_all_preprocessing() and calculate_gaps_and_mark_absorption() added
    **create_streamlined_output() added
    8. removed? create_preprocessed_dataframe() - Main function that:
        - Creates drugs_active list
        - Creates tail_drugs_active (active drugs from most recently absorbed sub-segment)
        - Identifies steroid-only, CAR-T, and ASCT segments
        - Assigns ASCT indices
        - Calculates gaps between segments
        - Marks steroid-only segments for absorption

Updates to be incorporated:
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

# ============================================================================
# SCRIPT 1: Reconstruct COTA family combinations
# original file: reconstruct_cota_family_combinations.py
# purpose: reconstruct COTA family combinations using Fiona's provided family-combination list
# ============================================================================

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


# ============================================================================
# SCRIPT 2: Split-row repair and transition logic
# original file: all_misclass_in_one_file.py
# purpose: repairs split/continuation rows from Script 1 output,
# comparing COTA LoT numbers against reviewer-assigned numbers to find discrepancies
# ============================================================================

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


# ============================================================================
# NEW PRE-PROCESSING FUNCTIONS BASED ON RULE ENGINE SPECIFICATION
# ============================================================================

def identify_active_drugs(drugs_list: list) -> list:
    """
    Identify active drugs (drugs_active) - excludes steroids and SCT/conditioning drugs.

    According to the spec:
    - drugs_active: Active drug set --- excludes steroids (DEX, PRED, etc.)
      and SCT/conditioning drugs (Melphalan, Busulfan)
    """
    steroids = {"dexamethasone", "prednisone", "prednisolone", "methylprednisolone"}
    conditioning_drugs = {"melphalan", "busulfan"}

    active_drugs = []
    for drug in drugs_list:
        drug_lower = str(drug).lower().strip()
        if drug_lower not in steroids and drug_lower not in conditioning_drugs:
            active_drugs.append(drug_lower)

    return active_drugs


def identify_steroid_only_segment(drugs_list: list) -> bool:
    """
    Check if a segment contains only steroid drugs.
    """
    steroids = {"dexamethasone", "prednisone", "prednisolone", "methylprednisolone"}

    if not drugs_list:
        return False

    all_steroids = all(str(drug).lower().strip() in steroids for drug in drugs_list)
    return all_steroids


def check_segment_type(regimen_text: str, segment_type: str) -> bool:
    """Check if segment matches a specific type. If a segment contains CAR-T therapy,
    or ASCT (Autologous Stem Cell Transplant)."""
    if pd.isna(regimen_text):
        return False

    text_lower = str(regimen_text).lower()

    indicators = {
        'car_t': ["cart", "car-t", "chimeric antigen receptor"],
        'asct': ["autologous", "sct", "transplant", "stem cell"]
    }

    if segment_type not in indicators:
        return False

    # For ASCT, exclude CAR-T
    if segment_type == 'asct' and any(t in text_lower for t in ['cart', 'car-t']):
        return False

    return any(indicator in text_lower for indicator in indicators[segment_type])


def calculate_gap_days(date1, date2) -> float:
    """Calculate gap in days between two dates. Returns NaN if dates are missing."""
    if pd.isna(date1) or pd.isna(date2):
        return np.nan

    date1 = pd.to_datetime(date1, errors='coerce')
    date2 = pd.to_datetime(date2, errors='coerce')

    if pd.isna(date1) or pd.isna(date2):
        return np.nan

    return float((date2 - date1).days)


def assign_asct_index(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign ASCT index to each ASCT event per patient.

    According to the spec:
        - Scan all segments chronologically and assign asct_index to each ASCT event
        - First ASCT encountered → asct_index = 1
        - Second ASCT → asct_index = 2, and so on
        - This index is independent of LOT numbering
    """
    # Sort directly
    df = df.sort_values(
        ['cpid', 'date_start_line_of_therapy'] if 'date_start_line_of_therapy' in df.columns else ['cpid']
    )

    df["asct_index"] = np.nan

    for patient in df["cpid"].unique():
        asct_counter = 0
        patient_mask = df["cpid"] == patient
        for idx in df[patient_mask].index:
            # Use the pre-computed is_asct column
            if df.loc[idx, "is_asct"]:
                asct_counter += 1
                df.loc[idx, "asct_index"] = asct_counter

    return df


def apply_all_preprocessing(df: pd.DataFrame) -> pd.DataFrame:
    """Apply all preprocessing steps in a single pipeline."""
    result = df.copy()

    # Parse drugs once and reuse
    result['parsed_drugs_list'] = result['line_of_therapy_name'].apply(parse_lot_drugs)

    # Create all derived columns in one pass
    result['drugs_active'] = result['parsed_drugs_list'].apply(identify_active_drugs)
    result['drugs_active_count'] = result['drugs_active'].apply(len)
    result['is_steroid_only'] = result['parsed_drugs_list'].apply(identify_steroid_only_segment)
    result['is_car_t'] = result['line_of_therapy_name'].apply(lambda x: check_segment_type(x, 'car_t'))
    result['is_asct'] = result['line_of_therapy_name'].apply(lambda x: check_segment_type(x, 'asct'))

    # Assign ASCT index in one pass
    result = assign_asct_index(result)

    # Calculate gaps and mark absorption in one pass
    result = calculate_gaps_and_mark_absorption(result)

    # Set tail drugs active (copy for now)
    result['tail_drugs_active'] = result['drugs_active']
    result['preprocessing_merged'] = False

    return result


def calculate_gaps_and_mark_absorption(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate gaps and mark steroid-only segments for absorption in one pass.


    Combined identify_steroid_only_absorption():
    Identify steroid-only first segments that should be absorbed.
    According to the spec:
    If the patient's first treatment segment contains only steroids and the
    gap to the first non-steroid segment is ≤ 7 days, it is absorbed into
    that first active segment (same LOT start).
    """

    df = df.sort_values(
        ['cpid', 'date_start_line_of_therapy'] if 'date_start_line_of_therapy' in df.columns else ['cpid'])

    df['segment_gap_to_previous'] = np.nan
    df['is_first_patient_segment'] = False
    df['absorb_steroid_only'] = False

    for patient in df['cpid'].unique():
        patient_mask = df['cpid'] == patient
        patient_indices = df[patient_mask].index

        for i, idx in enumerate(patient_indices):
            if i == 0:
                df.loc[idx, 'is_first_patient_segment'] = True
            else:
                prev_idx = patient_indices[i - 1]
                gap = calculate_gap_days(
                    df.loc[prev_idx, 'date_end_line_of_therapy'],
                    df.loc[idx, 'date_start_line_of_therapy']
                )
                df.loc[idx, 'segment_gap_to_previous'] = gap

                # Check absorption for first segment
                if i == 1 and df.loc[prev_idx, 'is_steroid_only'] and not df.loc[idx, 'is_steroid_only']:
                    if not pd.isna(gap) and gap <= 7:
                        df.loc[prev_idx, 'absorb_steroid_only'] = True

    return df


def create_streamlined_output(df: pd.DataFrame, output_path: Path):
    """Create fewer, more informative output sheets."""
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:

        # 1. Main data with all columns
        df.to_excel(writer, sheet_name='All_Data', index=False)

        # 2. Misclassified rows (combine misclassified and needs_context)
        issues_df = df[df['transition_alignment_status'].str.startswith('candidate_') |
                       df['transition_alignment_status'].eq('needs_context_missing_or_unparseable')]
        issues_df.to_excel(writer, sheet_name='Issues_Review', index=False)

        # 3. Summary (combine multiple summaries)
        summary_data = {
            'Category': [],
            'Count': []
        }

        # Pattern summary
        pattern_summary = df.groupby('reconstructed_family_combination').size()
        for pattern, count in pattern_summary.items():
            summary_data['Category'].append(f'Pattern: {pattern}')
            summary_data['Count'].append(count)

        # Transition status
        for status, count in df['transition_alignment_status'].value_counts().items():
            summary_data['Category'].append(f'Transition: {status}')
            summary_data['Count'].append(count)

        # ASCT summary
        asct_count = df[df['is_asct']].shape[0]
        summary_data['Category'].append('ASCT Events')
        summary_data['Count'].append(asct_count)

        # Steroid-only absorption
        steroid_abs = df[df['absorb_steroid_only']].shape[0]
        summary_data['Category'].append('Steroid-only Absorbed')
        summary_data['Count'].append(steroid_abs)

        pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)

        # 4. Unmapped drugs (still useful for review)
        unknown_drugs = df[df['unknown_or_unmapped_drugs'].astype(bool)]['unknown_or_unmapped_drugs']
        if not unknown_drugs.empty:
            drug_counts = unknown_drugs.str.split(', ').explode().value_counts()
            drug_counts.to_frame('Count').to_excel(writer, sheet_name='Unmapped_Drugs')


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
    # Load data
    raw_path = ".." + os.sep + "data" + os.sep + "LoT Adjudication Datasets.xlsx"
    fiona_path = "COTA misclassification.xlsx"

    print("Reading input files...")
    cota = pd.read_excel(raw_path, sheet_name="Cota")
    fiona_line = pd.read_excel(fiona_path, sheet_name="Line_by_Line")

    # Step 1: Reconstruct family combinations
    print("Reconstructing family combinations...")
    cota, fiona_categories, fiona_category_set = reconstruct_family_combinations(cota, fiona_line)

    # Step 2: Repair split rows
    print("\nRepairing split rows...")
    cota.columns = cota.columns.str.strip()
    print(f"Rows before split-row repair: {len(cota)}")
    cota = repair_split_rows(cota)
    print(f"Rows after split-row repair: {len(cota)}")

    # Step 3: Apply combined preprocessing
    print("\nApplying combined preprocessing...")
    cota = apply_all_preprocessing(cota)

    # Step 4: Add transition logic
    print("\nAdding transition logic...")
    cota = add_transition_logic(cota)

    # Step 5: Create streamlined output
    output_path = Path("Output/COTA_cleaned.xlsx")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    create_streamlined_output(cota, output_path)

    print(f"\nFinal output saved: {output_path}")
    print("\nSummary:")
    print(f"Total rows: {len(cota)}")
    print(f"Transition status counts:")
    print(cota["transition_alignment_status"].value_counts(dropna=False))


if __name__ == "__main__":
    main()


"""
conditioning_regimens = [
    {"melphalan"},  # High-dose melphalan alone
    {"bortezomib", "melphalan"},  # Bor-HDM
    {"busulfan", "melphalan"},  # BUMEL
    {"busulfan", "melphalan", "bortezomib"},  # BuMelVel
    # ... keep your existing BEAM definition if needed
    {"carmustine", "cytarabine", "etoposide", "melphalan"},  # BEAM
]
"""
