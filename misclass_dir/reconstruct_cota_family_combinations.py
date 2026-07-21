"""
Reconstruct COTA family combinations using Fiona's provided family-combination list.

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

Run example:
  python reconstruct_cota_family_combinations.py \
    --raw "LoT Adjudication Datasets(3).xlsx" \
    --fiona "COTA misclassification.xlsx" \
    --output "COTA_family_reconstructed.xlsx"
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd


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

#------------------------------------------------------------------------------------------------------------

raw_path = "LoT Adjudication Datasets.xlsx"
fiona_path = "COTA misclassification.xlsx"
output_path = Path("misclass_dir/Output/COTA_family_reconstructed.xlsx") #added this
# change the output name once confirming the path is correct
#originally: output_path = Path("MisclassificationInBestCases/Output/COTA_family_reconstructed.xlsx")
output_path.parent.mkdir(parents=True, exist_ok=True)

cota = pd.read_excel(raw_path, sheet_name="Cota")
fiona_line = pd.read_excel(fiona_path, sheet_name="Line_by_Line")

if "Family Combination" not in fiona_line.columns:
    raise ValueError("Fiona Line_by_Line sheet must contain a 'Family Combination' column.")
if "line_of_therapy_name" not in cota.columns:
    raise ValueError("Raw COTA sheet must contain a 'line_of_therapy_name' column.")

fiona_categories = sorted(
    {normalize_text(x) for x in fiona_line["Family Combination"].dropna() if normalize_text(x)}
)
fiona_category_set = set(fiona_categories)
class_order = build_class_order(fiona_categories)

parsed_drugs_col = []
family_col = []
fiona_match_col = []
in_fiona_col = []
unknown_drugs_col = []
mapped_classes_col = []

for regimen in cota["line_of_therapy_name"]:
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

cota.insert(len(cota.columns), "parsed_drugs", parsed_drugs_col)
cota.insert(len(cota.columns), "mapped_family_classes", mapped_classes_col)
cota.insert(len(cota.columns), "reconstructed_family_combination", family_col)
cota.insert(len(cota.columns), "in_provided_fiona_category_list", in_fiona_col)
cota.insert(len(cota.columns), "fiona_category_match", fiona_match_col)
cota.insert(len(cota.columns), "unknown_or_unmapped_drugs", unknown_drugs_col)

# Summary by reconstructed family output string.
summary = (
    cota.groupby(["reconstructed_family_combination", "in_provided_fiona_category_list"], dropna=False)
    .size()
    .reset_index(name="raw_cota_total_rows")
    .sort_values(["in_provided_fiona_category_list", "raw_cota_total_rows"], ascending=[False, False])
)

# Fiona categories and their line-by-line counts, based on Line_by_Line only.
fiona_counts = Counter(normalize_text(x) for x in fiona_line["Family Combination"].dropna())
fiona_categories_df = pd.DataFrame(
    [{"provided_fiona_family_combination": cat, "fiona_line_by_line_count": fiona_counts[cat]} for cat in fiona_categories]
)

unknown_summary = (
    cota.loc[cota["unknown_or_unmapped_drugs"].astype(bool), "unknown_or_unmapped_drugs"]
    .str.split(", ")
    .explode()
    .value_counts()
    .rename_axis("unknown_or_unmapped_drug")
    .reset_index(name="row_count")
)

with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
    cota.to_excel(writer, sheet_name="COTA_Reconstructed", index=False)
    summary.to_excel(writer, sheet_name="Pattern_Summary", index=False)
    fiona_categories_df.to_excel(writer, sheet_name="Fiona_Category_List", index=False)
    unknown_summary.to_excel(writer, sheet_name="Unmapped_Drugs", index=False)

print(f"Saved: {output_path}")
print(f"COTA rows processed: {len(cota)}")
print(f"Fiona categories loaded from Line_by_Line: {len(fiona_categories)}")
print(f"Rows matching Fiona categories: {sum(in_fiona_col)}")
print(f"Rows not matching Fiona categories: {len(cota) - sum(in_fiona_col)}")