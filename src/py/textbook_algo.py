"""
LOT (Lines of Therapy) Algorithm for Multiple Myeloma 
========================================================================================
Based on Rajkumar guidelines

Outputs a comparison table to stdout and saves results to lot_results.xlsx
"""

import pandas as pd

# ---------------------------------------------------------------------------
# Drug synonym table
# ---------------------------------------------------------------------------
DRUG_SYNONYMS = {
    # Daratumumab variants
    "Daratumumab/Hyaluronidase-Fihj": "Daratumumab",

    # Isatuximab variants
    "Isatuximab-Irfc": "Isatuximab",

    # Common formatting / minor variants (defensive normalization)
    "Bortezomib SC": "Bortezomib",
    "Bortezomib IV": "Bortezomib",

    "Carfilzomib IV": "Carfilzomib",

    "Lenalidomide (Revlimid)": "Lenalidomide",
    "Pomalidomide (Pomalyst)": "Pomalidomide",

    # Dexamethasone sometimes appears abbreviated or combined
    "Dex": "Dexamethasone",

    # Cyclophosphamide variants
    "Cyclophosphamide IV": "Cyclophosphamide",
    "Cyclophosphamide PO": "Cyclophosphamide",

    # Melphalan variants
    "Melphalan IV": "Melphalan",
    "Melphalan PO": "Melphalan",

    # Prednisone variants
    "Prednisone PO": "Prednisone",
}

# Tokens that represent non-treatment rows and should be excluded
NON_TREATMENT_TOKENS = {"Transplant", "Line Zero"}


def normalize_drug(drug: str) -> str:
    """Return the canonical name for a drug, resolving known synonyms."""
    return DRUG_SYNONYMS.get(drug.strip(), drug.strip())


def parse_regimen(regimen_str: str) -> frozenset[str]:
    """Split a comma-separated regimen string into a frozenset of normalized drug names."""
    return frozenset(normalize_drug(d) for d in regimen_str.split(",") if d.strip())


# ---------------------------------------------------------------------------
# Fixed LOT algorithm
# ---------------------------------------------------------------------------

def lot_algorithm_fixed(trt_list: list[str]) -> int:
    """
    Count lines of therapy for a single patient.

    Parameters
    ----------
    trt_list : Ordered list of treatment regimen strings (comma-separated drugs each).

    Returns
    -------
    Number of lines of therapy.
    """
    # Remove non-treatment rows (Transplant, Line Zero)
    trts = [t for t in trt_list if t not in NON_TREATMENT_TOKENS]

    if not trts:
        return 0

    lines = 1
    prev = parse_regimen(trts[0])   # drugs in the immediately preceding regimen

    for i in range(1, len(trts)):
        curr = parse_regimen(trts[i])
        overlap = curr & prev          # drugs shared with the immediate previous regimen
        added   = curr - prev          # drugs new relative to the immediate previous regimen

        if overlap:
            # Rule 1: there is still an overlapping drug, so the previous regimen has NOT been fully discontinued.
            if added:
                # Rule 2: one or more new drugs were added (unplanned addition/substitution)
                # -> always a new line of therapy.
                # (FIX: removed the original trt_prev subset check here)
                lines += 1
            # in the other case, only drugs were removed (partial discontinuation) -> not a new line
        else:
            # No drug in common with the immediate previous regimen.
            # The previous regimen has been fully discontinued and a new one started.
            # -> new line of therapy (Rajkumar rule 1).
            lines += 1

        prev = curr

    return lines


# ---------------------------------------------------------------------------
# Original algorithm (Python translation of the R code)
# ---------------------------------------------------------------------------

def lot_algorithm_original(trt_list: list[str]) -> int:
    """
    Translation of the original R algorithm (unfixed) for comparison.
    """
    trts = [t for t in trt_list if t != "Transplant"]
    if not trts:
        return 1

    lines = 1
    trt_prev = set(trts[0].split(","))

    for i in range(1, len(trts)):
        trt2 = set(trts[i].split(","))
        trt1 = set(trts[i - 1].split(","))

        if trt2 & trt1:                       
            if trt2 - trt1:                   
                if not trt2.issubset(trt_prev):   
                    lines += 1
                    trt_prev = trt2.copy()
        else:                                 
            lines += 1
            trt_prev = trt2.copy()

    return lines


# ---------------------------------------------------------------------------
# Data loading (mirrors the R cleaning logic)
# ---------------------------------------------------------------------------

def load_patient_data(filepath: str, sheet: int = 0) -> pd.DataFrame:
    """
    Load and reshape the Flatiron adjudication spreadsheet into one row per patient.

    Returns a DataFrame with columns:
        f, trt (list), rev1_lot, rev2_lot, flatiron_lot
    """
    df_raw = pd.read_excel(filepath, sheet_name=sheet, header=None, skiprows=3)

    cols = [
        "f", "diagnosis_date", "death_flag", "death_date", "treatment",
        "rev1_name", "rev1_lot", "rev2_name", "rev2_lot", "flatiron_lot",
        "startdate", "enddate", "response", "assessment_date", "specimen_type",
        "residual", "col17", "col18",
    ]
    df_raw.columns = cols

    df = df_raw[["f", "treatment", "rev1_lot", "rev2_lot", "flatiron_lot"]].copy()

    patients, current_f, current_trts = [], None, []
    current_rev1 = current_rev2 = current_flat = None

    for _, row in df.iterrows():
        f   = row["f"]
        trt = row["treatment"]

        if pd.isna(f):
            # Continuation row: append to the last treatment of the current patient
            if not pd.isna(trt) and current_trts:
                current_trts[-1] = str(current_trts[-1]) + str(trt)
        elif f == current_f:
            current_rev1 = row["rev1_lot"]
            current_rev2 = row["rev2_lot"]
            current_flat = row["flatiron_lot"]
            if not pd.isna(trt):
                current_trts.append(str(trt))
        else:
            if current_f is not None:
                patients.append({
                    "f": current_f,
                    "trt": current_trts,
                    "rev1_lot": current_rev1,
                    "rev2_lot": current_rev2,
                    "flatiron_lot": current_flat,
                })
            current_f    = f
            current_rev1 = row["rev1_lot"]
            current_rev2 = row["rev2_lot"]
            current_flat = row["flatiron_lot"]
            current_trts = [] if pd.isna(trt) else [str(trt)]

    if current_f is not None:
        patients.append({
            "f": current_f,
            "trt": current_trts,
            "rev1_lot": current_rev1,
            "rev2_lot": current_rev2,
            "flatiron_lot": current_flat,
        })

    df_patient = pd.DataFrame(patients)
    df_patient["rev1_lot"]     = pd.to_numeric(df_patient["rev1_lot"],     errors="coerce")
    df_patient["rev2_lot"]     = pd.to_numeric(df_patient["rev2_lot"],     errors="coerce")
    df_patient["flatiron_lot"] = pd.to_numeric(df_patient["flatiron_lot"], errors="coerce")
    return df_patient


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    filepath = "LoT Adjudication Datasets.xlsx"

    print(f"Loading data from '{filepath}' ...\n")
    df = load_patient_data(filepath)
    n  = len(df)
    print(f"  {n} patients loaded.\n")

    # Run both algorithms
    df["lot_original"] = df["trt"].apply(lot_algorithm_original)
    df["lot_fixed"]    = df["trt"].apply(lot_algorithm_fixed)

    # Deltas vs reviewer 1
    df["orig_err"]  = df["lot_original"] - df["rev1_lot"]
    df["fixed_err"] = df["lot_fixed"]    - df["rev1_lot"]

    # -----------------------------------------------------------------------
    # Accuracy summary
    # -----------------------------------------------------------------------
    def accuracy_row(errors, label):
        correct = (errors == 0).sum()
        over    = (errors  > 0).sum()
        under   = (errors  < 0).sum()
        print(f"  {label}")
        print(f"    Correct : {correct}/{n}  ({correct/n:.1%})")
        print(f"    Over    : {over}/{n}     ({over/n:.1%})")
        print(f"    Under   : {under}/{n}    ({under/n:.1%})")
        print()

    print("=" * 60)
    print("ACCURACY vs Reviewer 1 LOT")
    print("=" * 60)
    accuracy_row(df["orig_err"],  "Original algorithm")
    accuracy_row(df["fixed_err"], "Fixed algorithm   ")

    # -----------------------------------------------------------------------
    # Cases where the fix changes the result
    # -----------------------------------------------------------------------
    changed = df[df["lot_original"] != df["lot_fixed"]].copy()
    print("=" * 60)
    print(f"CASES CHANGED BY THE FIX  ({len(changed)} patients)")
    print("=" * 60)
    for _, row in changed.iterrows():
        orig_correct  = row["orig_err"]  == 0
        fixed_correct = row["fixed_err"] == 0
        verdict = (
            "improved" if not orig_correct and fixed_correct else
            "worsened" if orig_correct and not fixed_correct else
            "both wrong (different error)"
        )
        print(
            f"  rev1={int(row['rev1_lot'])}  "
            f"orig={int(row['lot_original'])}  "
            f"fixed={int(row['lot_fixed'])}  "
            f"-> {verdict}"
        )
        print(f"    Regimen sequence: {row['trt']}")
        print()

    # -----------------------------------------------------------------------
    # Remaining misclassifications after fix
    # -----------------------------------------------------------------------
    still_wrong = df[df["fixed_err"] != 0]
    print("=" * 60)
    print(f"STILL MISCLASSIFIED AFTER FIX  ({len(still_wrong)} patients)")
    print("=" * 60)
    print("These require AI-medium or human review (drug rechallenge context,")
    print("maintenance vs active therapy, solo-drug intent, etc.)\n")
    for _, row in still_wrong.iterrows():
        direction = "over" if row["fixed_err"] > 0 else "under"
        print(
            f"  rev1={int(row['rev1_lot'])}  "
            f"fixed={int(row['lot_fixed'])}  "
            f"({direction} by {abs(int(row['fixed_err']))})"
        )
        print(f"    {row['trt']}")
        print()

    # -----------------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------------
    output_cols = ["f", "rev1_lot", "rev2_lot", "flatiron_lot",
                   "lot_original", "orig_err", "lot_fixed", "fixed_err", "trt"]
    out = df[output_cols].copy()
    out["trt"] = out["trt"].apply(lambda x: " | ".join(x))
    out.to_excel("lot_results_flatiron.xlsx", index=False)
    print("Results saved to lot_results.xlsx")


if __name__ == "__main__":
    main()