"""
LOT (Lines of Therapy) Algorithm for Multiple Myeloma - COTA Dataset
=====================================================================

COTA-specific data issues 
-------------------------
  ISSUE 1 - Duplicate consecutive lines (data artefact).
            The same regimen sometimes appears in back-to-back COTA rows for
            the same patient. Per Rajkumar, restarting the same regimen without
            any intervening regimen = same line. These duplicates inflate COTA's
            own line count above both the algorithm and the reviewer.
            Fix: un-duplicate consecutive identical drug sets before applying rules.

  ISSUE 2 - SCT variants embedded in line entries.
            COTA bundles autologous SCT into the same line entry as the drugs
            (e.g. "[bortezomib, dexamethasone], [autologous SCT]"). Multiple
            variant spellings exist in the data.
            Fix: strip all SCT tokens before drug comparison (same logic as
            removing "Transplant" in the Flatiron script).

  ISSUE 3 - Partial drug removal split as a new line by COTA.
            COTA's pre-adjudication sometimes records a drug being dropped as a
            separate line entry (e.g. BDL -> BDL -> DL as three rows). Per Rajkumar
            rule 1, partial discontinuation is NOT a new line. The algorithm
            correctly merges these - but this means the algorithm will under-count
            relative to reviewers who follow COTA's structure. These are flagged
            in the output as AI-medium cases for human review.

  ISSUE 4 - Investigational regimens and CAR-T therapies.
            Some entries have no parseable drug names (e.g. "investigational -
            regimen") or represent CAR-T cell therapies bundled with lymphodepletion
            regimens. These cannot be adjudicated by rule-based logic.
            Fix: detect and flag these as requiring human review.

Reviewer columns
----------------
COTA uses two reviewers identified as Alpesh (col: Alpesh 1 LoT) and Alberto
(col: Alberto LOT). Both store a per-line LOT number; the patient's total LOT
is the maximum value across all their lines. This script compares the algorithm
against Alpesh as the primary reviewer (same convention as the Flatiron script
using Reviewer 1).

Outputs a comparison table to stdout and saves results to lot_results_cota.xlsx.
"""

import re
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# SCT token variants found in the COTA data - all excluded before comparison
SCT_TOKENS: set[str] = {
    "autologous sct",
    "autologous sct1",
    "autologous sct2",
    "autologoussct",
    "autologoussct2",
    "allogeneic sct",
}

# Drug synonym table - extend as new formulation variants are discovered
DRUG_SYNONYMS = {
    # Names are already lowercase and consistent but dict is here for future additions.
}

# Strings that indicate a line cannot be adjudicated by rule-based logic
NON_PARSEABLE_MARKERS = [
    "investigational",
]

# CAR-T marker - entries containing this are flagged for human review
CART_MARKER = "cart"


# ---------------------------------------------------------------------------
# Drug name helpers
# ---------------------------------------------------------------------------

def normalize_drug(drug: str) -> str:
    """Lowercase, strip whitespace, and resolve known synonyms."""
    d = drug.strip().lower()
    return DRUG_SYNONYMS.get(d, d)


def parse_cota_line(line_str: str) -> frozenset[str]:
    """
    Extract all drug names from a COTA line string.

    COTA format:  "[drug1, drug2], [drug3]"  (one or more bracketed groups)
    Continuation rows may be partial:  "drug3]"  or  ", drug4]"

    Returns a frozenset of normalized drug names, with SCT tokens removed.
    Returns an empty frozenset if the line is entirely SCT or empty.
    """
    if not line_str or not line_str.strip():
        return frozenset()

    # Extract content from all [...] groups
    groups = re.findall(r"\[([^\]]+)\]", line_str)

    if not groups:
        # No brackets found - treat the whole string as a single drug group
        # (handles partial continuation rows)
        cleaned = re.sub(r"[\[\]]", "", line_str).strip()
        groups = [cleaned] if cleaned else []

    drugs: set[str] = set()
    for group in groups:
        for raw in group.split(","):
            d = normalize_drug(raw)
            if d and d not in SCT_TOKENS:
                drugs.add(d)

    return frozenset(drugs)


def is_non_parseable(drug_set: frozenset[str]) -> bool:
    """Return True if any drug name contains a non-parseable marker."""
    return any(
        marker in drug
        for drug in drug_set
        for marker in NON_PARSEABLE_MARKERS
    )


def contains_cart(drug_set: frozenset[str]) -> bool:
    """Return True if the line contains a CAR-T therapy entry."""
    return any(CART_MARKER in drug for drug in drug_set)


# ---------------------------------------------------------------------------
# Data loading and cleaning
# ---------------------------------------------------------------------------

def load_patient_data(filepath: str, sheet: str = "Cota") -> pd.DataFrame:
    """
    Load and reshape the COTA sheet into one row per patient.

    Handles:
      - Header row in row 0 (no skip needed, unlike Flatiron)
      - NA rows where cpid is blank but line_of_therapy_name continues
        the previous row's treatment string
      - Multiple reviewer LOT columns (Alpesh, Alberto)

    Returns a DataFrame with columns:
        cpid, trt_sequence (list of frozensets), rev_alpesh, rev_alberto, cota_lot
    """
    df_raw = pd.read_excel(filepath, sheet_name=sheet, header=0)

    # Forward-fill cpid and line number into continuation rows
    current_cpid = None
    current_line = None
    for i, row in df_raw.iterrows():
        if pd.notna(row["cpid"]):
            current_cpid = row["cpid"]
            current_line = row["line_of_therapy_c"]
        else:
            df_raw.at[i, "cpid"] = current_cpid
            df_raw.at[i, "line_of_therapy_c"] = current_line

    df_raw["line_of_therapy_name"] = df_raw["line_of_therapy_name"].fillna("")
    df_raw["line_of_therapy_c"]    = pd.to_numeric(df_raw["line_of_therapy_c"],  errors="coerce")
    df_raw["Alpesh 1 LoT"]         = pd.to_numeric(df_raw["Alpesh 1 LoT"],        errors="coerce")
    df_raw["Alberto LOT "]         = pd.to_numeric(df_raw["Alberto LOT "],         errors="coerce")

    # Concatenate treatment name fragments that span multiple raw rows
    # (same cpid + same line number -> join the strings)
    grouped = (
        df_raw
        .groupby(["cpid", "line_of_therapy_c"], sort=False)["line_of_therapy_name"]
        .apply(lambda parts: "".join(str(p) for p in parts))
        .reset_index()
        .rename(columns={"line_of_therapy_name": "trt_full"})
    )

    # Reviewer totals per patient = max line number assigned by each reviewer
    reviewer = (
        df_raw
        .groupby("cpid", sort=False)
        .agg(
            rev_alpesh=("Alpesh 1 LoT",      lambda s: pd.to_numeric(s, errors="coerce").max()),
            rev_alberto=("Alberto LOT ",     lambda s: pd.to_numeric(s, errors="coerce").max()),
            cota_lot=("line_of_therapy_c",   lambda s: pd.to_numeric(s, errors="coerce").max()),
        )
        .reset_index()
    )

    # Parse drug sets and sort lines within each patient
    grouped["drugs"] = grouped["trt_full"].apply(parse_cota_line)
    grouped = grouped.sort_values(["cpid", "line_of_therapy_c"])

    # Build one row per patient with ordered treatment sequence
    patients = []
    for cpid, grp in grouped.groupby("cpid", sort=False):
        patients.append({
            "cpid":         cpid,
            "trt_sequence": grp["drugs"].tolist(),
            "line_nums":    grp["line_of_therapy_c"].tolist(),
            "trt_strings":  grp["trt_full"].tolist(),
        })

    df_patient = pd.DataFrame(patients).merge(reviewer, on="cpid")
    return df_patient


# ---------------------------------------------------------------------------
# LOT algorithm (COTA version)
# ---------------------------------------------------------------------------

def lot_algorithm_cota(trt_sequence: list[frozenset]) -> tuple[int, list[str]]:
    """
    Count lines of therapy for a single COTA patient and return flags.

    Parameters
    ----------
    trt_sequence : Ordered list of drug sets, one per COTA line entry (SCT tokens already removed).

    Returns
    -------
    lot : int
        Number of lines of therapy per Rajkumar rules.
    flags : list of str
        Per-patient flags for cases needing review:
          'duplicate_line'     - consecutive identical regimens (data artefact)
          'partial_removal'    - drug dropped without addition (merged by algo,
                                 may differ from reviewer who follows COTA structure)
          'non_parseable'      - contains investigational/unnamed regimen
          'cart_therapy'       - contains CAR-T entry
    """
    flags = []

    # --- Data quality flags (before deduplication) ---
    for i in range(1, len(trt_sequence)):
        if trt_sequence[i] == trt_sequence[i - 1] and trt_sequence[i]:
            flags.append("duplicate_line")
            break

    # --- Remove consecutive duplicate lines (data artefact fix) ---
    deduped = []
    for drugs in trt_sequence:
        if not deduped or drugs != deduped[-1]:
            deduped.append(drugs)

    # --- Remove empty sets (lines that were entirely SCT) ---
    cleaned = [d for d in deduped if d]

    if not cleaned:
        return 0, flags

    # --- Flag non-parseable and CAR-T entries ---
    for drugs in cleaned:
        if is_non_parseable(drugs):
            if "non_parseable" not in flags:
                flags.append("non_parseable")
        if contains_cart(drugs):
            if "cart_therapy" not in flags:
                flags.append("cart_therapy")

    # --- Apply Rajkumar rules ---
    lines = 1
    prev = cleaned[0]

    for i in range(1, len(cleaned)):
        curr = cleaned[i]
        overlap = curr & prev
        added   = curr - prev

        if overlap:
            if added:
                # Rule 2: unplanned addition/substitution -> new line
                lines += 1
            else:
                # Partial discontinuation only -> same line (Rajkumar rule 1)
                # Flag: COTA may have counted this as a separate line
                if "partial_removal" not in flags:
                    flags.append("partial_removal")
        else:
            # Full switch -> new line (Rajkumar rule 1)
            lines += 1

        prev = curr

    return lines, flags


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    filepath = "LoT Adjudication Datasets.xlsx"

    print(f"Loading COTA data from '{filepath}' ...\n")
    df = load_patient_data(filepath)
    n  = len(df)
    print(f"  {n} patients loaded.\n")

    # Run algorithm
    results           = df["trt_sequence"].apply(lot_algorithm_cota)
    df["lot_alg"]     = results.apply(lambda x: x[0])
    df["flags"]       = results.apply(lambda x: x[1])
    df["rev_alpesh"]  = df["rev_alpesh"].astype(int)
    df["cota_lot"]    = df["cota_lot"].astype(int)
    df["alg_err"]     = df["lot_alg"] - df["rev_alpesh"]
    df["cota_err"]    = df["cota_lot"] - df["rev_alpesh"]

    # -----------------------------------------------------------------------
    # Accuracy summary
    # -----------------------------------------------------------------------
    def accuracy_row(errors: pd.Series, label: str) -> None:
        correct = (errors == 0).sum()
        over    = (errors  > 0).sum()
        under   = (errors  < 0).sum()
        print(f"  {label}")
        print(f"    Correct : {correct}/{n}  ({correct/n:.1%})")
        print(f"    Over    : {over}/{n}     ({over/n:.1%})")
        print(f"    Under   : {under}/{n}    ({under/n:.1%})")
        print()

    print("=" * 60)
    print("ACCURACY vs Reviewer (Alpesh)")
    print("=" * 60)
    accuracy_row(df["alg_err"],  "LOT algorithm          ")
    accuracy_row(df["cota_err"], "COTA dataset (baseline)")

    # -----------------------------------------------------------------------
    # Flag summary
    # -----------------------------------------------------------------------
    from collections import Counter
    all_flags: list[str] = [f for flags in df["flags"] for f in flags]
    flag_counts = Counter(all_flags)

    print("=" * 60)
    print("FLAG SUMMARY (patients)")
    print("=" * 60)
    flag_labels = {
        "duplicate_line":  "Duplicate consecutive lines (data artefact - deduplicated before algo)",
        "partial_removal": "Partial drug removal merged by algo (COTA/reviewer may split as new line)",
        "non_parseable":   "Contains investigational/unnamed regimen -> human review recommended",
        "cart_therapy":    "Contains CAR-T therapy -> human review recommended",
    }
    for flag, count in flag_counts.most_common():
        print(f"  {flag_labels.get(flag, flag)}: {count} patients")
    print()

    # -----------------------------------------------------------------------
    # Misclassified cases
    # -----------------------------------------------------------------------
    wrong = df[df["alg_err"] != 0].copy()
    wrong["direction"] = wrong["alg_err"].apply(lambda x: "over" if x > 0 else "under")

    print("=" * 60)
    print(f"MISCLASSIFIED  ({len(wrong)} patients)")
    print("=" * 60)

    for direction in ("over", "under"):
        subset = wrong[wrong["direction"] == direction]
        print(f"\n--- {direction.upper()}-COUNTED ({len(subset)}) ---\n")
        for _, row in subset.iterrows():
            print(
                f"  rev={int(row['rev_alpesh'])}  "
                f"alg={int(row['lot_alg'])}  "
                f"cota={int(row['cota_lot'])}  "
                f"({direction} by {abs(int(row['alg_err']))})"
                f"  flags={row['flags']}"
            )
            for line_str in row["trt_strings"]:
                print(f"    {line_str.strip()}")
            print()

    # -----------------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------------
    output_cols = [
        "cpid", "rev_alpesh", "rev_alberto", "cota_lot",
        "lot_alg", "alg_err", "cota_err", "flags",
    ]
    out = df[output_cols].copy()
    out["trt_sequence"] = df["trt_strings"].apply(lambda lines: " | ".join(str(l).strip() for l in lines))
    out = out[output_cols + ["trt_sequence"]]
    out["flags"] = out["flags"].apply(lambda f: ", ".join(f) if f else "")
    out.to_excel("lot_results_cota.xlsx", index=False)
    print("Results saved to lot_results_cota.xlsx")


if __name__ == "__main__":
    main()