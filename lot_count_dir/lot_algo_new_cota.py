"""
Deterministic Line-of-Therapy (LoT) counter for the new COTA export.

Input : Output/COTA_cleaned.xlsx (written by preprocessing_new_cota.py)
Output: Output/COTA_cleaned_with_LOT.xlsx

How it works
------------
Every preprocessed row is one vendor line of therapy. The engine walks each
patient's rows in order and, for every consecutive pair, decides whether the
later row continues the current LoT (MERGE) or starts a new one (NEW). Rules
are applied first-match-wins, in this order:

    P1 (confirmed progression)      -> NEW
    P1b (planned sequential)        -> MERGE
    Rule 1 (steroid-only segment)   -> MERGE
    Rule 1 (steroid-only 1st seg.)  -> MERGE
    Mandatory-drug planned triplet  -> MERGE
    Rules 2+3 (identical drugs)     -> MERGE
    Pre-ASCT re-induction           -> MERGE
    ASCT rule                       -> MERGE
    CAR-T rule                      -> MERGE
    P3 (maintenance after combo)    -> MERGE
    Default                         -> NEW  (flagged for review)

Because the engine can only ever merge vendor lines, "agreement with the
vendor" means keeping the vendor split. The rules that merge are calibrated
against the reviewer-adjudicated subset (see score_vs_adjudication.py).

Change log (Sept 2026, after review of the adjudicated data)
------------------------------------------------------------
* P1 no longer defers to the ASCT rule when the prior line contains a
  transplant, and no longer defers identical-drug retreatment to Rules 2+3.
  A documented progression is a new LoT unless it is a <=7-day renewal.
* ASCT rule rewritten. It used to merge any transition adjacent to a
  transplant line within 180 days regardless of drugs or progression; in this
  export the transplant line usually already contains the maintenance, so the
  following vendor line is a genuinely new LoT. It now merges only
  (a) induction -> transplant line whose pre-transplant drugs come from the
      induction, or
  (b) transplant line -> maintenance whose drugs are a subset of that line,
  and never across a documented progression or a gap > 180 days.
* CAR-T rule direction fixed. The "post-CAR-T exclusion" was being applied to
  the segment BEFORE CAR-T (bridging therapy), which reviewers count as its
  own line. Lymphodepletion right before infusion still merges.
* Rules 2+3 compare the full active drug set as well as the tail chunk, so
  [DRd],[Dd] -> [DRd] (a drug held for toxicity and resumed) merges.
* Counting pass now tracks the current LoT start date (the mandatory-triplet
  rule measured from the patient's very first line, so it never fired). The
  pre-ASCT re-induction rule used to compare the prior row to itself (so it
  never fired); it now fires at the induction -> re-induction transition by
  looking one row ahead for the transplant.
* Steroid-only FIRST segment absorption (spec rule 1) is handled here; the
  preprocessing step only aligns the start date.
* Output workbook is written once (it was written twice; the second write
  dropped the Patient_Summary sheet and the pivot's regimen index).
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "Output" / "COTA_cleaned.xlsx"
LEGACY_INPUT_PATH = BASE_DIR / "Output" / "COTA_misclassified_rows_UPD_newrules.xlsx"
OUTPUT_PATH = INPUT_PATH.with_name("COTA_cleaned_with_LOT.xlsx")

# -----------------------------------------------------------------------------
# Rule parameters (days unless stated)
# -----------------------------------------------------------------------------
P1_PD_WINDOW_DAYS = 60             # PD must be recorded within this window before the next line
PRESCRIPTION_RENEWAL_GAP_DAYS = 7  # same drugs restarted within this gap = renewal, not a new line
IDENTICAL_DRUGS_MAX_GAP_DAYS = 180
IDENTICAL_DRUGS_USE_FULL_SET = True  # False = tail chunk only (original spec wording)
ASCT_MAX_GAP_DAYS = 180
PRE_ASCT_REINDUCTION_MAX_DAYS = 90
CAR_T_CONDITIONING_MAX_DAYS = 14   # lymphodepletion start -> CAR-T start
POST_CAR_T_ABSORB_DAYS = 30        # segment starting within this many days after CAR-T is absorbed
MANDATORY_TRIPLET_WINDOW_DAYS = 45 # measured from the CURRENT LoT start
STEROID_FIRST_SEGMENT_GAP_DAYS = 7
P1B_MAX_GAP_DAYS = 3
P3_MAX_GAP_DAYS = 30

STEROIDS = {"dexamethasone", "prednisone", "prednisolone", "methylprednisolone"}
CONDITIONING = {"melphalan", "busulfan", "carmustine"}
LYMPHODEPLETION = {"cyclophosphamide", "fludarabine", "bendamustine"}
MANDATORY_CLASS_DRUGS = {"daratumumab", "isatuximab", "carfilzomib", "ixazomib"}
MAINTENANCE_DRUGS = {"lenalidomide", "thalidomide", "daratumumab", "ixazomib"}
PD_REASON_TOKENS = {"progression", "progressed", "refractory", "pd"}

# -----------------------------------------------------------------------------
# Colors / flags (kept for the colour-coded review workbook)
# -----------------------------------------------------------------------------
PATTERN_COLORS = {
    "p1_confirmed_progression": "FF9999",
    "p1b_planned_sequential": "99CCFF",
    "rule1_steroid_only_absorbed": "CCFFCC",
    "rule1_first_segment_steroid_absorbed": "CCFFCC",
    "mandatory_drug_planned_triplet": "99FF99",
    "rules_2_3_identical_active_drugs": "CCCCFF",
    "pre_asct_reinduction": "FFCC99",
    "asct_rule_post_transplant": "FFCCCC",
    "car_t_rule": "FF99CC",
    "p3_maintenance_after_combination": "CCFF99",
    "default_new_lot": "FFE6CC",
    "START": "FFFFFF",
}

# Patterns whose decision is considered safe (review_flag = 1). Everything
# else, including the default NEW-LoT fallback, is flagged for review (2).
FIXABLE_PATTERNS = {
    "p1_confirmed_progression",
    "p1b_planned_sequential",
    "rule1_steroid_only_absorbed",
    "rule1_first_segment_steroid_absorbed",
    "mandatory_drug_planned_triplet",
    "rules_2_3_identical_active_drugs",
    "pre_asct_reinduction",
    "asct_rule_post_transplant",
    "car_t_rule",
    "p3_maintenance_after_combination",
}

# -----------------------------------------------------------------------------
# Text / regimen helpers
# -----------------------------------------------------------------------------

def normalize_text(value) -> str:
    """Removes extra whitespace and normalizes text strings."""
    if value is None or (isinstance(value, float) and np.isnan(value)) or (not isinstance(value, str) and pd.isna(value)):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def clean_category(value) -> str:
    text = normalize_text(value)
    return text.replace(" (Not in provided Fiona's category list)", "").strip()


def clean_token(value: str) -> str:
    token = normalize_text(value).lower()
    token = token.strip().strip("[]'").strip('"')
    token = token.replace("[", "").replace("]", "")
    return re.sub(r"\s+", " ", token).strip()


def split_list_string(value) -> list[str]:
    """Split values like 'a, b, c' or 'A + B + C' into clean lowercase tokens."""
    text = normalize_text(value)
    if not text:
        return []
    text = text.replace(" (Not in provided Fiona's category list)", "")
    raw_parts = text.split(" + ") if " + " in text else text.split(",")
    parts: list[str] = []
    for part in raw_parts:
        token = clean_token(part)
        if token and token not in parts:
            parts.append(token)
    return parts


def normalize_regimen_string(value) -> str:
    text = normalize_text(value).lower()
    if not text:
        return ""
    text = text.replace("[ ", "[").replace(" ]", "]")
    text = re.sub(r"\s*,\s*", ", ", text)
    text = re.sub(r"\]\s*,\s*\[", "], [", text)
    return re.sub(r"\s+", " ", text).strip()


def parse_regimen_chunks(value) -> list[tuple[str, ...]]:
    """Parse bracket-level regimen chunks: '[a, b], [c]' -> [('a','b'), ('c',)]."""
    text = normalize_regimen_string(value)
    if not text:
        return []
    bracket_chunks = re.findall(r"\[([^\]]+)\]", text)
    if not bracket_chunks:
        bracket_chunks = [text.replace("[", "").replace("]", "")]
    chunks: list[tuple[str, ...]] = []
    for chunk in bracket_chunks:
        drugs: list[str] = []
        for raw in chunk.split(","):
            drug = clean_token(raw)
            if drug and drug not in drugs:
                drugs.append(drug)
        if drugs:
            chunks.append(tuple(drugs))
    return chunks


def chunk_list_to_str(chunks: list[tuple[str, ...]]) -> str:
    return " | ".join("[" + ", ".join(chunk) + "]" for chunk in chunks)


def get_drug_set_from_chunks(chunks: list[tuple[str, ...]]) -> set[str]:
    return {drug for chunk in chunks for drug in chunk}


def is_transplant_token(token: str) -> bool:
    return "sct" in token or "transplant" in token


def get_drug_set(row: pd.Series) -> set[str]:
    chunks = parse_regimen_chunks(row.get("line_of_therapy_name"))
    if chunks:
        return get_drug_set_from_chunks(chunks)
    if "parsed_drugs" in row and normalize_text(row.get("parsed_drugs")):
        return set(split_list_string(row.get("parsed_drugs")))
    return set()


def _active_filter(drugs: set[str]) -> set[str]:
    return {d for d in drugs if d not in STEROIDS and d not in CONDITIONING and not is_transplant_token(d)}


def _list_column(row: pd.Series, column: str) -> set[str] | None:
    """Read a pre-processing list column. Returns None only when the column is
    absent. An empty cell (''/NaN after an Excel round trip) is an empty set,
    NOT a signal to fall back to re-parsing the regimen."""
    if column not in row:
        return None
    value = row.get(column)
    if isinstance(value, (list, set, tuple)):
        return set(value)
    if isinstance(value, str):
        return set(split_list_string(value)) if value.strip() else set()
    return set()  # NaN / None


def get_active_drugs(row: pd.Series) -> set[str]:
    """drugs_active from pre-processing (steroids, conditioning and procedures removed)."""
    active = _list_column(row, "drugs_active")
    if active is not None:
        return active
    return _active_filter(get_drug_set(row))


def get_tail_active_drugs(row: pd.Series) -> set[str]:
    """tail_drugs_active from pre-processing (active drugs of the last bracket chunk)."""
    tail = _list_column(row, "tail_drugs_active")
    if tail is not None:
        return tail
    chunks = parse_regimen_chunks(row.get("line_of_therapy_name"))
    return _active_filter(set(chunks[-1])) if chunks else set()


def get_pre_asct_active_drugs(row: pd.Series) -> set[str]:
    """Active drugs in the bracket chunks that come BEFORE the transplant chunk
    (mobilisation / continued induction). Conditioning agents are excluded."""
    drugs: set[str] = set()
    for chunk in parse_regimen_chunks(row.get("line_of_therapy_name")):
        if any(is_transplant_token(t) for t in chunk):
            break
        drugs |= _active_filter(set(chunk))
    return drugs


def get_family_set(row: pd.Series) -> set[str]:
    if "mapped_family_classes" in row and normalize_text(row.get("mapped_family_classes")):
        return set(split_list_string(row.get("mapped_family_classes")))
    return set(split_list_string(clean_category(row.get("reconstructed_family_combination"))))


def set_to_str(values: set[str]) -> str:
    return ", ".join(sorted(values))


def get_regimen_signature(row: pd.Series) -> dict:
    chunks = parse_regimen_chunks(row.get("line_of_therapy_name"))
    return {
        "normalized_regimen": normalize_regimen_string(row.get("line_of_therapy_name")),
        "chunks": chunks,
        "drug_set": get_drug_set(row),
        "active_drugs": get_active_drugs(row),
        "tail_active_drugs": get_tail_active_drugs(row),
        "family_set": get_family_set(row),
    }


def calculate_gap_days(date1, date2) -> float:
    """Days from date1 to date2 (negative if date2 is earlier). NaN if either is missing."""
    if date1 is None or date2 is None:
        return np.nan
    try:
        d1 = pd.to_datetime(date1, errors="coerce")
        d2 = pd.to_datetime(date2, errors="coerce")
    except Exception:
        return np.nan
    if pd.isna(d1) or pd.isna(d2):
        return np.nan
    return float((d2 - d1).days)


def split_discontinue_reasons(value) -> list[str]:
    """The vendor stores one reason per regimen chunk, separated by backslashes,
    e.g. 'toxicity\\toxicity\\progression'. Whitespace inside a token
    ('doctor _preference') is export noise and is removed."""
    text = normalize_text(value).lower()
    if not text:
        return []
    return [re.sub(r"\s+", "", part) for part in text.split("\\") if part.strip()]


def last_discontinue_reason(row: pd.Series) -> str:
    reasons = split_discontinue_reasons(row.get("discontinue_reason") if row is not None else None)
    return reasons[-1] if reasons else ""


def has_pd_signal(row: pd.Series | None) -> bool:
    """True if any recorded discontinuation reason on the row is a progression /
    refractory signal. Token-based, so 'pandemic_epidemic' cannot match 'pd'."""
    if row is None:
        return False
    return bool(set(split_discontinue_reasons(row.get("discontinue_reason"))) & PD_REASON_TOKENS)


def segment_gap_days(row: pd.Series, prev_row: pd.Series) -> float:
    return calculate_gap_days(prev_row.get("date_end_line_of_therapy"), row.get("date_start_line_of_therapy"))


# -----------------------------------------------------------------------------
# Rules. Each returns (fired, explanation).
# -----------------------------------------------------------------------------

def rule_p1_confirmed_progression(row: pd.Series, prev_row: pd.Series, gap_days: float) -> tuple[bool, str]:
    """
    P1: Confirmed Progression -> New LoT.

    Trigger: the prior line records progression / refractory disease and the
    next line starts within P1_PD_WINDOW_DAYS of it ending.

    Exception: identical active drugs restarted within
    PRESCRIPTION_RENEWAL_GAP_DAYS (a renewal, not a new line).

    Removed exceptions (they hid real progressions in the adjudicated data):
    - "prior line contains ASCT -> leave to the ASCT rule". The transplant
      line here usually already contains the maintenance, so the next vendor
      line after a documented progression is a new LoT.
    - "identical drugs with gap <= 180 -> leave to Rules 2+3". Retreatment
      with the same regimen after documented progression is a new LoT.
    """
    if prev_row is None:
        return False, "no previous row to check for PD"
    if not has_pd_signal(prev_row):
        return False, "no PD/refractory signal in previous segment's discontinuation reason"
    if pd.isna(gap_days):
        return False, "gap days unknown - cannot verify PD window"
    if gap_days > P1_PD_WINDOW_DAYS:
        return False, f"PD signal detected but gap of {gap_days:.0f} days exceeds {P1_PD_WINDOW_DAYS}-day window"

    prev_active = get_active_drugs(prev_row)
    curr_active = get_active_drugs(row)
    if gap_days <= PRESCRIPTION_RENEWAL_GAP_DAYS and prev_active == curr_active:
        return False, f"prescription renewal: same active drugs, gap {gap_days:.0f} days <= {PRESCRIPTION_RENEWAL_GAP_DAYS}"

    return True, f"P1 fired: PD/refractory in previous segment within {gap_days:.0f} days (<= {P1_PD_WINDOW_DAYS}-day window)"


def rule_p1b_planned_sequential(row: pd.Series, prev_row: pd.Series, gap_days: float) -> tuple[bool, str]:
    """
    P1b: Planned Sequential Therapy -> Same LoT.
    Gap <= 3 days, next segment is single-agent LEN/THAL, prior segment had
    >= 2 active drugs and lasted 60-365 days, and no PD signal.
    """
    if pd.isna(gap_days) or gap_days > P1B_MAX_GAP_DAYS:
        return False, f"gap > {P1B_MAX_GAP_DAYS} days"
    curr_active = get_active_drugs(row)
    if len(curr_active) != 1:
        return False, f"not single-agent: {len(curr_active)} active drugs"
    curr_drug = next(iter(curr_active))
    if curr_drug not in {"lenalidomide", "thalidomide"}:
        return False, f"not LEN or THAL: {curr_drug}"
    prev_active = get_active_drugs(prev_row)
    if len(prev_active) < 2:
        return False, f"prior segment has {len(prev_active)} active drugs (< 2)"
    duration = prev_row.get("segment_duration_days")
    if not pd.isna(duration) and (duration < 60 or duration > 365):
        return False, f"prior segment duration {duration} days (not 60-365)"
    if has_pd_signal(row) or has_pd_signal(prev_row):
        return False, "PD signal present - cannot treat as planned sequential"
    return True, "planned induction -> maintenance sequence (e.g., VRd -> LEN)"


def rule_steroid_only_absorbed(row: pd.Series, prev_row: pd.Series) -> tuple[bool, str]:
    """Rule 1: a steroid-only segment is absorbed into the prior LoT (unless the prior segment is CAR-T)."""
    if not bool(row.get("is_steroid_only", False)):
        return False, "not steroid-only segment"
    if bool(prev_row.get("is_car_t", False)):
        return False, "prior segment is CAR-T"
    return True, "steroid-only segment absorbed into prior LOT"


def rule_first_segment_steroid_absorbed(row: pd.Series, prev_row: pd.Series, gap_days: float,
                                        prev_is_first_segment: bool) -> tuple[bool, str]:
    """
    Rule 1 (first segment): a steroid-only FIRST segment followed within
    STEROID_FIRST_SEGMENT_GAP_DAYS by an active segment is part of that
    segment's LoT (same LoT start). preprocessing_new_cota.py aligns the
    start date; the merge itself happens here.
    """
    if not prev_is_first_segment:
        return False, "prior segment is not the patient's first segment"
    if not bool(prev_row.get("is_steroid_only", False)):
        return False, "prior segment is not steroid-only"
    if bool(row.get("is_steroid_only", False)):
        return False, "current segment is also steroid-only"
    if pd.isna(gap_days) or gap_days > STEROID_FIRST_SEGMENT_GAP_DAYS:
        return False, f"gap {gap_days} days > {STEROID_FIRST_SEGMENT_GAP_DAYS}"
    return True, "steroid-only first segment absorbed into first active segment"


def rule_mandatory_drug_planned_triplet(row: pd.Series, prev_row: pd.Series,
                                        lot_start_date, gap_from_lot_start: float) -> tuple[bool, str]:
    """
    Mandatory-Drug Planned Triplet -> Same LoT.
    The next segment adds exactly one drug from {dara, isa, carfilzomib,
    ixazomib}, removes nothing, and does so within MANDATORY_TRIPLET_WINDOW_DAYS
    of the CURRENT LoT's start.
    """
    prev_active = get_active_drugs(prev_row)
    curr_active = get_active_drugs(row)
    if not prev_active.issubset(curr_active):
        return False, "drugs removed from prior segment"
    added = curr_active - prev_active
    if len(added) != 1:
        return False, f"{len(added)} drugs added"
    added_drug = next(iter(added))
    if added_drug not in MANDATORY_CLASS_DRUGS:
        return False, f"added drug {added_drug} not in mandatory class"
    if pd.isna(gap_from_lot_start) or gap_from_lot_start > MANDATORY_TRIPLET_WINDOW_DAYS:
        return False, f"addition at day {gap_from_lot_start} > {MANDATORY_TRIPLET_WINDOW_DAYS} days from LoT start"
    return True, f"planned triplet: added {added_drug} within {MANDATORY_TRIPLET_WINDOW_DAYS} days of LoT start"


def rule_rules_2_3_identical_active_drugs(row: pd.Series, prev_row: pd.Series,
                                          gap_days: float) -> tuple[bool, str]:
    """
    Rules 2+3: Identical Active Drugs -> Same LoT.

    Trigger: the active drugs are identical and the gap is
    <= IDENTICAL_DRUGS_MAX_GAP_DAYS. "Identical" means the tail chunks match,
    or (IDENTICAL_DRUGS_USE_FULL_SET) the full active drug sets of the two
    lines match. The second form covers the common vendor shape
    [DRd],[Dd] -> [DRd], where one drug was held for toxicity and resumed;
    reviewers merged that shape ~80% of the time when no progression was
    recorded. (A documented progression is caught by P1 first.)
    """
    prev_tail = get_tail_active_drugs(prev_row)
    curr_tail = get_tail_active_drugs(row)
    prev_active = get_active_drugs(prev_row)
    curr_active = get_active_drugs(row)

    tail_match = bool(prev_tail) and prev_tail == curr_tail
    full_match = IDENTICAL_DRUGS_USE_FULL_SET and bool(curr_active) and prev_active == curr_active
    if not (tail_match or full_match):
        return False, "active drugs not identical"
    if pd.isna(gap_days) or gap_days > IDENTICAL_DRUGS_MAX_GAP_DAYS:
        return False, f"gap {gap_days} days > {IDENTICAL_DRUGS_MAX_GAP_DAYS}"
    if bool(prev_row.get("is_asct", False)):
        return False, "prior segment contains ASCT (handled by ASCT rule)"
    how = "tail chunk" if tail_match else "full active drug set"
    return True, f"identical active drugs ({how}) with gap <= {IDENTICAL_DRUGS_MAX_GAP_DAYS} days"


def rule_pre_asct_reinduction(row: pd.Series, prev_row: pd.Series,
                              next_row: pd.Series | None) -> tuple[bool, str]:
    """
    Pre-ASCT Re-induction -> Same LoT.

    Spec: when the segment immediately before a transplant (re-induction) uses
    a proper subset of the induction drugs and started within
    PRE_ASCT_REINDUCTION_MAX_DAYS of the transplant, induction, re-induction
    and transplant are one LoT.

    Because the engine decides transitions in order, this rule is evaluated at
    the induction -> re-induction transition and looks one row AHEAD for the
    transplant. The re-induction -> transplant transition is then merged by
    the ASCT rule (shape a).
    """
    if next_row is None or not bool(next_row.get("is_asct", False)):
        return False, "next segment is not ASCT"
    if bool(row.get("is_asct", False)) or bool(prev_row.get("is_asct", False)):
        return False, "current or prior segment already contains ASCT"
    reinduction_active = get_active_drugs(row)
    induction_active = get_active_drugs(prev_row)
    if not reinduction_active or not reinduction_active.issubset(induction_active) or reinduction_active == induction_active:
        return False, "re-induction drugs not a proper subset of induction drugs"
    if has_pd_signal(prev_row) or has_pd_signal(row):
        return False, "PD recorded before transplant"
    gap = calculate_gap_days(row.get("date_start_line_of_therapy"), next_row.get("date_start_line_of_therapy"))
    if pd.isna(gap) or gap > PRE_ASCT_REINDUCTION_MAX_DAYS:
        return False, f"re-induction start -> ASCT gap {gap} days > {PRE_ASCT_REINDUCTION_MAX_DAYS}"
    return True, "pre-ASCT re-induction (drugs dropped before transplant, transplant follows within 90 days)"


def rule_asct_rule(row: pd.Series, prev_row: pd.Series, asct_index=np.nan) -> tuple[bool, str]:
    """
    ASCT rule -> Same LoT, for two shapes only:

    (a) induction line -> transplant line. Merge when every active drug that
        appears BEFORE the transplant chunk (mobilisation, continued induction)
        was already part of the induction line. A transplant line that starts a
        new regimen before conditioning is salvage, i.e. a new LoT.
    (b) transplant line -> maintenance line. Merge when the maintenance line's
        active drugs are a (non-empty) subset of the transplant line's drugs.
        Anything new after transplant is a new LoT.

    Both require: no documented progression on the prior line, and a gap
    <= ASCT_MAX_GAP_DAYS. Transplant -> transplant never merges here.

    Evidence (adjudicated COTA): the previous version merged every
    ASCT-adjacent transition within 180 days and was wrong 17 of 24 times.
    """
    cur_asct = bool(row.get("is_asct", False))
    prev_asct = bool(prev_row.get("is_asct", False))
    if not cur_asct and not prev_asct:
        return False, "not ASCT-related transition"
    gap = segment_gap_days(row, prev_row)
    if pd.isna(gap) or gap > ASCT_MAX_GAP_DAYS:
        return False, f"gap {gap} days > {ASCT_MAX_GAP_DAYS} or unknown"
    if has_pd_signal(prev_row):
        return False, "progression recorded on the prior line"

    prev_active = get_active_drugs(prev_row)
    curr_active = get_active_drugs(row)
    idx = "" if pd.isna(asct_index) else f" (asct_index={int(asct_index)})"

    if cur_asct and not prev_asct:
        pre_asct = get_pre_asct_active_drugs(row)
        if pre_asct.issubset(prev_active):
            return True, f"induction -> transplant line{idx}: pre-ASCT drugs {set_to_str(pre_asct) or 'conditioning only'} come from induction"
        return False, f"transplant line starts new drugs before ASCT: {set_to_str(pre_asct - prev_active)}"

    if prev_asct and not cur_asct:
        if curr_active and curr_active.issubset(prev_active):
            return True, f"post-transplant maintenance{idx}: {set_to_str(curr_active)} subset of transplant line"
        return False, "post-ASCT line adds drugs not in the transplant line"

    return False, "transplant -> transplant"


def rule_car_t_rule(row: pd.Series, prev_row: pd.Series) -> tuple[bool, str]:
    """
    CAR-T rule -> Same LoT, for two shapes only:

    (a) lymphodepletion -> CAR-T: the prior line is only lymphodepleting
        chemotherapy (cyclophosphamide / fludarabine / bendamustine) that
        started within CAR_T_CONDITIONING_MAX_DAYS of the CAR-T line.
        Bridging therapy before CAR-T is its own line (reviewers split 4/5).
    (b) CAR-T -> next segment starting within POST_CAR_T_ABSORB_DAYS of the
        CAR-T line ending, with no progression recorded.

    The previous version applied (b) to the segment BEFORE CAR-T, merging the
    bridging regimen into the CAR-T line whenever the gap was <= 30 days.
    """
    cur_cart = bool(row.get("is_car_t", False))
    prev_cart = bool(prev_row.get("is_car_t", False))
    if not cur_cart and not prev_cart:
        return False, "not CAR-T transition"

    if cur_cart and not prev_cart:
        prev_active = get_active_drugs(prev_row)
        start_gap = calculate_gap_days(prev_row.get("date_start_line_of_therapy"), row.get("date_start_line_of_therapy"))
        if prev_active and prev_active.issubset(LYMPHODEPLETION) and not pd.isna(start_gap) and start_gap <= CAR_T_CONDITIONING_MAX_DAYS:
            return True, f"lymphodepletion merged into CAR-T (started {start_gap:.0f} days before infusion)"
        return False, "bridging therapy before CAR-T is its own line"

    if prev_cart and not cur_cart:
        gap = segment_gap_days(row, prev_row)
        if has_pd_signal(prev_row):
            return False, "progression recorded on the CAR-T line"
        if not pd.isna(gap) and gap <= POST_CAR_T_ABSORB_DAYS:
            return True, f"post-CAR-T exclusion (starts {gap:.0f} days after CAR-T, <= {POST_CAR_T_ABSORB_DAYS})"
        return False, f"segment starts {gap} days after CAR-T"

    return False, "CAR-T -> CAR-T"


def rule_p3_maintenance_after_combination(row: pd.Series, prev_row: pd.Series, gap_days: float) -> tuple[bool, str]:
    """
    P3: Maintenance After Combination -> Same LoT.
    Next segment is a single maintenance drug (LEN/THAL/DARA/IXA, or bortezomib
    with its own criteria) that was in the prior combination's tail, the prior
    combination lasted >= 30 days, no PD, and the gap is within limits.
    """
    prev_active = get_active_drugs(prev_row)
    curr_active = get_active_drugs(row)
    if not curr_active.issubset(prev_active) or curr_active == prev_active:
        return False, "current active drugs not proper subset of prior"
    if has_pd_signal(prev_row):
        return False, "PD in prior segment - this is not maintenance"
    if len(curr_active) != 1:
        return False, f"not single-agent: {len(curr_active)} active drugs"
    curr_drug = next(iter(curr_active))
    is_bortezomib = curr_drug == "bortezomib"
    if not is_bortezomib and curr_drug not in MAINTENANCE_DRUGS:
        return False, f"not a maintenance drug: {curr_drug}"
    if curr_drug not in get_tail_active_drugs(prev_row):
        return False, f"{curr_drug} not in prior tail_active_drugs"
    if len(prev_active) < 2:
        return False, "prior segment not a combination"
    duration = prev_row.get("segment_duration_days")
    if not pd.isna(duration) and duration < 30:
        return False, f"prior segment duration {duration} days (< 30)"

    if is_bortezomib:
        if bool(prev_row.get("is_asct", False)):
            if pd.isna(gap_days) or gap_days <= 30:
                return True, "post-ASCT bortezomib maintenance (<= 30 days, DEX allowed)"
        elif not pd.isna(gap_days) and gap_days <= 14:
            if "dexamethasone" not in get_drug_set(row):
                return True, "non-ASCT bortezomib maintenance (<= 14 days, no DEX)"
            return False, "non-ASCT bortezomib with DEX (not allowed)"
    elif not pd.isna(gap_days) and gap_days <= P3_MAX_GAP_DAYS:
        return True, f"maintenance {curr_drug} (gap <= {P3_MAX_GAP_DAYS} days)"
    return False, "gap criteria not met"


# -----------------------------------------------------------------------------
# Transition decision (first match wins)
# -----------------------------------------------------------------------------

def _decision(decision: str, rule: str, pattern: str, explanation: str, needs_review: bool = False) -> dict:
    return {"decision": decision, "rule_applied": rule, "pattern": pattern,
            "explanation": explanation, "needs_review": needs_review}


def determine_transition_type(current_row: pd.Series, previous_row: pd.Series,
                              next_row: pd.Series | None, lot_start_date,
                              gap_from_lot_start: float,
                              previous_is_first_segment: bool = False) -> dict:
    """
    Decide whether current_row continues the LoT of previous_row (MERGE) or
    starts a new one (NEW).

    next_row           : the row after current_row (None if there is none);
                         only the pre-ASCT re-induction rule looks ahead
    lot_start_date     : start date of the LoT previous_row belongs to
    gap_from_lot_start : days from lot_start_date to current_row's start
    """
    gap_days = segment_gap_days(current_row, previous_row)
    asct_index = current_row.get("asct_index", np.nan)

    fired, why = rule_p1_confirmed_progression(current_row, previous_row, gap_days)
    if fired:
        return _decision("NEW", "P1", "p1_confirmed_progression", why)

    fired, why = rule_p1b_planned_sequential(current_row, previous_row, gap_days)
    if fired:
        return _decision("MERGE", "P1b", "p1b_planned_sequential", why)

    fired, why = rule_steroid_only_absorbed(current_row, previous_row)
    if fired:
        return _decision("MERGE", "Rule 1", "rule1_steroid_only_absorbed", why)

    fired, why = rule_first_segment_steroid_absorbed(current_row, previous_row, gap_days, previous_is_first_segment)
    if fired:
        return _decision("MERGE", "Rule 1 (first segment)", "rule1_first_segment_steroid_absorbed", why)

    fired, why = rule_mandatory_drug_planned_triplet(current_row, previous_row, lot_start_date, gap_from_lot_start)
    if fired:
        return _decision("MERGE", "Mandatory Triplet", "mandatory_drug_planned_triplet", why)

    fired, why = rule_rules_2_3_identical_active_drugs(current_row, previous_row, gap_days)
    if fired:
        return _decision("MERGE", "Rules 2+3", "rules_2_3_identical_active_drugs", why)

    fired, why = rule_pre_asct_reinduction(current_row, previous_row, next_row)
    if fired:
        return _decision("MERGE", "Pre-ASCT", "pre_asct_reinduction", why)

    fired, why = rule_asct_rule(current_row, previous_row, asct_index)
    if fired:
        return _decision("MERGE", "ASCT Rule", "asct_rule_post_transplant", why)

    fired, why = rule_car_t_rule(current_row, previous_row)
    if fired:
        return _decision("MERGE", "CAR-T Rule", "car_t_rule", why)

    fired, why = rule_p3_maintenance_after_combination(current_row, previous_row, gap_days)
    if fired:
        return _decision("MERGE", "P3", "p3_maintenance_after_combination", why)

    # A regimen change with no documented progression and no matching merge
    # rule defaults to a NEW LoT (standard convention) and is flagged for
    # review because it is the least-specified case.
    return _decision("NEW", "Default", "default_new_lot",
                     "Default: New LOT - regimen changed, no specific merge rule matched", needs_review=True)


# -----------------------------------------------------------------------------
# Counting pass
# -----------------------------------------------------------------------------

def count_lots_with_rule_engine(df_patient_rows: pd.DataFrame) -> dict:
    """
    Walk one patient's rows in order, deciding MERGE/NEW for each transition
    and numbering LoTs as we go. Because decisions are sequential, the current
    LoT start date is always known when a transition is evaluated (the
    mandatory-triplet rule needs it); the pre-ASCT rule additionally looks
    one row ahead for the transplant.
    """
    sorted_rows = df_patient_rows.sort_values("_original_row_order")
    if len(sorted_rows) == 0:
        return {"total_lots": 0, "lot_assignments": [], "transitions": [], "flags": [], "new_lot_triggers": []}

    lot_assignments: list[dict] = []
    transitions: list[dict] = []
    flags: list[dict] = []
    new_lot_triggers: list[dict] = []

    current_lot = 1
    lot_start_date = None
    prev_row = None

    for idx in range(len(sorted_rows)):
        row = sorted_rows.iloc[idx]
        row_id = row["_original_row_order"]
        next_row = sorted_rows.iloc[idx + 1] if idx + 1 < len(sorted_rows) else None

        if idx == 0:
            lot_start_date = row.get("date_start_line_of_therapy")
            result = _decision("START", "START", "START", "First segment for patient")
            is_new_lot = False
        else:
            gap_from_lot_start = calculate_gap_days(lot_start_date, row.get("date_start_line_of_therapy"))
            result = determine_transition_type(row, prev_row, next_row, lot_start_date,
                                               gap_from_lot_start, previous_is_first_segment=(idx == 1))
            is_new_lot = result["decision"] == "NEW"
            if is_new_lot:
                current_lot += 1
                lot_start_date = row.get("date_start_line_of_therapy")
                new_lot_triggers.append({"row_id": row_id, "lot_number": current_lot, "reason": result["pattern"],
                                         "rule": result["rule_applied"], "explanation": result["explanation"]})
            if result.get("needs_review"):
                flags.append({"row_id": row_id, "reason": result["pattern"], "explanation": result["explanation"]})

        transitions.append({"row_id": row_id, "decision": result["decision"], "rule_applied": result["rule_applied"],
                            "pattern": result["pattern"], "explanation": result["explanation"],
                            "needs_review": result.get("needs_review", False)})
        lot_assignments.append({"row_id": row_id, "lot_number": current_lot,
                                "regimen": row.get("line_of_therapy_name", ""),
                                "is_new_lot": is_new_lot, "lot_start_row": idx == 0 or is_new_lot,
                                "rule_applied": result["rule_applied"], "pattern": result["pattern"],
                                "explanation": result["explanation"]})
        prev_row = row

    return {"total_lots": current_lot, "lot_assignments": lot_assignments, "transitions": transitions,
            "flags": flags, "new_lot_triggers": new_lot_triggers}


def get_lot_count_for_regimen(df_all_patients: pd.DataFrame):
    """Run the counter for every patient. Returns (per-row lot_df, regimen_summary, patient_df)."""
    all_lot_data = []
    all_patient_data = []

    for patient in df_all_patients["cpid"].unique():
        patient_rows = df_all_patients[df_all_patients["cpid"] == patient]
        result = count_lots_with_rule_engine(patient_rows)
        for assignment in result["lot_assignments"]:
            all_lot_data.append({
                "cpid": patient,
                "row_id": assignment["row_id"],
                "lot_number": assignment["lot_number"],
                "total_lots_for_patient": result["total_lots"],
                "regimen": assignment["regimen"],
                "is_new_lot": assignment["is_new_lot"],
                "lot_start_row": assignment["lot_start_row"],
                "rule_applied": assignment["rule_applied"],
                "pattern": assignment["pattern"],
                "explanation": assignment["explanation"],
            })
        all_patient_data.append({
            "cpid": patient,
            "total_lots": result["total_lots"],
            "num_segments": len(patient_rows),
            "flags": len(result["flags"]),
            "new_lot_triggers": len(result["new_lot_triggers"]),
        })

    lot_df = pd.DataFrame(all_lot_data)
    patient_df = pd.DataFrame(all_patient_data)
    regimen_summary = (
        lot_df.groupby(["regimen", "lot_number"])
        .agg({"cpid": "nunique", "row_id": "count"})
        .rename(columns={"cpid": "patients_on_regimen", "row_id": "total_occurrences"})
        .reset_index()
        .sort_values(["regimen", "lot_number"])
    )
    return lot_df, regimen_summary, patient_df


# -----------------------------------------------------------------------------
# Per-row diagnostics (Misclassified_Analysis sheet)
# -----------------------------------------------------------------------------

def build_transition_diagnostics(row: pd.Series, prev_row: pd.Series | None, decision: dict) -> dict:
    """Columns describing the transition into `row`, using the decision the
    counting pass already made, so the two sheets can never disagree."""
    cur_sig = get_regimen_signature(row)
    if prev_row is not None:
        prev_sig = get_regimen_signature(prev_row)
    else:
        prev_sig = {"normalized_regimen": "", "chunks": [], "drug_set": set(), "active_drugs": set(),
                    "tail_active_drugs": set(), "family_set": set()}

    pattern = decision.get("pattern", "START")
    review_flag = 1 if pattern in FIXABLE_PATTERNS else 2
    return {
        "misclassification_pattern": pattern,
        "pattern_explanation": decision.get("explanation", ""),
        "rule_applied": decision.get("rule_applied", "START"),
        "review_flag": review_flag,
        "suggested_corrected_lot": row.get("doctor_lot_numeric_for_transition") if review_flag == 1 else "",
        "previous_line_of_therapy_name": normalize_text(prev_row.get("line_of_therapy_name")) if prev_row is not None else "",
        "previous_discontinue_reason": normalize_text(prev_row.get("discontinue_reason")) if prev_row is not None else "",
        "current_regimen_normalized": cur_sig["normalized_regimen"],
        "previous_regimen_normalized": prev_sig["normalized_regimen"],
        "current_regimen_chunks": chunk_list_to_str(cur_sig["chunks"]),
        "previous_regimen_chunks": chunk_list_to_str(prev_sig["chunks"]),
        "current_drugs_clean": set_to_str(cur_sig["drug_set"]),
        "previous_drugs_clean": set_to_str(prev_sig["drug_set"]),
        "current_active_drugs": set_to_str(cur_sig["active_drugs"]),
        "previous_active_drugs": set_to_str(prev_sig["active_drugs"]),
        "added_drugs_vs_previous": set_to_str(cur_sig["drug_set"] - prev_sig["drug_set"]),
        "removed_drugs_vs_previous": set_to_str(prev_sig["drug_set"] - cur_sig["drug_set"]),
        "current_families_clean": set_to_str(cur_sig["family_set"]),
        "previous_families_clean": set_to_str(prev_sig["family_set"]),
        "added_families_vs_previous": set_to_str(cur_sig["family_set"] - prev_sig["family_set"]),
        "removed_families_vs_previous": set_to_str(prev_sig["family_set"] - cur_sig["family_set"]),
        "gap_days_to_previous": segment_gap_days(row, prev_row) if prev_row is not None else np.nan,
    }


# -----------------------------------------------------------------------------
# Pipeline
# -----------------------------------------------------------------------------

def load_input(input_path: Path | None = None) -> tuple[pd.DataFrame, Path]:
    path = input_path or (INPUT_PATH if INPUT_PATH.exists() else LEGACY_INPUT_PATH)
    if not path.exists():
        raise FileNotFoundError(
            f"No preprocessing output found. Expected {INPUT_PATH} "
            f"(run preprocessing_new_cota.py first) or {LEGACY_INPUT_PATH}."
        )
    xl = pd.ExcelFile(path)
    sheet = next((n for n in ("All_Data", "All_COTA_With_Transitions") if n in xl.sheet_names), None)
    if sheet is None:
        raise ValueError(f"{path.name} must contain an 'All_Data' or 'All_COTA_With_Transitions' sheet. Found: {xl.sheet_names}")
    all_rows = pd.read_excel(path, sheet_name=sheet)
    all_rows.columns = all_rows.columns.str.strip()
    print(f"Input: {path.name} (all rows from '{sheet}')")
    print(f"Loaded {len(all_rows)} total rows")

    missing = [c for c in ("cpid", "_original_row_order") if c not in all_rows.columns]
    if missing:
        raise ValueError(f"Missing required columns in input workbook: {missing}")
    if "doctor_lot_numeric_for_transition" not in all_rows.columns:
        all_rows["doctor_lot_numeric_for_transition"] = 1  # legacy placeholder, unused by the rule engine
    all_rows = all_rows.sort_values(["cpid", "_original_row_order"]).reset_index(drop=True)
    all_rows["_row_key"] = all_rows["cpid"].astype(str) + "__" + all_rows["_original_row_order"].astype(str)
    return all_rows, path


def run_pipeline(input_path: Path | None = None, output_path: Path = OUTPUT_PATH) -> pd.DataFrame:
    all_rows, _ = load_input(input_path)

    lot_df, regimen_summary, patient_df = get_lot_count_for_regimen(all_rows)

    lot_by_row = lot_df.set_index("row_id")
    all_rows_with_lot = all_rows.copy()
    all_rows_with_lot["lot_number"] = all_rows_with_lot["_original_row_order"].map(lot_by_row["lot_number"])
    all_rows_with_lot["total_lots_patient"] = all_rows_with_lot["cpid"].map(
        lot_df.groupby("cpid")["total_lots_for_patient"].first())
    for col in ("rule_applied", "pattern", "explanation"):
        all_rows_with_lot[col] = all_rows_with_lot["_original_row_order"].map(lot_by_row[col])

    # Per-row diagnostics, reusing the counting pass decisions.
    diagnostics = []
    prev_row = None
    prev_cpid = None
    for _, row in all_rows.iterrows():
        if row["cpid"] != prev_cpid:
            prev_row = None
            prev_cpid = row["cpid"]
        decision = lot_by_row.loc[row["_original_row_order"]]
        diagnostics.append(build_transition_diagnostics(row, prev_row, decision))
        prev_row = row
    misclassified_analysis = pd.concat([all_rows.reset_index(drop=True), pd.DataFrame(diagnostics)], axis=1)

    print("\n=== LOT Assignment Summary ===")
    print(f"Total patients processed: {all_rows_with_lot['cpid'].nunique()}")
    print(f"Total rows processed: {len(all_rows_with_lot)}")
    print(f"Rows with LOT numbers assigned: {all_rows_with_lot['lot_number'].notna().sum()}")
    print(f"Transitions by rule:\n{all_rows_with_lot['rule_applied'].value_counts().to_string()}")

    pivot = (regimen_summary.pivot(index="regimen", columns="lot_number", values="patients_on_regimen")
             .fillna(0).astype(int))

    print(f"\nSaving results to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        all_rows_with_lot.to_excel(writer, sheet_name="All_Data_With_LOT", index=False)
        misclassified_analysis.to_excel(writer, sheet_name="Misclassified_Analysis", index=False)
        lot_df.to_excel(writer, sheet_name="Detailed_LoT_Assignments", index=False)
        regimen_summary.to_excel(writer, sheet_name="Regimen_LoT_Summary", index=False)
        pivot.to_excel(writer, sheet_name="Regimen_LoT_Pivot")
        patient_df.sort_values("cpid").to_excel(writer, sheet_name="Patient_Summary", index=False)
    print(f"Successfully saved to: {output_path}")
    print("Sheets: All_Data_With_LOT, Misclassified_Analysis, Detailed_LoT_Assignments, "
          "Regimen_LoT_Summary, Regimen_LoT_Pivot, Patient_Summary")
    return all_rows_with_lot


if __name__ == "__main__":
    run_pipeline()
