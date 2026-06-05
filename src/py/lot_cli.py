"""
LOT Adjudication Assistant  -  CLI
==================================

Take a Multiple-Myeloma patient's ordered treatment-regimen sequence, count the
number of Lines of Therapy (LOT) with the rule-based algorithm, attach a
*calibrated* confidence, and make a triage decision:

    AUTO-ACCEPT  (high confidence)   vs   NEEDS DOCTOR REVIEW  (low confidence)

CONFIDENCE COMES FROM THE REGIMEN ITSELF (the deployable, works-on-new-data path).
A real new patient has no answer key and no second opinion, so by default the
confidence is read from the patient's own *flag profile* - does the sequence
contain anything the rules struggle with?  This needs nothing but the regimen
you type in:

    no edge-case flags                      -> HIGH   (auto-accept)
    a soft flag (partial drug removal etc.) -> MEDIUM (recommend review)
    a hard flag (CAR-T / investigational)   -> LOW    (needs a doctor)

    python3 src/py/lot_cli.py                      # interactive
    python3 src/py/lot_cli.py --regimens "[bortezomib, dexamethasone, lenalidomide] | [bortezomib, dexamethasone]"
    python3 src/py/lot_cli.py --regimens "[a, b] | [c]" --vendor 2   # optional agreement check
    echo "[a, b] | [c]" | python3 src/py/lot_cli.py

Regimen input format: one "line entry" per patient line, separated by " | ".
Each line entry lists its drugs, optionally in [brackets], e.g.
    "[bortezomib, dexamethasone] | [lenalidomide]"
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CALIB_PATH = os.path.join(HERE, "calibration.json")

# ---------------------------------------------------------------------------
# Algorithm  (faithful copy of the pure functions in textbook_algo_cota.py so
# this CLI stays dependency-free and standard-library only)
# ---------------------------------------------------------------------------

SCT_TOKENS = {
    "autologous sct", "autologous sct1", "autologous sct2",
    "autologoussct", "autologoussct2", "allogeneic sct",
}
DRUG_SYNONYMS: dict[str, str] = {}
NON_PARSEABLE_MARKERS = ["investigational"]
CART_MARKER = "cart"
HARD_FLAGS = {"non_parseable", "cart_therapy"}  # rules genuinely cannot adjudicate

FLAG_TEXT = {
    "partial_removal": "a drug was dropped without adding one (algorithm keeps it the same line; vendor/reviewer may split it)",
    "duplicate_line":  "a regimen repeats back-to-back (data artefact)",
    "non_parseable":   "contains an investigational / unnamed regimen the rules can't read",
    "cart_therapy":    "contains CAR-T therapy the rules can't adjudicate",
}


def normalize_drug(drug: str) -> str:
    d = drug.strip().lower()
    return DRUG_SYNONYMS.get(d, d)


def parse_cota_line(line_str: str) -> frozenset[str]:
    """Extract normalized drug names from one line entry; drop SCT tokens."""
    if not line_str or not line_str.strip():
        return frozenset()
    groups = re.findall(r"\[([^\]]+)\]", line_str)
    if not groups:
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
    return any(m in drug for drug in drug_set for m in NON_PARSEABLE_MARKERS)


def contains_cart(drug_set: frozenset[str]) -> bool:
    return any(CART_MARKER in drug for drug in drug_set)


def lot_algorithm_cota(trt_sequence: list[frozenset]) -> tuple[int, list[str]]:
    """Count lines of therapy (Rajkumar rules) + return review flags."""
    flags: list[str] = []
    for i in range(1, len(trt_sequence)):
        if trt_sequence[i] == trt_sequence[i - 1] and trt_sequence[i]:
            flags.append("duplicate_line")
            break

    deduped: list[frozenset] = []
    for drugs in trt_sequence:
        if not deduped or drugs != deduped[-1]:
            deduped.append(drugs)
    cleaned = [d for d in deduped if d]
    if not cleaned:
        return 0, flags

    for drugs in cleaned:
        if is_non_parseable(drugs) and "non_parseable" not in flags:
            flags.append("non_parseable")
        if contains_cart(drugs) and "cart_therapy" not in flags:
            flags.append("cart_therapy")

    lines = 1
    prev = cleaned[0]
    for i in range(1, len(cleaned)):
        curr = cleaned[i]
        overlap = curr & prev
        added = curr - prev
        if overlap:
            if added:
                lines += 1
            elif "partial_removal" not in flags:
                flags.append("partial_removal")
        else:
            lines += 1
        prev = curr
    return lines, flags


def parse_regimen_arg(raw: str) -> tuple[list[frozenset], list[str]]:
    """Split a ' | '-separated patient string into (drug-sets, display-strings)."""
    line_strings = [seg.strip() for seg in raw.split("|") if seg.strip()]
    seq = [parse_cota_line(seg) for seg in line_strings]
    return seq, line_strings


# ---------------------------------------------------------------------------
# Confidence  (calibrated buckets read from calibration.json)
# ---------------------------------------------------------------------------

def load_calibration() -> dict:
    if not os.path.exists(CALIB_PATH):
        sys.exit(
            f"ERROR: {CALIB_PATH} not found.\n"
            "Generate it once with:\n"
            "  .venv/bin/python src/py/calibrate_confidence.py"
        )
    with open(CALIB_PATH) as fh:
        return json.load(fh)


DECISION = {
    "HIGH":   "AUTO-ACCEPT  (no doctor needed)",
    "MEDIUM": "RECOMMEND REVIEW",
    "LOW":    "NEEDS DOCTOR REVIEW",
}
ICON = {"HIGH": "[OK ]", "MEDIUM": "[ ? ]", "LOW": "[ ! ]"}


def _bucket(c: dict, group: str, key: str) -> dict:
    """The calibrated stats dict for one bucket (accuracy, n, ci_low, ci_high, within1)."""
    return dict(c[group].get(key, {}))


def confidence_from_flags(flags: list[str], calib: dict) -> dict:
    """Ad-hoc input (no second count): confidence from the flag profile."""
    c = calib["cota"]
    if any(f in HARD_FLAGS for f in flags):
        return {"bucket": "LOW", **_bucket(c, "flag_buckets", "hardflag")}
    if flags:
        return {"bucket": "MEDIUM", **_bucket(c, "flag_buckets", "softflag")}
    return {"bucket": "HIGH", **_bucket(c, "flag_buckets", "noflag")}


def confidence_from_agreement(algo: int, vendor: int, flags: list[str], calib: dict) -> dict:
    """Independent second count available: stronger confidence from agreement."""
    c = calib["cota"]
    if algo != vendor:
        return {"bucket": "LOW", **_bucket(c, "agreement_buckets", "LOW")}
    if flags:
        return {"bucket": "MEDIUM", **_bucket(c, "agreement_buckets", "MEDIUM")}
    return {"bucket": "HIGH", **_bucket(c, "agreement_buckets", "HIGH")}


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _pct(x) -> str:
    return "n/a" if x is None else f"{x:.0%}"


def render(line_strings, count, flags, conf, *, vendor=None, mode="flags"):
    bucket = conf["bucket"]
    print("LOT Adjudication Assistant")
    print("-" * 58)
    print(f"Input: {len(line_strings)} regimen line(s)")
    for i, s in enumerate(line_strings, 1):
        print(f"  {i}. {s}")
    print()
    print(f"  Predicted LOT count : {count}")
    print(f"  Confidence          : {bucket} ({_pct(conf['accuracy'])})   {ICON[bucket]} -> {DECISION[bucket]}")
    print()
    print("  Why:")
    if mode == "agreement":
        if vendor is not None:
            verb = "agrees with" if count == vendor else "DISAGREES with"
            print(f"   - algorithm {verb} the COTA vendor count ({count} vs {vendor})")
    else:
        print("   - confidence is from the regimen's own flags - the realistic mode for a")
        print("     new patient, who usually has no independent second count")
    if flags:
        for f in flags:
            print(f"   - flag: {FLAG_TEXT.get(f, f)}")
    else:
        print("   - no edge-case flags fired")
    ci = ""
    if conf.get("ci_low") is not None:
        ci = f", 95% CI {conf['ci_low']:.0%}-{conf['ci_high']:.0%}"
    print(f"   - on the labelled data, cases like this matched human reviewers "
          f"{_pct(conf['accuracy'])} of the time{ci} (n={conf['n']})")
    if conf.get("within1") is not None:
        print(f"   - and landed within +/-1 line {conf['within1']:.0%} of the time "
              f"(an off-by-one count is often clinically close)")
    print()


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_predict(raw: str, vendor, calib: dict) -> None:
    seq, line_strings = parse_regimen_arg(raw)
    if not seq:
        sys.exit("No regimens parsed. Example:\n"
                 '  --regimens "[bortezomib, dexamethasone] | [lenalidomide]"')
    count, flags = lot_algorithm_cota(seq)
    if vendor is not None:
        conf = confidence_from_agreement(count, vendor, flags, calib)
        render(line_strings, count, flags, conf, vendor=vendor, mode="agreement")
    else:
        conf = confidence_from_flags(flags, calib)
        render(line_strings, count, flags, conf, mode="flags")


# ---------------------------------------------------------------------------
# Interactive mode  (array-form input)
# ---------------------------------------------------------------------------

def _to_set(drugs) -> frozenset:
    """Normalize a list of drug names into a frozenset (drop SCT tokens / blanks)."""
    out = set()
    for d in drugs:
        nd = normalize_drug(str(d))
        if nd and nd not in SCT_TOKENS:
            out.add(nd)
    return frozenset(out)


def _fmt(fs: frozenset) -> str:
    return "[" + ", ".join(sorted(fs)) + "]" if fs else "[]"


def parse_array_input(raw: str):
    """
    Parse a regimen sequence typed in array form. Each inner group = one line of therapy.

    Accepts any of:
        [[bortezomib, dexamethasone], [lenalidomide]]
        [["bortezomib","dexamethasone"], ["lenalidomide"]]
        [bortezomib, dexamethasone] | [lenalidomide]
        bortezomib, dexamethasone            (a single line)

    Returns (sequence_of_drugsets, display_strings).
    """
    raw = raw.strip()

    # 1) A quoted Python/JSON list (list-of-lists, or list-of-line-strings).
    try:
        obj = ast.literal_eval(raw)
    except Exception:
        obj = None
    if isinstance(obj, (list, tuple)) and obj:
        seq = []
        for item in obj:
            if isinstance(item, (list, tuple)):
                seq.append(_to_set(item))
            else:
                seq.append(_to_set(str(item).split(",")))
        return seq, [_fmt(s) for s in seq]

    # 2) Bracket groups: each [...] is one line.
    groups = re.findall(r"\[([^\[\]]+)\]", raw)
    if groups:
        seq = [_to_set(g.split(",")) for g in groups]
        return seq, [_fmt(s) for s in seq]

    # 3) Pipe-separated lines.
    if "|" in raw:
        seq = [_to_set(seg.split(",")) for seg in raw.split("|") if seg.strip()]
        return seq, [_fmt(s) for s in seq]

    # 4) A single comma-separated line.
    s = _to_set(raw.split(","))
    return ([s] if s else []), ([_fmt(s)] if s else [])


def cmd_interactive(calib: dict) -> None:
    print("LOT Adjudication Assistant  -  interactive mode")
    print("-" * 58)
    print("Enter a regimen sequence as an array. Each inner group = one line of therapy.")
    print("Examples:")
    print("  [[bortezomib, dexamethasone, lenalidomide], [carfilzomib, dexamethasone, daratumumab]]")
    print("  [bortezomib, dexamethasone] | [lenalidomide]")
    print()
    print("After each entry you'll be asked (optionally) for a 'second-opinion count':")
    print("an independent LOT number for the SAME patient from another source - e.g. the")
    print("count the COTA data vendor ships alongside their data. If you have one, it gives")
    print("a stronger confidence; if not (the usual case), just press Enter to skip.")
    print()
    print("Press Enter on an empty regimens line (or type q) to quit.")
    print()
    while True:
        try:
            raw = input("regimens> ").strip()
        except EOFError:
            break
        if not raw or raw.lower() in {"q", "quit", "exit"}:
            print("bye")
            break
        seq, display = parse_array_input(raw)
        if not seq:
            print("  could not parse - try  [[druga, drugb], [drugc]]\n")
            continue
        count, flags = lot_algorithm_cota(seq)
        try:
            v = input("Second-opinion LOT count? (optional - e.g. COTA's own number; Enter to skip)> ").strip()
        except EOFError:
            v = ""
        vendor = int(v) if v.lstrip("-").isdigit() else None
        print()
        if vendor is not None:
            conf = confidence_from_agreement(count, vendor, flags, calib)
            render(display, count, flags, conf, vendor=vendor, mode="agreement")
        else:
            conf = confidence_from_flags(flags, calib)
            render(display, count, flags, conf, mode="flags")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Predict Lines of Therapy (LOT) + calibrated confidence + auto-accept/needs-doctor decision.")
    ap.add_argument("--regimens", help='Patient regimen sequence, lines separated by " | ".')
    ap.add_argument("--vendor", type=int, default=None,
                    help="(optional) COTA vendor LOT count for this patient -> enables agreement-based confidence.")
    ap.add_argument("-i", "--interactive", action="store_true",
                    help="Interactive mode: type a regimen array, get a result. (Default when run with no args.)")
    args = ap.parse_args()

    calib = load_calibration()

    if args.interactive:
        cmd_interactive(calib)
    elif args.regimens:
        cmd_predict(args.regimens, args.vendor, calib)
    elif not sys.stdin.isatty():
        raw = sys.stdin.read().strip()
        if raw:
            cmd_predict(raw, args.vendor, calib)
        else:
            cmd_interactive(calib)
    else:
        cmd_interactive(calib)


if __name__ == "__main__":
    main()
