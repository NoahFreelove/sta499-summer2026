"""
Calibrate confidence for the LOT Adjudication Assistant.
=======================================================

Idea (matches the professor-approved plan):
  1. Predict LOT with the rule-based algorithm.
  2. Check whether it AGREES with the COTA data-vendor count.
  3. Measure, on the labelled data, how often each "agreement bucket" actually
     matches the human reviewer (the ground truth). That measured match-rate IS
     the calibrated confidence.
  4. Auto-accept the high-confidence bucket; route the rest to a doctor.

This script reuses the existing algorithm/loaders unchanged and writes
`calibration.json` (consumed by lot_cli.py) plus a human-readable summary.

Run with a python that has pandas + openpyxl, e.g.:
    .venv/bin/python src/py/calibrate_confidence.py
"""

from __future__ import annotations

import json
import os
import sys

import pandas as pd

# Import the proven algorithm + loaders from the existing scripts (same folder).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from textbook_algo import load_patient_data as load_flatiron, lot_algorithm_fixed
from textbook_algo_cota import load_patient_data as load_cota, lot_algorithm_cota

HARD_FLAGS = {"non_parseable", "cart_therapy"}  # rules genuinely cannot adjudicate these

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data", "LoT Adjudication Datasets.xlsx")
OUT = os.path.join(HERE, "calibration.json")


def _wilson(k: int, n: int, z: float = 1.96) -> tuple:
    """95% Wilson score interval for a binomial proportion -> (low, high)."""
    if n == 0:
        return (None, None)
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def _stats(correct: pd.Series, within1: "pd.Series | None" = None) -> dict:
    """{n, accuracy, ci_low, ci_high[, within1]} for a boolean 'is correct' series."""
    n = int(len(correct))
    k = int(correct.sum())
    acc = (k / n) if n else None
    lo, hi = _wilson(k, n)
    out = {"n": n, "accuracy": acc, "ci_low": lo, "ci_high": hi}
    if within1 is not None:
        out["within1"] = float(within1.mean()) if n else None
    return out


def calibrate_cota() -> dict:
    df = load_cota(DATA)
    res = df["trt_sequence"].apply(lot_algorithm_cota)
    df["lot_alg"] = res.apply(lambda x: x[0])
    df["flags"] = res.apply(lambda x: x[1])

    for c in ("rev_alpesh", "rev_alberto", "cota_lot"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Ground truth = primary reviewer (Alpesh), same convention as the existing script.
    df = df[df["rev_alpesh"].notna()].copy()

    df["correct"] = df["lot_alg"] == df["rev_alpesh"]
    df["within1"] = (df["lot_alg"] - df["rev_alpesh"]).abs() <= 1   # off-by-one accuracy
    df["has_flag"] = df["flags"].apply(lambda f: len(f) > 0)
    df["hard_flag"] = df["flags"].apply(lambda f: any(x in HARD_FLAGS for x in f))

    n = len(df)

    # --- Flag buckets: the DEPLOYABLE signal. ----------------------------------
    # Computed from the regimen sequence alone (no second count), so they apply
    # to a brand-new patient. This is what the app uses by default.
    noflag = df[~df["has_flag"]]
    softflag = df[df["has_flag"] & ~df["hard_flag"]]
    hardflag = df[df["hard_flag"]]

    # --- Agreement buckets: VALIDATION ONLY. -----------------------------------
    # These need an independent vendor count, which a real new patient won't
    # have. Restrict to the rows that actually carry one so re-calibrating on a
    # new cohort can't silently dump vendor-less patients into "disagree".
    has_vendor = df["cota_lot"].notna()
    df["agree_vendor"] = df["lot_alg"] == df["cota_lot"]   # NaN vendor -> False
    av = df[has_vendor]
    n_vendor = int(len(av))
    high = av[av["agree_vendor"] & ~av["has_flag"]]
    med = av[av["agree_vendor"] & av["has_flag"]]
    low = av[~av["agree_vendor"]]

    # Human-vs-human agreement = the "even experts disagree" floor (confidence ceiling).
    both = df[df["rev_alberto"].notna()]
    human_agreement = (
        float((both["rev_alpesh"] == both["rev_alberto"]).mean()) if len(both) else None
    )

    auto_cov = len(high) / n_vendor if n_vendor else None
    auto_acc = float(high["correct"].mean()) if len(high) else None

    return {
        "n": n,
        "ground_truth": "Alpesh (primary reviewer)",
        "overall_accuracy_algo": float(df["correct"].mean()),
        "overall_within1_algo": float(df["within1"].mean()),
        "overall_ci_algo": list(_wilson(int(df["correct"].sum()), n)),
        "overall_accuracy_vendor": (float((av["cota_lot"] == av["rev_alpesh"]).mean()) if n_vendor else None),
        "n_with_vendor": n_vendor,
        "human_agreement_alpesh_vs_alberto": human_agreement,
        "human_pair_n": int(len(both)),
        "agreement_buckets": {
            "HIGH": {**_stats(high["correct"], high["within1"]), "rule": "algorithm == vendor AND no edge-case flags"},
            "MEDIUM": {**_stats(med["correct"], med["within1"]), "rule": "algorithm == vendor BUT an edge-case flag fired"},
            "LOW": {**_stats(low["correct"], low["within1"]), "rule": "algorithm != vendor"},
        },
        "flag_buckets": {
            "noflag": {**_stats(noflag["correct"], noflag["within1"]), "rule": "no edge-case flags"},
            "softflag": {**_stats(softflag["correct"], softflag["within1"]), "rule": "soft flag (e.g. partial_removal / duplicate_line)"},
            "hardflag": {**_stats(hardflag["correct"], hardflag["within1"]), "rule": "hard flag (CAR-T / investigational) - rules cannot adjudicate"},
        },
        # Deployment triage = the no-flag auto-accept bucket. Uses ONLY the
        # flag signal, so it is the headline a new patient actually gets.
        "deployment_triage": {
            "auto_accept_coverage": (len(noflag) / n) if n else None,
            "auto_accept_accuracy": (float(noflag["correct"].mean()) if len(noflag) else None),
            "auto_accept_within1": (float(noflag["within1"].mean()) if len(noflag) else None),
            "review_coverage": (1 - len(noflag) / n) if n else None,
        },
        # Agreement triage = validation only (needs a vendor count to exist).
        "triage": {
            "auto_accept_coverage": auto_cov,
            "auto_accept_accuracy": auto_acc,
            "auto_accept_within1": float(high["within1"].mean()) if len(high) else None,
            "review_coverage": (1 - auto_cov) if auto_cov is not None else None,
        },
    }


def calibrate_flatiron() -> dict:
    """Secondary check on the Flatiron cohort (fixed algorithm, no flag system)."""
    df = load_flatiron(DATA)
    df["lot_alg"] = df["trt"].apply(lot_algorithm_fixed)
    for c in ("rev1_lot", "flatiron_lot"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[df["rev1_lot"].notna()].copy()

    df["correct"] = df["lot_alg"] == df["rev1_lot"]
    df["within1"] = (df["lot_alg"] - df["rev1_lot"]).abs() <= 1
    df["agree_vendor"] = df["lot_alg"] == df["flatiron_lot"]

    n = len(df)
    agree = df[df["agree_vendor"]]
    disagree = df[~df["agree_vendor"]]
    return {
        "n": n,
        "ground_truth": "Reviewer 1",
        "overall_accuracy_algo": float(df["correct"].mean()),
        "overall_within1_algo": float(df["within1"].mean()),
        "overall_accuracy_vendor": float((df["flatiron_lot"] == df["rev1_lot"]).mean()),
        "agreement_buckets": {
            "AGREE": {**_stats(agree["correct"], agree["within1"]), "rule": "algorithm == vendor (Flatiron line number)"},
            "DISAGREE": {**_stats(disagree["correct"], disagree["within1"]), "rule": "algorithm != vendor"},
        },
        "triage": {
            "auto_accept_coverage": (len(agree) / n) if n else None,
            "auto_accept_accuracy": float(agree["correct"].mean()) if len(agree) else None,
            "auto_accept_within1": float(agree["within1"].mean()) if len(agree) else None,
        },
    }


def _pct(x):
    return "  n/a" if x is None else f"{x:6.1%}"


def _ci(b):
    lo, hi = b.get("ci_low"), b.get("ci_high")
    return "" if lo is None else f"  95%CI[{lo:.0%}-{hi:.0%}]"


def _w1(b):
    w = b.get("within1")
    return "" if w is None else f"  +/-1={w:.0%}"


def main() -> None:
    cal = {"cota": calibrate_cota(), "flatiron": calibrate_flatiron()}
    with open(OUT, "w") as fh:
        json.dump(cal, fh, indent=2)

    c = cal["cota"]
    print("=" * 66)
    print("CONFIDENCE CALIBRATION  -  COTA cohort  (truth = Alpesh)")
    print("=" * 66)
    ci = c["overall_ci_algo"]
    print(f"  Patients ...................... {c['n']}")
    print(f"  Algorithm overall accuracy .... {_pct(c['overall_accuracy_algo'])}"
          f"  (95% CI {ci[0]:.0%}-{ci[1]:.0%};  within +/-1 = {c['overall_within1_algo']:.0%})")
    print(f"  Vendor   overall accuracy ..... {_pct(c['overall_accuracy_vendor'])}")
    print(f"  Human vs human agreement ...... {_pct(c['human_agreement_alpesh_vs_alberto'])}"
          f"   (n={c['human_pair_n']})  <- confidence ceiling")
    print()
    print("  AGREEMENT BUCKETS  (independent second count available)")
    for name, b in c["agreement_buckets"].items():
        print(f"    {name:7s} acc={_pct(b['accuracy'])}{_ci(b)}{_w1(b)}  n={b['n']:<4d}  [{b['rule']}]")
    print()
    print("  FLAG BUCKETS  (ad-hoc input, no second count)")
    for name, b in c["flag_buckets"].items():
        print(f"    {name:9s} acc={_pct(b['accuracy'])}{_ci(b)}{_w1(b)}  n={b['n']:<4d}  [{b['rule']}]")
    print()
    dt = c["deployment_triage"]
    print("  TRIAGE HEADLINE  (deployment - no second opinion needed)")
    print(f"    Auto-accept the no-flag cases . {_pct(dt['auto_accept_coverage'])} of patients")
    print(f"    ...and on those, accuracy is .. {_pct(dt['auto_accept_accuracy'])}"
          f"  (within +/-1 = {_pct(dt['auto_accept_within1'])})")
    print(f"    Routed to a doctor ............ {_pct(dt['review_coverage'])} of patients")
    t = c["triage"]
    print(f"    [validation] given a vendor 2nd count, agreement reproduces this split:"
          f" {_pct(t['auto_accept_coverage'])} auto-accept @ {_pct(t['auto_accept_accuracy'])}")
    print()

    f = cal["flatiron"]
    print("=" * 66)
    print("CROSS-CHECK  -  Flatiron cohort  (truth = Reviewer 1)")
    print("=" * 66)
    print(f"  Patients={f['n']}  algo={_pct(f['overall_accuracy_algo'])} (within +/-1 = {f['overall_within1_algo']:.0%})  vendor={_pct(f['overall_accuracy_vendor'])}")
    for name, b in f["agreement_buckets"].items():
        print(f"    {name:9s} acc={_pct(b['accuracy'])}{_ci(b)}{_w1(b)}  n={b['n']:<4d}  [{b['rule']}]")
    print()
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
