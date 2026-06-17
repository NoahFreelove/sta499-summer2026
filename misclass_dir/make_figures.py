"""
Reproduce the 5 recommended visualizations from the COTA misclassification
pattern hunt (FINDINGS_2026-06-17.md). Read-only on the workbook; writes PNGs
to misclass_dir/figs/.

Run: /Users/freelove/STA499/.venv/bin/python make_figures.py
"""
from __future__ import annotations
from pathlib import Path
import re, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import fisher_exact

BASE = Path(__file__).resolve().parent
FIG = BASE / "figs"; FIG.mkdir(exist_ok=True)
DATA = BASE / "COTA_misclassified_rows_UPD.xlsx"

OVER_C = "#d1495b"   # over-split (error)
ALIGN_C = "#3a6ea5"  # aligned (control)
ACCENT = "#2a9d8f"
GREY = "#9aa0a6"

# ---------------------------------------------------------------- load
a = pd.read_excel(DATA, sheet_name="All_COTA_With_Transitions")
a.columns = a.columns.str.strip()
a = a.sort_values(["cpid", "_original_row_order"]).copy()
S = a["transition_alignment_status"]
a["is_over"] = (S == "candidate_misclassification_cota_over_split")
TR = a[S != "baseline_first_lot_row"].copy()   # 780 transitions
N_OVER = int(a["is_over"].sum())

# ---------------------------------------------------------------- helpers
def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return p, max(0, c-h), min(1, c+h)

def or_ci(a_, b_, c_, d_):
    if min(a_, b_, c_, d_) == 0:
        a_, b_, c_, d_ = a_+0.5, b_+0.5, c_+0.5, d_+0.5
    orr = (a_*d_)/(b_*c_)
    se = math.sqrt(1/a_+1/b_+1/c_+1/d_)
    return orr, math.exp(math.log(orr)-1.96*se), math.exp(math.log(orr)+1.96*se)

def blank(x):
    return pd.isna(x) or str(x).strip() == ""

def drugset(name):
    if pd.isna(name):
        return frozenset()
    toks = re.findall(r"[a-zA-Z][a-zA-Z0-9\- ]+", str(name).lower())
    return frozenset(t.strip() for t in toks if t.strip())

def n_phases(name):
    if pd.isna(name):
        return 0
    g = re.findall(r"\[([^\]]+)\]", str(name))
    return len(g) if g else (1 if str(name).strip() else 0)

# ================================================================ FIG 1: decay hazard
def fig_decay():
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for line_col, color, lab, mark in [
        ("cota_lot_numeric", OVER_C, "COTA line number (as-reported)", "o"),
        ("doctor_lot_numeric_for_transition", ACCENT, "Doctor ground-truth line", "s"),
    ]:
        rows = []
        for k in range(2, 11):
            sub = TR[TR[line_col] == k]
            if len(sub) < 1:
                continue
            p, lo, hi = wilson(int(sub["is_over"].sum()), len(sub))
            rows.append((k, p, lo, hi, len(sub)))
        ks = [r[0] for r in rows]
        ps = [r[1]*100 for r in rows]
        lo = [(r[1]-r[2])*100 for r in rows]
        hi = [(r[3]-r[1])*100 for r in rows]
        ax.errorbar(ks, ps, yerr=[lo, hi], marker=mark, color=color, lw=2,
                    capsize=3, label=lab, alpha=0.9)
    ax.set_xlabel("Line of therapy")
    ax.set_ylabel("COTA over-split rate (%)  ± 95% Wilson CI")
    ax.set_title("Finding 7: the decay is a real monotone hazard, not survivorship\n"
                 "(steeper on the doctor's line → not an x-axis artifact; late lines L11+ omitted = small-n)")
    ax.set_xticks(range(2, 11))
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    ax.annotate("peak ≈14% at line 2", xy=(2, 14), xytext=(3.2, 17),
                fontsize=9, color=GREY,
                arrowprops=dict(arrowstyle="->", color=GREY))
    fig.tight_layout()
    fig.savefig(FIG/"fig1_decay_hazard.png", dpi=150)
    plt.close(fig)

# ================================================================ FIG 2: discontinue reason flip
def fig_reason():
    a["blank_cur"] = a["discontinue_reason"].apply(blank)
    a["blank_prev"] = a.groupby("cpid")["discontinue_reason"].shift(1).apply(blank)
    T = a[S != "baseline_first_lot_row"]
    o, al = T[T.is_over], T[~T.is_over]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 5.5))

    # Panel A: the flip
    attribs = [("Previous-row\n(line that ended)", "blank_prev"),
               ("Current-row\n(line that began)", "blank_cur")]
    xs = np.arange(len(attribs)); w = 0.36
    over_v, align_v, labels = [], [], []
    for name, col in attribs:
        ov = o[col].mean()*100; av = al[col].mean()*100
        over_v.append(ov); align_v.append(av)
        a11 = int(o[col].sum()); a01 = int(al[col].sum())
        orr, p = fisher_exact([[a11, len(o)-a11], [a01, len(al)-a01]])
        labels.append(f"OR={orr:.2f}\np={p:.1e}")
    axL.bar(xs-w/2, over_v, w, color=OVER_C, label=f"Over-split (n={len(o)})")
    axL.bar(xs+w/2, align_v, w, color=ALIGN_C, label=f"Aligned (n={len(al)})")
    for i, t in enumerate(labels):
        axL.text(i, max(over_v[i], align_v[i])+1.5, t, ha="center", fontsize=9,
                 color=("black" if "e-0" not in t else OVER_C))
    axL.set_xticks(xs); axL.set_xticklabels([n for n, _ in attribs])
    axL.set_ylabel("% with blank / undocumented discontinue reason")
    axL.set_title("Finding 1: blank reason on the PREVIOUS line\ntriples over-split odds — but only with the right row")
    axL.set_ylim(0, 42); axL.legend(frameon=False); axL.grid(axis="y", alpha=0.25)

    # Panel B: prev-row reason composition
    def bucket(x):
        if blank(x):
            return "blank/undocumented"
        t = str(x).split("\\")[0].strip().lower()
        return "progression" if t == "progression" else "toxicity/other"
    T2 = a[S != "baseline_first_lot_row"].copy()
    T2["prev_reason"] = a.groupby("cpid")["discontinue_reason"].shift(1).reindex(T2.index)
    T2["pb"] = T2["prev_reason"].apply(bucket)
    cats = ["blank/undocumented", "progression", "toxicity/other"]
    xs2 = np.arange(len(cats))
    ov = [ (T2[T2.is_over].pb == c).mean()*100 for c in cats ]
    av = [ (T2[~T2.is_over].pb == c).mean()*100 for c in cats ]
    axR.bar(xs2-w/2, ov, w, color=OVER_C, label="Over-split")
    axR.bar(xs2+w/2, av, w, color=ALIGN_C, label="Aligned")
    axR.set_xticks(xs2); axR.set_xticklabels(cats, rotation=12)
    axR.set_ylabel("% of transitions (previous-line reason)")
    axR.set_title("Previous-line reason composition\n(over-splits are enriched for undocumented endings)")
    axR.legend(frameon=False); axR.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG/"fig2_discontinue_reason_flip.png", dpi=150)
    plt.close(fig)

# ================================================================ FIG 3: family forest
def fig_family():
    fam = "reconstructed_family_combination"
    T = a[S != "baseline_first_lot_row"].copy()
    T[fam] = (T[fam].astype(str)
              .str.replace(r"\s*\(Not in provided Fiona's category list\)", "", regex=True)
              .str.strip())
    rows = []
    for combo, sub in T.groupby(fam):
        n = len(sub)
        if n < 5 or combo in ("", "nan"):
            continue
        a_ = int(sub.is_over.sum()); b_ = n - a_
        rest = T[T[fam] != combo]
        c_ = int(rest.is_over.sum()); d_ = len(rest) - c_
        orr, lo, hi = or_ci(a_, b_, c_, d_)
        _, p = fisher_exact([[a_, b_], [c_, d_]])
        rows.append((combo, orr, lo, hi, a_, n, p))
    tested = len(rows)
    bonf = 0.05 / tested
    rows.sort(key=lambda r: r[1])
    fig, ax = plt.subplots(figsize=(10.5, max(4, 0.42*len(rows)+2)))
    for i, (combo, orr, lo, hi, k, n, p) in enumerate(rows):
        is_pi = combo == "Proteasome inhibitors + Steroids"
        passes = p < bonf
        col = ACCENT if is_pi else (OVER_C if orr > 1 else ALIGN_C)
        ax.plot([lo, hi], [i, i], color=col, lw=2.2 if passes else 1.4,
                alpha=0.95 if (is_pi or passes) else 0.6)
        ax.plot(orr, i, "o", color=col, ms=9 if is_pi else 6,
                mec="black" if passes else col, mew=1.4 if passes else 0)
        ax.text(hi*1.05, i, f"{k}/{n}" + ("  ★Bonf" if passes else ""),
                va="center", fontsize=8, color="black" if passes else GREY)
    ax.axvline(1, color="black", lw=1)
    ax.set_xscale("log")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0][:46] for r in rows], fontsize=8)
    ax.set_xlabel("Over-split odds ratio vs all other transitions (log scale)  —  <1 = protected,  >1 = error-prone")
    ax.set_title(f"Finding 3: per-family over-split OR (current-row family, n≥5 opportunities; {tested} families)\n"
                 f"PI+Steroids (teal) is the only Bonferroni survivor (α={bonf:.4f}) — and it is PROTECTED")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG/"fig3_family_forest.png", dpi=150)
    plt.close(fig)

# ================================================================ FIG 4: identical-set waterfall
def fig_waterfall():
    a["ds"] = a["line_of_therapy_name"].apply(drugset)
    a["ds_prev"] = a.groupby("cpid")["ds"].shift(1)
    a["nph"] = a["line_of_therapy_name"].apply(n_phases)
    a["nph_prev"] = a.groupby("cpid")["nph"].shift(1)
    o = a[a.is_over]
    same = o[[len(c) > 0 and c == p for c, p in zip(o.ds, o.ds_prev)]]
    phase_only = same[same.nph != same.nph_prev]
    steps = [("All COTA\nover-splits", N_OVER, OVER_C),
             ("Identical drug set\nto prior row", len(same), "#e08a3c"),
             ("…differ only in\nbracket-phase count", len(phase_only), ACCENT)]
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for i, (lab, v, c) in enumerate(steps):
        ax.bar(i, v, 0.62, color=c)
        ax.text(i, v+0.8, f"{v}\n({v/N_OVER*100:.0f}% of all)", ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(range(len(steps)))
    ax.set_xticklabels([s[0] for s in steps])
    ax.set_ylabel("number of over-split rows")
    ax.set_ylim(0, N_OVER+8)
    ax.set_title("Finding 2: the dominant mechanism is a bracket-phase FORMATTING artifact\n"
                 "COTA opens a new line where the drug set is unchanged (OR ≈39 vs aligned)")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG/"fig4_identical_set_waterfall.png", dpi=150)
    plt.close(fig)
    return len(same), len(phase_only)

# ================================================================ FIG 5: severity = accumulation
def fig_severity():
    o = a[a.is_over].copy()
    per_pt = o.groupby("cpid").size()
    o = o.merge(per_pt.rename("pt_count"), left_on="cpid", right_index=True)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    # Panel A: cumulative shift values among the 74 rows, split by repeat-offender
    shifts = sorted(o["cota_minus_doctor_lot_shift"].unique())
    single = [((o.cota_minus_doctor_lot_shift == s) & (o.pt_count == 1)).sum() for s in shifts]
    repeat = [((o.cota_minus_doctor_lot_shift == s) & (o.pt_count > 1)).sum() for s in shifts]
    axL.bar(shifts, single, color=ALIGN_C, label="patient has 1 over-split")
    axL.bar(shifts, repeat, bottom=single, color=OVER_C, label="repeat-offender patient (>1)")
    axL.set_xlabel("cumulative cota−doctor shift at the over-split row")
    axL.set_ylabel("number of over-split rows")
    axL.set_title("Every per-transition jump is +1 (doc_inc=0).\nCumulative ≥2 only appears in repeat-offender patients")
    axL.set_xticks(shifts); axL.legend(frameon=False); axL.grid(axis="y", alpha=0.25)

    # Panel B: per-patient over-split count distribution
    dist = per_pt.value_counts().sort_index()
    axR.bar(dist.index, dist.values, color=ACCENT)
    for x, y in zip(dist.index, dist.values):
        axR.text(x, y+0.4, str(y), ha="center", fontsize=9)
    axR.set_xlabel("over-splits per patient")
    axR.set_ylabel("number of patients")
    axR.set_title(f"Finding 4: severity is opportunity, not super-spreaders\n"
                  f"{len(per_pt)} patients carry {N_OVER} over-splits (no overdispersion beyond exposure)")
    axR.set_xticks(sorted(per_pt.unique())); axR.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG/"fig5_severity_accumulation.png", dpi=150)
    plt.close(fig)

# ---------------------------------------------------------------- run
if __name__ == "__main__":
    fig_decay()
    fig_reason()
    fig_family()
    n_same, n_phase = fig_waterfall()
    fig_severity()
    print("figs written to", FIG)
    print(f"  fig4 sanity: identical-set={n_same}/{N_OVER}, bracket-phase-only={n_phase}")
    for p in sorted(FIG.glob("*.png")):
        print("  ", p.name)
