#!/usr/bin/env python3
"""Figure S8: re-validation against functionally screened epitope labels.

    python3 scripts/fig_s8_tcrvdb.py [out_dir]

Draws Supplementary Figure S8 from the summary tables that
scripts/tcrvdb_revalidation.py writes. The TCRvdb labels (Messemaker et al.,
Immunity 59:2775-2788, 2026) are available on request only and may not be
redistributed, so this script reads just the committed summaries:

  results/tables/extensions/tcrvdb_cluster_enrichment.csv   panels A, C, D
  results/tables/extensions/tcrvdb_purity_recomputed.csv    panel B
  results/tables/extensions/tcrvdb_enrichment_stratified.csv panels C, D

A  Retention: the share of receptors each setting places in a cluster. The
   first three settings are scored on the same 1,151 CDR3-beta keys; the
   1.25 A paired setting is keyed on the paired receptor (1,504 keys).
B  Cluster purity with labels as published and with every non-validating
   receptor reassigned to a class of its own.
C  SPACE-TCR 0.75 A: share of labels that validated, among receptors it
   clustered and those it did not, within strata of CDR3-beta sequence
   neighbours (Levenshtein distance <= 2).
D  The same comparison for ClusTCR.

Faded bars are strata where the difference is not significant (one-sided
Fisher, p >= 0.05, uncorrected). ClusTCR clusters too few receptors with no
sequence neighbour to test that stratum, and the panel says so.

The script checks the plotted values against the tables before drawing and
stops on any disagreement.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(REPO, "results", "tables", "extensions")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "results", "figures", "rebuild")
os.makedirs(OUT, exist_ok=True)

RED, PINK, DARK, GREY = "#e31a1c", "#fb9a99", "#6b6b6b", "#d9d9d9"
TXT, SUB = "#222222", "#888888"
if any("Inter" in f.name for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#999999", "axes.titleweight": "bold",
                     "axes.titlesize": 13, "axes.titlelocation": "left",
                     "axes.labelcolor": "#444444", "xtick.color": "#444444",
                     "ytick.color": "#444444", "hatch.linewidth": 1.2})

enr = pd.read_csv(os.path.join(TAB, "tcrvdb_cluster_enrichment.csv"))
pur = pd.read_csv(os.path.join(TAB, "tcrvdb_purity_recomputed.csv"))
strat = pd.read_csv(os.path.join(TAB, "tcrvdb_enrichment_stratified.csv"),
                    dtype={"stratum": str})

failures = []
def check(label, got, want, tol):
    ok = abs(got - want) <= tol
    print(f"    {label:<52} got {got:.4f}  want {want:.4f}  {'OK' if ok else 'MISMATCH'}")
    if not ok:
        failures.append(label)

# ---- panel A: retention ------------------------------------------------------
SETTINGS = [("SPACE-TCR RMSD (0.75A)", "focal_0.75A_cdr3b", "SPACE-TCR\n0.75 Å", RED),
            ("SPACE-TCR RMSD (1A)", "focal_0.75A_cdr3b", "SPACE-TCR\n1.0 Å", PINK),
            ("clusTCR", "focal_0.75A_cdr3b", "ClusTCR\nsequence", DARK),
            ("SPACE-TCR RMSD (1.25A)", "permissive_1.25A_paired", "SPACE-TCR\n1.25 Å paired", RED)]

ret = []
for m, an, lab, col in SETTINGS:
    r = enr[(enr.method == m) & (enr.analysis == an)].iloc[0]
    k = int(r.clustered_validating + r.clustered_nonvalidating)
    n = k + int(r.unclustered_validating + r.unclustered_nonvalidating)
    ret.append((lab, col, k, n))

print("Figure S8 - re-validation against functionally screened labels")
check("A focal retention 238/1151", ret[0][2] / ret[0][3], 238 / 1151, 1e-9)
check("A ClusTCR retention 118/1151", ret[2][2] / ret[2][3], 118 / 1151, 1e-9)
check("A 1.25 A paired retention 698/1504", ret[3][2] / ret[3][3], 698 / 1504, 1e-9)
if len({r[3] for r in ret[:3]}) != 1:
    failures.append("A: the first three settings do not share one key universe")

# ---- panel B: purity ---------------------------------------------------------
def purity(m, an, basis):
    return float(pur[(pur.method == m) & (pur.analysis == an) & (pur.basis == basis)].purity.iloc[0])

pb = [(lab, col, purity(m, an, "as_published"), purity(m, an, "non-validating_relabelled"))
      for m, an, lab, col in SETTINGS]
check("B focal purity, as published", pb[0][2], 0.9963, 0.001)
check("B focal purity, non-validating relabelled", pb[0][3], 0.9572, 0.001)
check("B 1.25 A paired purity, relabelled", pb[3][3], 0.8775, 0.001)

# ---- panels C, D: stratified validation --------------------------------------
STRATA = ["0", "1-2", "3-9", "10+"]
ld = strat[strat.control == "LD<=2 sequence degree"]
def rows(m):
    d = ld[ld.method == m].set_index("stratum")
    return {s: d.loc[s] for s in STRATA if s in d.index}

sp, cl = rows("SPACE-TCR RMSD (0.75A)"), rows("clusTCR")
check("C SPACE-TCR stratum 0, clustered", sp["0"].rate_clustered, 0.7895, 0.0005)
check("C SPACE-TCR stratum 0, unclustered", sp["0"].rate_unclustered, 0.2798, 0.0005)
check("D ClusTCR stratum 3-9 odds ratio", cl["3-9"].odds_ratio, 0.9358, 0.0005)
if "0" in cl:
    failures.append("D: ClusTCR now has a degree-0 stratum; the 'too few' label is wrong")

if failures:
    sys.exit("Error: figure value disagrees with the tables: " + "; ".join(failures))

# ---- draw ----------------------------------------------------------------------
fig = plt.figure(figsize=(14, 10.3))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.18, left=0.06, right=0.98,
                      top=0.95, bottom=0.12)

ax = fig.add_subplot(gs[0, 0])
x = range(len(ret))
ax.bar(x, [k / n for _, _, k, n in ret], width=0.55, color=[c for _, c, _, _ in ret])
for i, (lab, col, k, n) in enumerate(ret):
    ax.text(i, k / n + 0.01, f"{100 * k / n:.1f}%", ha="center", va="bottom",
            fontsize=12, fontweight="bold", color=TXT)
    ax.text(i, 0.012, f"{k}/{n}", ha="center", va="bottom", fontsize=9,
            color="white" if col != PINK else "#ffffff")
ax.set_xticks(list(x), [r[0] for r in ret], fontsize=10)
ax.set_ylim(0, 0.55)
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
ax.set_ylabel("receptors placed in a cluster")
ax.set_title("A · Retention")
ax.text(1, 0.535, f"same {ret[0][3]:,} keys for the first three", ha="center",
        fontsize=9, color=SUB)
ax.grid(axis="y", color="#eeeeee"); ax.set_axisbelow(True)

ax = fig.add_subplot(gs[0, 1])
w = 0.36
for i, (lab, col, p0, p1) in enumerate(pb):
    ax.bar(i - w / 2, p0, w, color=col)
    ax.bar(i + w / 2, p1, w, color=col, alpha=0.75, hatch="////", edgecolor="white")
    ax.text(i - w / 2, p0 + 0.002, f"{100 * p0:.1f}", ha="center", va="bottom", fontsize=10)
    ax.text(i + w / 2, p1 + 0.002, f"{100 * p1:.1f}", ha="center", va="bottom", fontsize=10)
ax.set_xticks(range(len(pb)), [r[0] for r in pb], fontsize=10)
ax.set_ylim(0.85, 1.015)
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
ax.set_ylabel("cluster purity")
ax.set_title("B · Purity")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=DARK, label="labels as published"),
                   Patch(facecolor=DARK, alpha=0.75, hatch="////", edgecolor="white",
                         label="non-validating relabelled")],
          loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False, fontsize=9)
ax.grid(axis="y", color="#eeeeee"); ax.set_axisbelow(True)

def strata_panel(ax, d, colour, title, method):
    w = 0.36
    for i, s in enumerate(STRATA):
        if s not in d:
            ax.text(i, 0.04, "too few\nclustered\nto test", ha="center", va="bottom",
                    fontsize=8, color=SUB, style="italic")
            continue
        r = d[s]
        sig = r.p_one_sided < 0.05
        n_un = int(r.unclustered_validating + r.unclustered_nonvalidating)
        n_cl = int(r.clustered_validating + r.clustered_nonvalidating)
        ax.bar(i - w / 2, r.rate_unclustered, w, color=GREY, edgecolor="#bbbbbb")
        ax.bar(i + w / 2, r.rate_clustered, w, color=colour, alpha=1.0 if sig else 0.35,
               edgecolor="#bbbbbb" if not sig else colour)
        ax.text(i - w / 2, r.rate_unclustered + 0.01, f"{100 * r.rate_unclustered:.0f}%",
                ha="center", va="bottom", fontsize=10, color="#555555")
        ax.text(i + w / 2, r.rate_clustered + 0.01, f"{100 * r.rate_clustered:.0f}%",
                ha="center", va="bottom", fontsize=10, fontweight="bold" if sig else "normal",
                color=colour if sig else "#999999")
        ax.text(i - w / 2, 0.012, f"n={n_un}", ha="center", fontsize=8, color="#777777")
        ax.text(i + w / 2, 0.012, f"n={n_cl}", ha="center", fontsize=8,
                color="white" if sig else "#777777")
        if not sig:
            top = max(r.rate_clustered, r.rate_unclustered)
            ax.text(i, top + 0.06, "no real\ndifference", ha="center", va="bottom",
                    fontsize=8, color=SUB)
    ax.set_xticks(range(len(STRATA)), STRATA)
    ax.set_xlabel("CDR3β sequence neighbours at LD ≤ 2", color=SUB)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    ax.set_ylabel("share of epitope labels that\nheld up when tested")
    ax.set_title(title)
    ax.legend(handles=[Patch(facecolor=GREY, edgecolor="#bbbbbb", label="receptors it did NOT cluster"),
                       Patch(color=colour, label="receptors it clustered")],
              loc="upper left", frameon=False, fontsize=9)
    ax.grid(axis="y", color="#eeeeee"); ax.set_axisbelow(True)

strata_panel(fig.add_subplot(gs[1, 0]), sp, RED,
             "C · Labels held up more often when SPACE-TCR clustered the receptor", "SPACE-TCR")
strata_panel(fig.add_subplot(gs[1, 1]), cl, DARK,
             "D · ClusTCR shows no such effect on the same receptors", "ClusTCR")

sig_p = [sp[s].p_one_sided for s in STRATA if s in sp and sp[s].p_one_sided < 0.05]
if any(p * len(STRATA) >= 0.05 for p in sig_p):
    sys.exit("Error: a significant SPACE-TCR stratum does not survive Bonferroni; fix the footnote")
fig.text(0.01, 0.015,
         "Purity axis starts at 85%. One-sided Fisher, uncorrected; SPACE-TCR's significant "
         "strata survive Bonferroni (p = " + ", ".join(f"{p:.1e}" for p in sig_p) + ").",
         fontsize=9, color=SUB)

path = os.path.join(OUT, "figS8_tcrvdb_revalidation.png")
fig.savefig(path, dpi=300)
print(f"  saved: {os.path.basename(path)}")
