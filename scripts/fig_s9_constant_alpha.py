#!/usr/bin/env python3
"""Figure S9: paired alpha-beta versus single-chain beta clustering.

    python3 scripts/fig_s9_constant_alpha.py [out_dir]

One CDR-length bin, CDRB(10,8,13), holding 1,523 receptors, clustered two ways:

  paired      all six CDR loops of the real paired models, framework-aligned.
              Computed here from the deposited distance tables, so it needs no
              structures.
  single      the beta loops only, from the same receptors re-modelled with
              ONE constant alpha partner (tcr1's), so any difference between
              receptors' beta geometry comes from the beta sequence alone. Read
              from the SPACE2TCR distance tables of that run, binned on the beta
              loops (results/tables/extensions/s9_constant_alpha/).

Panel D reads the CDR3 dihedrals of the constant-alpha models (Stage D of the
same run) from results/tables/extensions/s9_constant_alpha_dihedrals.tsv.gz.
scripts/s9_constant_alpha.sh rebuilds the run from the Zenodo structures.

Clustering only compares receptors inside a CDR-length bin, so a whole bin
gives exactly the clusters a full-cohort run would give for those receptors.

A  Purity against retention.
B  Clusters at purity >= 0.9 against threshold.
C  Retention against threshold. Each arm is scored on its own geometry (six
   loops for paired, beta loops on the beta framework for single), so the
   thresholds in B and C are not on a common scale.
D  Pairwise CDR3-beta backbone dihedral distance (mean circular distance over
   phi/psi, IMGT 105-117), real-alpha against constant-alpha models. Dihedrals
   are internal coordinates and so chain-separable; if the spread is unchanged,
   the constant partner did not homogenise the beta conformations.

A constant alpha is an approximation to a single-chain predictor, which does
not exist in this toolchain (TCRBuilder2 cannot model one chain).

The script checks the plotted values against the published ones before drawing
and stops on any disagreement.
"""
import collections
import csv
import glob
import gzip
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = os.path.join(REPO, "results", "tables", "extensions")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "results", "figures", "rebuild")
os.makedirs(OUT, exist_ok=True)

IDS = os.path.join(EXT, "s9_bin_CDRB_10_8_13_ids.txt")
SINGLE = os.path.join(EXT, "s9_constant_alpha")
DIH_CONST = os.path.join(EXT, "s9_constant_alpha_dihedrals.tsv.gz")
DIH_REAL = os.path.join(REPO, "data", "cohort_structural", "dihedrals.tsv.gz")
DIST = os.path.join(REPO, "data", "cohort_structural", "distance_tables")
LABELS = os.path.join(REPO, "data", "jouannet_labelled_input.tsv")
THRESHOLDS = [0.4, 0.5, 0.6, 0.75, 0.9, 1.0]
FOCAL = 0.75

for f in (SINGLE, DIH_CONST):
    if not os.path.exists(f):
        sys.exit(f"missing {os.path.relpath(f, REPO)}\n"
                 "This is output of the constant-alpha run: the Stage C distance tables "
                 "and the Stage D dihedral table.\nRebuild it with "
                 "bash scripts/s9_constant_alpha.sh (needs the Zenodo structures and "
                 "ImmuneBuilder).")

RED, GREY = "#e31a1c", "#6b6b6b"
if any("Inter" in f.name for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#999999", "axes.titleweight": "bold",
                     "axes.titlesize": 13, "axes.titlelocation": "left",
                     "axes.labelcolor": "#444444", "xtick.color": "#444444",
                     "ytick.color": "#444444"})

ids = [l.strip() for l in open(IDS) if l.strip()]
idset = set(ids)
lab = {r["id"]: r["epitope"] for r in csv.DictReader(open(LABELS), delimiter="\t")}

# ---- distance tables -----------------------------------------------------------
def read_pairs(folder):
    out = []
    for f in sorted(glob.glob(os.path.join(folder, "*.pairwise.tsv.gz"))):
        with gzip.open(f, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                a = os.path.basename(r["id_i"]).rsplit(".pdb", 1)[0]
                b = os.path.basename(r["id_j"]).rsplit(".pdb", 1)[0]
                if a in idset and b in idset and r["rmsd"] not in ("", "nan", "NA"):
                    out.append((a, b, float(r["rmsd"])))
    return out

pairs = read_pairs(DIST)            # paired arm: the published cohort run
pairs_single = read_pairs(SINGLE)   # single arm: the constant-alpha run

# The constant-alpha run must cover exactly this bin.
members = set()
for f in glob.glob(os.path.join(SINGLE, "bin_*.csv")):
    for r in csv.DictReader(open(f)):
        members.add(os.path.basename(r["ID"]).rsplit(".pdb", 1)[0])
if members != idset:
    sys.exit(f"Error: the constant-alpha run holds {len(members)} receptors, "
             f"not the {len(idset)} of the bin")

def components(edges):
    parent = {}
    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x]); x = parent[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    groups = collections.defaultdict(list)
    for x in set(parent) | set(parent.values()): groups[find(x)].append(x)
    return [g for g in groups.values() if len(g) >= 2]

def score(clusters, n):
    hi = [collections.Counter(lab[x] for x in c).most_common(1)[0][1] / len(c) >= 0.9
          for c in clusters]
    return np.mean(hi) if hi else np.nan, sum(map(len, clusters)) / n

paired = pd.DataFrame(
    [(t, *score(components([(a, b) for a, b, d in pairs if d <= t]), len(ids)))
     for t in THRESHOLDS], columns=["threshold", "hi", "retention"])

single = pd.DataFrame(
    [(t, *score(components([(a, b) for a, b, d in pairs_single if d <= t]), len(ids)))
     for t in THRESHOLDS], columns=["threshold", "hi", "retention"])

# ---- panel D: pairwise CDR3-beta dihedral distance ---------------------------
def beta_matrix(path):
    d = pd.read_csv(path, sep="\t", dtype={"resi": str})
    d = d[(d.chain == "B") & d.tcr_id.isin(idset)].copy()
    d["r"] = pd.to_numeric(d.resi, errors="coerce")
    d = d[(d.r >= 105) & (d.r <= 117)].sort_values(["tcr_id", "r"])
    m = {k: np.ravel(np.column_stack([g.phi_rad.values, g.psi_rad.values]))
         for k, g in d.groupby("tcr_id")}
    return np.stack([m[k] for k in ids if k in m])

def pairwise(M):
    out = []
    for i in range(len(M) - 1):
        x = np.abs(M[i + 1:] - M[i]); x = np.minimum(x, 2 * np.pi - x)
        out.append(x.mean(1))
    return np.concatenate(out)

d_real, d_const = pairwise(beta_matrix(DIH_REAL)), pairwise(beta_matrix(DIH_CONST))

# ---- checks ------------------------------------------------------------------
failures = []
def check(label, got, want, tol):
    ok = abs(got - want) <= tol
    print(f"    {label:<48} got {got:.4f}  want {want:.4f}  {'OK' if ok else 'MISMATCH'}")
    if not ok: failures.append(label)

at = lambda df, col: float(df.loc[np.isclose(df.threshold, FOCAL), col].iloc[0])
print("Figure S9 - paired versus single-chain beta, bin CDRB(10,8,13)")
check("bin size", len(ids), 1523, 0)
check("paired, clusters at purity >= 0.9, 0.75 A", at(paired, "hi"), 0.964, 0.0005)
check("paired, retention, 0.75 A", at(paired, "retention"), 0.054, 0.001)
check("single, clusters at purity >= 0.9, 0.75 A", at(single, "hi"), 0.432, 0.0005)
check("single, retention, 0.75 A", at(single, "retention"), 0.235, 0.0005)
check("real-alpha dihedral median (rad)", np.median(d_real), 0.835, 0.0005)
check("constant-alpha dihedral median (rad)", np.median(d_const), 0.831, 0.0005)
if failures:
    sys.exit("Error: figure value disagrees with the paper: " + "; ".join(failures))

# ---- draw ----------------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(14, 11.6))
fig.subplots_adjust(left=0.08, right=0.97, top=0.86, bottom=0.07, hspace=0.32, wspace=0.22)
fig.suptitle("Paired αβ versus single-chain β — bin CDRB(10,8,13), "
             f"n = {len(ids):,}, single linkage", x=0.02, ha="left",
             fontsize=15, fontweight="bold")
ARMS = [(paired, RED, "Paired αβ  (all six CDR loops)"),
        (single, GREY, "Single-chain β  (β loops, β framework)")]
fig.legend(handles=[plt.Line2D([], [], color=c, lw=3, label=l) for _, c, l in ARMS],
           loc="upper left", bbox_to_anchor=(0.06, 0.94), ncol=2, frameon=False, fontsize=12)

ax = axs[0, 0]
for df, col, _ in ARMS:
    ax.plot(df.retention, df.hi, "-o", color=col, lw=2.5, ms=7)
    f = df[np.isclose(df.threshold, FOCAL)].iloc[0]
    ax.scatter(f.retention, f.hi, s=260, facecolors="none", edgecolors=col, lw=2)
    ax.annotate("0.75 Å", (f.retention, f.hi), xytext=(14, 6), textcoords="offset points",
                color=col, fontweight="bold", fontsize=11)
ax.set(xlabel="retention", ylabel="clusters with purity ≥ 0.9", ylim=(0, 1.06), xlim=(0, None))
ax.set_title("A · Purity vs retention")

for ax, col_name, ttl, ylab_, ylim in [
        (axs[0, 1], "hi", "B · Purity", "clusters with purity ≥ 0.9", (0, 1.06)),
        (axs[1, 0], "retention", "C · Retention", f"fraction of the {len(ids):,} clustered", (0, None))]:
    for df, col, _ in ARMS:
        ax.plot(df.threshold, df[col_name], "-o", color=col, lw=2.5, ms=7)
        v = at(df, col_name)
        ax.annotate(f"{v:.3f}", (FOCAL, v), xytext=(10, 8 if df is paired or col_name == "hi" else -16),
                    textcoords="offset points", color=col, fontweight="bold", fontsize=11)
    ax.axvline(FOCAL, color="#999999", ls="--", lw=1)
    ax.set(xlabel="threshold (Å)", ylabel=ylab_, ylim=ylim)
    ax.set_title(ttl)

ax = axs[1, 1]
bins = np.linspace(0.15, 1.5, 55)
ax.hist(d_real, bins=bins, histtype="step", color=RED, lw=2.2)
ax.hist(d_const, bins=bins, histtype="step", color=GREY, lw=2.2)
ax.axvline(np.median(d_real), color=RED, ls="--", lw=1.5)
ax.axvline(np.median(d_const), color=GREY, ls="--", lw=1.5)
ax.text(0.04, 0.93, f"medians {np.median(d_real):.3f} vs {np.median(d_const):.3f} rad\n"
        f"over {len(d_real) / 1e6:.2f}M pairs", transform=ax.transAxes, va="top",
        fontsize=11, color="#444444", linespacing=1.6)
ax.set(xlabel="pairwise CDR3β dihedral distance (rad)", ylabel="pairs")
ax.set_title("D · β-only spread")

for ax in axs.flat:
    ax.grid(color="#eeeeee"); ax.set_axisbelow(True)

path = os.path.join(OUT, "figS9_paired_vs_singlechain.png")
fig.savefig(path, dpi=300)
print(f"  saved: {os.path.basename(path)}")
