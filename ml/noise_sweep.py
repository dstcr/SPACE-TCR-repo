#!/usr/bin/env python3
"""Apply a distilled model to the Benchmark A noise ladder.

    python3 ml/noise_sweep.py \
        --train-assignments cluster_assignments.tsv \
        --train-labels data/jouannet_labelled_input.tsv \
        --ladder-seqs  <dir with donor1_subsetN_sptce.tsv> \
        --ladder-labels labelled_cdr3b.csv \
        --out ml_noise

The model is trained on the curated CD8 cohort's structural clustering and
applied unchanged to the 10X donor-1 noise ladder. The two datasets share no
receptors, so no split is needed, and the setting matches the intended use:
sequence-only screening of a new repertoire.

Scoring follows Benchmark A for SPACE-TCR and the published methods (same
labels, purity definition and minimum cluster size), so the curve can be drawn
on the same axes as Figure 9.

Rung 6 holds more than 20,000 receptors, so neighbours within the threshold are
found by chunked matrix multiplication on the L2-normalised embedding and
components by scipy's connected_components, which takes seconds per rung.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics as M                                       # noqa: E402
from distill import Cohort, embed, load_cohort, strip_allele, train  # noqa: E402


def radius_components(Z, thr, chunk=2000):
    """Connected components of the graph joining points within `thr`.

    Z must be L2-normalised, so squared euclidean distance is 2 - 2*dot and a
    dot-product threshold is equivalent to a distance threshold.
    """
    n = len(Z)
    dot_min = 1.0 - (thr ** 2) / 2.0
    rows, cols = [], []
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        sim = Z[s:e] @ Z.T                       # (chunk, n)
        i, j = np.nonzero(sim >= dot_min)
        i = i + s
        keep = i < j                             # upper triangle only
        rows.append(i[keep]); cols.append(j[keep])
    if rows:
        rows = np.concatenate(rows); cols = np.concatenate(cols)
    else:
        rows = cols = np.array([], dtype=int)
    g = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    return lab


def score_rung(df, labels, cluster_ids, denom_labelled):
    """Benchmark A's scoring: purity over clusters of size >= 2, and the share
    of labelled receptors recovered. Unlabelled members count towards cluster
    size but not towards the dominant count, as .purity_multi() does in R."""
    a = pd.DataFrame({"receptor_id": df["receptor_id"].to_numpy(),
                      "cluster_id": np.asarray(cluster_ids, dtype=object),
                      "cdr3b": df["cdr3b"].to_numpy()})
    sizes = a.groupby("cluster_id")["receptor_id"].nunique()
    keep = set(sizes[sizes >= 2].index)
    a = a[a["cluster_id"].isin(keep)]
    if not len(a):
        return dict(n_clusters=0, n_clustered=0, purity=float("nan"),
                    mean_size=float("nan"), signal_recovery=0.0)
    lab = a.merge(labels, on="cdr3b", how="inner")
    dom = (lab.groupby(["cluster_id", "epitope"]).size()
              .groupby(level=0).max())
    tot = a.groupby("cluster_id")["receptor_id"].nunique()
    purity = float(dom.reindex(tot.index).fillna(0).sum() / tot.sum())
    return dict(n_clusters=int(tot.size),
                n_clustered=int(tot.sum()),
                purity=purity,
                mean_size=float(tot.mean()),
                signal_recovery=float(lab["cdr3b"].nunique() / max(denom_labelled, 1)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train-assignments", required=True)
    ap.add_argument("--train-labels", required=True)
    ap.add_argument("--ladder-seqs", required=True,
                    help="directory holding donor1_subsetN_sptce.tsv")
    ap.add_argument("--ladder-labels", required=True,
                    help="TSV/CSV with cdr3b and epitope columns")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--thresholds", default="0.012,0.020,0.028,0.043")
    ap.add_argument("--out", default="ml_noise")
    a = ap.parse_args(argv)

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    # ---- train on the curated cohort ---------------------------------------
    coh = load_cohort(a.train_assignments, a.train_labels)
    print(f"  training cohort: {len(coh.df)} receptors, "
          f"{coh.df['cluster_id'].nunique()} clusters")
    model, hist = train(coh, coh.df, epochs=a.epochs, seed=a.seed, quiet=False)

    # ---- the ladder ---------------------------------------------------------
    lab = pd.read_csv(a.ladder_labels, sep=None, engine="python")
    lab.columns = [c.lower() for c in lab.columns]
    lab = lab.rename(columns={c: "cdr3b" for c in lab.columns if "cdr3" in c})
    lab = lab[["cdr3b", "epitope"]].dropna().drop_duplicates()
    n_lab_total = lab["cdr3b"].nunique()
    print(f"  ladder labels: {n_lab_total} labelled CDR3b\n")

    thrs = [float(t) for t in a.thresholds.split(",")]
    rows = []
    for i in range(1, 7):
        f = Path(a.ladder_seqs) / f"donor1_subset{i}_sptce.tsv"
        if not f.exists():
            print(f"  [skip] rung {i}: {f.name} absent")
            continue
        s = pd.read_csv(f, sep="\t")
        df = pd.DataFrame({
            "receptor_id": s["id"].astype(str),
            "cdr3b": s["cdr3_aa_beta"].astype(str),
            "cdr3a": s["cdr3_aa_alpha"].astype(str),
            "epitope": "",
            "v_gene": s["v_call_beta"].map(strip_allele),
            "j_gene": s["j_call_beta"].map(strip_allele)})
        df = df[df["cdr3b"].str.len().between(6, 22)].reset_index(drop=True)
        n_here = df["cdr3b"].isin(set(lab["cdr3b"])).sum()
        noise = 100.0 * (1 - n_here / len(df))

        Z = embed(model, df, coh.gene_vocab)
        Z = Z / np.clip(np.linalg.norm(Z, axis=1, keepdims=True), 1e-9, None)
        Z = np.ascontiguousarray(Z, dtype=np.float32)

        for thr in thrs:
            cl = radius_components(Z, thr)
            r = score_rung(df, lab, cl, n_lab_total)
            r.update(subset=i, n_input=len(df), n_labelled=int(n_here),
                     pct_noise=round(noise, 1), threshold=thr,
                     method="distilled model")
            rows.append(r)
            print(f"  rung {i} ({noise:5.1f}% noise, n={len(df):6d})  thr {thr:.3f}  "
                  f"clusters {r['n_clusters']:6d}  purity {r['purity']:.3f}  "
                  f"mean size {r['mean_size']:.2f}  recovery {100*r['signal_recovery']:.1f}%")
        print()

    tab = pd.DataFrame(rows)
    tab.to_csv(out / "noise_sweep_model.csv", index=False)
    (out / "run.json").write_text(json.dumps(
        {"train_assignments": a.train_assignments, "epochs": a.epochs,
         "seed": a.seed, "thresholds": thrs,
         "n_train_receptors": int(len(coh.df)),
         "n_train_clusters": int(coh.df["cluster_id"].nunique())}, indent=2))
    print(f"  wrote {out/'noise_sweep_model.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
