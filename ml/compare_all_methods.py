#!/usr/bin/env python3
"""Score the distilled model against SPACE-TCR and all nine published methods
on the same receptors, with the same equations.

    python3 ml/compare_all_methods.py \
        --assignments cluster_assignments.tsv \
        --labels data/jouannet_labelled_input.tsv \
        --published ml/data/published_clusters.tsv.gz \
        --epochs 40 --seed 0 --out ml_compare

`distill.py` compares the model only with Levenshtein distance. This places it
in the full method comparison.

The comparison is restricted to the held-out test receptors, whose epitopes the
model has not seen, and every method is scored on that same receptor set with
the same denominator. This also means:

  - Units: five published methods were run on a chain-pooled list, so their
    output mixes CDR3alpha and CDR3beta. Intersecting on the test receptors'
    CDR3beta drops the alpha rows, so every method is scored per test receptor.
  - Thresholds: the published methods give one fixed clustering, while the
    model and Levenshtein have a free threshold, so their sweep is reported with
    the operating point marked.

See ml/README.md for the caveats. On this cohort the model's advantage over
edit distance is +0.041 +/- 0.041 AUROC across five seeds (p = 0.092), so
single-seed cluster metrics vary; pass --seeds to average several.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import metrics as M                                      # noqa: E402
from distill import (cluster_by_threshold, embed, epitope_split,  # noqa: E402
                     ld_distance_fn, load_cohort, score_clustering,
                     strip_allele, train)

# Five of the nine were run on a chain-pooled single-chain list. Recorded so
# the output can flag them, even though the
# restriction below means they are scored in the same unit as everyone else.
CHAIN_POOLED = {"GIANA", "iSMART", "TCRMatch",
                "Hamming_Distance", "Levenshtein_Distance"}


def score(df, cluster_ids, denom, label, extra=None):
    """One row of the table.

    Delegates to distill.score_clustering so that every number here is produced
    by the same code distill.py reports, which in turn is ml/metrics.py, which
    is checked against the R by ml/tests/test_metrics_match_r.py. A second
    implementation would be a second thing to keep in step.

    `cluster_ids` is a Series indexed by receptor_id, as cluster_by_threshold
    returns; NaN means unclustered.
    """
    if not isinstance(cluster_ids, pd.Series):
        cluster_ids = pd.Series(list(cluster_ids), index=df["receptor_id"])
    out = score_clustering(df, cluster_ids, denom, label)
    row = {
        "method": label,
        "n_clusters": out["n_clusters"],
        "clustered": out["n_units_clustered"],
        "retention": out["retention"],
        "purity": out["purity"],
        "pct_clusters_hi": out["pct_clusters_hi"],
        "pct_units_hi": out["pct_units_hi"],
        "polyspecificity": out["polyspecificity"],
        # the paper's distinctive claim: the share of
        # high-purity clusters too CDR3b-divergent for an edit distance to have
        # formed them.
        "pct_hi_divergent": out.get("pct_hi_clusters_sequence_diverse",
                                    float("nan")),
        "median_ld_hi": out.get("median_ld_of_hi_clusters", float("nan")),
    }
    row.update(extra or {})
    return row


def sweep(df, dist_fn, denom, label, lo, hi, n):
    """Sweep a distance threshold; return every row and the operating point."""
    ids = df["receptor_id"].tolist()
    rows = []
    for t in np.linspace(lo, hi, n):
        cl = cluster_by_threshold(ids, dist_fn, float(t))
        rows.append(score(df, cl, denom, label, {"threshold": round(float(t), 4)}))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--published", required=True,
                    help="TSV of method, cluster_id, cdr3b (see ml/data/)")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=1,
                    help="average the model over this many seeds from --seed")
    ap.add_argument("--test-frac", type=float, default=0.3)
    ap.add_argument("--universe", choices=("labelled", "cohort"), default="labelled",
                    help="which receptors to score on. 'labelled' (default) is "
                         "EVERY receptor in --labels carrying a held-out epitope. "
                         "'cohort' is only those SPACE-TCR already clustered, "
                         "which is circular -- see the note in the source.")
    ap.add_argument("--out", default="ml_compare")
    # The model's useful threshold range depends on the scale of its embedding,
    # which is not known in advance. The default 0.05-1.6 was chosen for a
    # well-spread encoder; a tighter one puts every interesting operating point
    # below 0.05 and the sweep returns one usable row. Widen or narrow here.
    ap.add_argument("--model-lo", type=float, default=0.05)
    ap.add_argument("--model-hi", type=float, default=1.6)
    ap.add_argument("--model-steps", type=int, default=12)
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cohort = load_cohort(a.assignments, a.labels)
    df = cohort.df
    print(f"  cohort: {len(df)} receptors, {df['cluster_id'].nunique()} clusters, "
          f"{df['epitope'].nunique()} epitopes")

    tr, te, split_info = epitope_split(df, frac_test=a.test_frac, seed=a.seed)

    # ---- the evaluation universe -------------------------------------------
    # load_cohort() keeps only receptors SPACE-TCR already placed in a cluster
    # of size >= 2 (1,555 of the cohort's 4,779). Scoring on that set makes
    # SPACE-TCR's retention ~100% BY CONSTRUCTION and asks every other method
    # to cluster a set SPACE-TCR chose. So by default the universe is every
    # receptor in --labels carrying a held-out epitope, clustered or not:
    # receptors SPACE-TCR missed now count against it, as they should.
    #
    # Restriction also inflates purity, and unevenly. Cutting a cluster down to
    # the held-out receptors removes its other-epitope members, so an
    # aggressive clusterer is flattered most: on this split DeepTCR's purity
    # goes 61.1% -> 96.6% and TCRdist3's 42.4% -> 73.0%, while clusTCR's moves
    # 98.0% -> 99.2%. That is unavoidable when evaluating on held-out epitopes,
    # but it means NO purity here is comparable to a full-cohort number, and
    # the ranking is not either. Read this table against itself only.
    test_epis = set(te["epitope"])
    if a.universe == "labelled":
        lab_all = pd.read_csv(a.labels, sep=None, engine="python")
        lab_all = lab_all.rename(columns={
            "id": "receptor_id", "cdr3_aa_beta": "cdr3b", "cdr3_aa_alpha": "cdr3a"})
        ev = lab_all[lab_all["epitope"].isin(test_epis)].copy()
        ev["receptor_id"] = ev["receptor_id"].astype(str)
        # drop the leak distill.py drops: test receptors sharing a CDR3b with train
        shared = set(tr["cdr3b"]) & set(ev["cdr3b"])
        n_leak = int(ev["cdr3b"].isin(shared).sum())
        ev = ev[~ev["cdr3b"].isin(shared)].reset_index(drop=True)
        # SPACE-TCR's own clustering, with NaN for what it never clustered
        ev = ev.merge(df[["receptor_id", "cluster_id"]], on="receptor_id", how="left")
        # the model needs the gene columns load_cohort() builds
        for src, dst in (("v_call_beta", "v_gene"), ("j_call_beta", "j_gene")):
            ev[dst] = ev[src].map(strip_allele) if src in ev.columns else "UNK"
        te = ev
        print(f"  universe: all {len(te)} receptors carrying a held-out epitope "
              f"({n_leak} dropped for sharing a CDR3b with train)")
        n_unclustered = int(te["cluster_id"].isna().sum())
        print(f"  of those, SPACE-TCR clustered {len(te) - n_unclustered} "
              f"and missed {n_unclustered}")
    else:
        print(f"  universe: the {len(te)} SPACE-TCR-clustered test receptors "
              f"(CIRCULAR -- SPACE-TCR retention is ~100% by construction)")
    denom = len(te)
    print(f"  train {len(tr)} / test {denom} receptors, "
          f"{te['epitope'].nunique()} held-out epitopes")
    print(f"  every method below is cut to those {denom} receptors "
          f"and divided by {denom}.\n")

    rows, sweeps = [], {}

    # ---- the structural clustering, restricted to the test receptors --------
    # This is the reference: what SPACE-TCR itself achieves on exactly these
    # receptors. Every other row is trying to reproduce it from sequence.
    rows.append(score(te, pd.Series(te["cluster_id"].to_numpy(),
                                    index=te["receptor_id"]),
                      denom, "SPACE-TCR (structural)", {"threshold": np.nan}))

    # ---- the published methods, restricted to the test receptors -----------
    pub = pd.read_csv(a.published, sep="\t")
    need = {"method", "cluster_id", "cdr3b"}
    if not need.issubset(pub.columns):
        raise SystemExit(f"--published needs {sorted(need)}; got {list(pub.columns)}")

    # Join on CDR3beta: a test receptor is "clustered" by a method if its beta
    # appears in one of that method's clusters. Alpha rows from the five
    # chain-pooled methods do not match, so all methods are scored per receptor.
    beta_to_rid = (te.drop_duplicates("cdr3b")
                     .set_index("cdr3b")["receptor_id"].to_dict())
    for name, g in pub.groupby("method", sort=False):
        g = g[g["cdr3b"].isin(beta_to_rid)].copy()
        if g.empty:
            print(f"  [skip] {name}: no test receptors")
            continue
        g["receptor_id"] = g["cdr3b"].map(beta_to_rid)
        # One cluster per receptor: a method that assigns a beta to several
        # clusters keeps the first, so no receptor is counted twice.
        g = g.drop_duplicates("receptor_id")
        cl = pd.Series(
            te["receptor_id"].map(g.set_index("receptor_id")["cluster_id"]).to_numpy(),
            index=te["receptor_id"])
        rows.append(score(te, cl, denom, name,
                          {"threshold": np.nan,
                           "chain_pooled_run": name in CHAIN_POOLED}))

    # ---- Levenshtein on CDR3beta -------------------------------------------
    ld_fn = ld_distance_fn(te["cdr3b"].tolist())
    sweeps["Levenshtein CDR3b"] = sweep(te, ld_fn, denom,
                                        "Levenshtein CDR3b", 1, 8, 8)

    # ---- the distilled model ------------------------------------------------
    model_sweeps = []
    for k in range(a.seeds):
        sd = a.seed + k
        print(f"\n  training model, seed {sd} ...")
        tr_k, _, _ = epitope_split(df, frac_test=a.test_frac, seed=a.seed)
        model, _ = train(cohort, tr_k, epochs=a.epochs, seed=sd, quiet=False)
        Z = embed(model, te, cohort.gene_vocab)
        Z = Z / np.clip(np.linalg.norm(Z, axis=1, keepdims=True), 1e-9, None)

        def mdist(i, j, Z=Z):
            return float(np.linalg.norm(Z[i] - Z[j]))

        model_sweeps.append(sweep(te, mdist, denom, "distilled model",
                                  a.model_lo, a.model_hi, a.model_steps))
    # average the sweep across seeds, threshold by threshold
    if model_sweeps:
        base = model_sweeps[0]
        avg = []
        for i in range(len(base)):
            r = dict(base[i])
            for k in ("retention", "purity", "pct_clusters_hi", "pct_units_hi",
                      "polyspecificity", "pct_hi_divergent"):
                vals = [s[i][k] for s in model_sweeps]
                vals = [v for v in vals if v == v]
                r[k] = float(np.mean(vals)) if vals else float("nan")
            r["n_clusters"] = int(np.mean([s[i]["n_clusters"] for s in model_sweeps]))
            r["n_seeds"] = len(model_sweeps)
            avg.append(r)
        sweeps["distilled model"] = avg

    # ---- operating point: loosest threshold still at the structural purity --
    truth_purity = rows[0]["purity"]
    for label, sw in sweeps.items():
        ok = [r for r in sw if r["purity"] >= truth_purity]
        pick = max(ok, key=lambda r: r["retention"]) if ok else \
            max(sw, key=lambda r: r["purity"])
        pick = dict(pick)
        pick["reached_truth_purity"] = bool(ok)
        rows.append(pick)

    tab = pd.DataFrame(rows)
    tab["n_clusters_hi"] = (tab["n_clusters"] * tab["pct_clusters_hi"]).round().astype(int)
    tab.to_csv(out / "comparison.csv", index=False)
    for label, sw in sweeps.items():
        pd.DataFrame(sw).to_csv(
            out / f"sweep_{label.split()[0].lower()}.csv", index=False)

    # ---- print ---------------------------------------------------------------
    print(f"\n  all methods on the same {denom} held-out receptors, "
          f"same denominator, same equations\n")
    # hiN is the BASE the div column is computed over. Without it a method
    # with one high-purity cluster can print a 100% divergent share and look
    # like a winner; the distilled model does exactly that.
    hdr = (f"  {'method':<26}{'clu':>5}{'kept':>6}{'reten':>7}{'purity':>8}"
           f"{'hiClu':>7}{'hiN':>5}{'poly':>7}{'div':>8}  thr")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    tab2 = tab.sort_values("purity", ascending=False)
    for _, r in tab2.iterrows():
        thr = "" if r["threshold"] != r["threshold"] else f"{r['threshold']:.2f}"
        star = "*" if r.get("chain_pooled_run") is True else " "
        div = "" if r["pct_hi_divergent"] != r["pct_hi_divergent"] \
            else f"{100*r['pct_hi_divergent']:5.1f}%"
        hi_n = int(round(r["n_clusters"] * r["pct_clusters_hi"]))
        # flag a divergent share resting on too few clusters to mean anything
        if div and hi_n < 5:
            div = div.strip() + "!"
        print(f"  {r['method']:<25}{star}{r['n_clusters']:>5}{r['clustered']:>6}"
              f"{100*r['retention']:>6.1f}%{100*r['purity']:>7.1f}%"
              f"{100*r['pct_clusters_hi']:>6.1f}%{hi_n:>5}"
              f"{100*r['polyspecificity']:>6.1f}%{div:>8}  {thr}")
    print("\n  * run on a chain-pooled list; the")
    print("    restriction to test receptors drops its alpha rows, so it is")
    print("    scored here in the same unit as everything else.")
    print("  div = share of high-purity clusters whose median within-cluster")
    print("        CDR3b Levenshtein distance is >= 4, i.e. groups no edit")
    print("        distance could have made. This is the paper's claim.")
    print("  hiN = how many high-purity clusters that share is computed over.")
    print("  !   = fewer than 5, so the div figure means nothing. Read hiN")
    print("        before reading div: a method with one high-purity cluster")
    print("        can print 100% and look like a winner.")
    print(f"\n  wrote {out/'comparison.csv'}\n")
    (out / "run.json").write_text(json.dumps(
        {"denom": denom, "n_train": len(tr), "seed": a.seed,
         "seeds": a.seeds, "epochs": a.epochs,
         "universe": a.universe, "split": split_info}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
