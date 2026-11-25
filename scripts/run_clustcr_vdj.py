#!/usr/bin/env python3
"""
run_clustcr_vdj.py  (v4 — correct ClusTCR API from source inspection)

ClusTCR fit() signature:
    fit(data, include_vgene=False, cdr3_col=None, v_gene_col=None, alpha=None)

Without vgene:
    fit(pd.Series of CDR3 strings)
    clusters_df columns: junction_aa, cluster

With vgene:
    fit(df, include_vgene=True, cdr3_col='my_cdr3_col', v_gene_col='my_v_col')
    clusters_df columns: junction_aa, v_call, cluster

With paired alpha+beta (no vgene):
    fit(beta_series, alpha=alpha_series)
    -- concatenates alpha to beta internally; not a true paired clustering

Additional Clustering() init params:
    chain='A' or 'B'  -- used for downstream analysis only, not clustering
    method='mcl'|'faiss'|'two-step'
    n_cpus            -- only used for two-step method
"""

from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import pandas as pd
from clustcr import Clustering

VFAM_RE = re.compile(r"^(TR[AB]V\d+)", re.IGNORECASE)

# ── Helpers ───────────────────────────────────────────────────────────────────

def v_family(v):
    if not v or pd.isna(v): return None
    v = str(v).split("*")[0].strip()
    m = VFAM_RE.match(v)
    return m.group(1).upper() if m else v.upper()

def read_input(path):
    df = pd.read_csv(path, dtype=str)
    df = df.loc[:, ~df.columns.str.match(r"^Unnamed")]
    df.columns = df.columns.str.strip()
    rename = {}
    for col in df.columns:
        cl = col.lower()
        if cl == "id":                       rename[col] = "id"
        elif "cdr3" in cl and "alpha" in cl: rename[col] = "cdr3_aa_alpha"
        elif "cdr3" in cl and "beta"  in cl: rename[col] = "cdr3_aa_beta"
        elif cl.startswith("v_") and "alpha" in cl: rename[col] = "v_call_alpha"
        elif cl.startswith("j_") and "alpha" in cl: rename[col] = "j_call_alpha"
        elif cl.startswith("v_") and "beta"  in cl: rename[col] = "v_call_beta"
        elif cl.startswith("j_") and "beta"  in cl: rename[col] = "j_call_beta"
    df = df.rename(columns=rename)
    for col in ["cdr3_aa_beta", "cdr3_aa_alpha"]:
        if col in df.columns:
            df[col] = df[col].str.strip()
    if "v_call_beta"  in df.columns: df["v_family_beta"]  = df["v_call_beta"].apply(v_family)
    if "v_call_alpha" in df.columns: df["v_family_alpha"] = df["v_call_alpha"].apply(v_family)
    return df

def choose_method(n_unique, requested):
    if requested != "auto": return requested
    return "mcl" if n_unique < 50_000 else "two-step"

# ── Mode definitions ──────────────────────────────────────────────────────────
# Each mode is a dict describing how to call ClusTCR fit()
#
# Modes:
#   beta_novgene   — CDR3β Series only
#   beta_vgene     — CDR3β + TRBV family stratification
#   alpha_novgene  — CDR3α Series only
#   alpha_vgene    — CDR3α + TRAV family stratification
#   paired         — CDR3β concatenated with CDR3α (fit(beta, alpha=alpha))
#                    This is a positional concatenation, not true paired clustering

def build_modes(method_override):
    modes = []
    for method in ("mcl", "faiss", "two-step"):
        if method_override not in ("auto", method):
            continue
        slug = method.replace("-", "_")
        for chain_mode in ("beta_novgene", "beta_vgene",
                           "alpha_novgene", "alpha_vgene",
                           "paired"):
            modes.append(dict(
                label      = f"{slug}_{chain_mode}",
                method     = method,
                chain_mode = chain_mode,
            ))
    return modes

# ── ClusTCR calls ─────────────────────────────────────────────────────────────

def fit_novgene(cl, cdr3_series):
    """fit(pd.Series)"""
    return cl.fit(cdr3_series)

def fit_vgene(cl, df, cdr3_col, v_gene_col):
    """
    fit(df, include_vgene=True, cdr3_col=..., v_gene_col=...)
    ClusTCR._vgene_clustering internally does:
        data["v_family"] = data[v_gene_col].apply(lambda x: x.split("-")[0].split("*")[0])
        remap = {cdr3_col: "junction_aa", v_gene_col: "v_call"}
    So v_gene_col values must already be family-level (TRBV6 not TRBV6-1*01)
    because it will try to split on "-" again — which is fine for family strings.
    """
    return cl.fit(df, include_vgene=True, cdr3_col=cdr3_col, v_gene_col=v_gene_col)

def fit_paired(cl, beta_series, alpha_series):
    """fit(beta_series, alpha=alpha_series) — concatenates α+β internally"""
    return cl.fit(beta_series, alpha=alpha_series)

# ── Parse clusters_df ─────────────────────────────────────────────────────────

def parse_clusters_df(clusters_df):
    """
    Normalise clusters_df to always have:
        junction_aa   — CDR3 sequence
        clustcr_cluster — cluster ID
        v_call (optional) — V gene (only in vgene mode)
    """
    df = clusters_df.copy()
    # column names from source: junction_aa, cluster, v_call (optional)
    if "cluster" in df.columns:
        df = df.rename(columns={"cluster": "clustcr_cluster"})
    return df

# ── Merge back to original rows ───────────────────────────────────────────────

def merge_back(df, clusters_df, cdr3_col, vgene_col, chain_mode, min_size):
    """
    Join cluster assignments back to original rows.
    Handles vgene / novgene / paired modes.
    """
    # Enforce min cluster size on junction_aa level
    sizes = clusters_df.groupby("clustcr_cluster")["junction_aa"].nunique()
    keep  = set(sizes[sizes >= min_size].index)
    c = clusters_df.copy()
    c["clustcr_cluster"] = c["clustcr_cluster"].where(
        c["clustcr_cluster"].isin(keep), other=-1)

    if "vgene" in chain_mode and vgene_col and "v_call" in c.columns:
        # Join on CDR3 + V gene
        merged = df.merge(
            c[["junction_aa", "v_call", "clustcr_cluster"]].drop_duplicates(),
            left_on=[cdr3_col, vgene_col],
            right_on=["junction_aa", "v_call"],
            how="left"
        )
    elif chain_mode == "paired":
        # Paired mode: ClusTCR concatenated β+α; the junction_aa in clusters_df
        # will be the concatenated sequence. Join on concatenated key.
        df = df.copy()
        df["_paired_key"] = df["cdr3_aa_beta"].fillna("") + df["cdr3_aa_alpha"].fillna("")
        c_paired = c[["junction_aa", "clustcr_cluster"]].drop_duplicates()
        merged = df.merge(
            c_paired,
            left_on="_paired_key",
            right_on="junction_aa",
            how="left"
        ).drop(columns=["_paired_key"])
    else:
        # Join on CDR3 only
        merged = df.merge(
            c[["junction_aa", "clustcr_cluster"]].drop_duplicates(),
            left_on=cdr3_col,
            right_on="junction_aa",
            how="left"
        )

    merged["clustcr_cluster"] = merged["clustcr_cluster"].fillna(-1).astype(int)
    return merged

# ── Save ──────────────────────────────────────────────────────────────────────

def save_mode(merged, clusters_df, output, cdr3_col, outdir, label):
    d = outdir / label
    d.mkdir(parents=True, exist_ok=True)
    merged.to_csv(d / "assignments.tsv", sep="\t", index=False)
    clusters_df.to_csv(d / "clusters.tsv", sep="\t", index=False)
    try:
        output.summary().to_csv(d / "summary.tsv", sep="\t", index=False)
    except Exception:
        pass
    cl = merged[merged["clustcr_cluster"] != -1]
    if not cl.empty:
        agg_col = cdr3_col if cdr3_col in cl.columns else "id"
        cl.groupby("clustcr_cluster").agg(
            n_rows      = ("id",     "size"),
            n_unique_id = ("id",     "nunique"),
        ).reset_index().sort_values("n_unique_id", ascending=False).to_csv(
            d / "stats.tsv", sep="\t", index=False)
    n_c = int((merged["clustcr_cluster"] != -1).sum())
    n_t = len(merged)
    print(f"  [{label}] clustered {n_c:,}/{n_t:,} ({n_c/max(n_t,1):.3f}) → {d}",
          file=sys.stderr)
    return n_c

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Run ClusTCR on all modes for a paired VDJ CSV."
    )
    ap.add_argument("--input",            required=True, type=Path)
    ap.add_argument("--outdir",           required=True, type=Path)
    ap.add_argument("--method",           default="auto",
                    choices=["auto", "mcl", "faiss", "two-step"])
    ap.add_argument("--n-cpus",           default="all")
    ap.add_argument("--min-cluster-size", type=int, default=2)
    args = ap.parse_args()

    n_cpus = "all" if str(args.n_cpus).lower() == "all" else int(args.n_cpus)

    print(f"\n[INFO] Reading {args.input}", file=sys.stderr)
    df = read_input(args.input)
    print(f"[INFO] Rows: {len(df):,}  Cols: {list(df.columns)}", file=sys.stderr)

    has_alpha  = "cdr3_aa_alpha"  in df.columns
    has_vbeta  = "v_family_beta"  in df.columns
    has_valpha = "v_family_alpha" in df.columns

    modes = build_modes(args.method)
    print(f"[INFO] {len(modes)} modes planned\n", file=sys.stderr)

    summary = []

    for m in modes:
        label      = m["label"]
        method     = m["method"]
        chain_mode = m["chain_mode"]

        # ── Skip if required columns are absent ───────────────────────────────
        if "alpha" in chain_mode and not has_alpha:
            print(f"  [{label}] SKIP — no cdr3_aa_alpha", file=sys.stderr)
            continue
        if chain_mode == "paired" and not has_alpha:
            print(f"  [{label}] SKIP — no cdr3_aa_alpha for paired mode", file=sys.stderr)
            continue
        if "vgene" in chain_mode and "beta" in chain_mode and not has_vbeta:
            print(f"  [{label}] SKIP — no v_call_beta", file=sys.stderr)
            continue
        if "vgene" in chain_mode and "alpha" in chain_mode and not has_valpha:
            print(f"  [{label}] SKIP — no v_call_alpha", file=sys.stderr)
            continue

        # ── Determine CDR3 and V gene columns ─────────────────────────────────
        if "alpha" in chain_mode:
            cdr3_col  = "cdr3_aa_alpha"
            vgene_col = "v_family_alpha"
            chain_arg = "A"
        elif chain_mode == "paired":
            cdr3_col  = "cdr3_aa_beta"
            vgene_col = None
            chain_arg = "B"
        else:
            cdr3_col  = "cdr3_aa_beta"
            vgene_col = "v_family_beta"
            chain_arg = "B"

        # ── Build input for ClusTCR ───────────────────────────────────────────
        if "vgene" in chain_mode:
            # DataFrame with cdr3_col and vgene_col columns
            sub = df[[cdr3_col, vgene_col]].dropna()
            sub = sub[(sub[cdr3_col] != "") & (sub[vgene_col] != "")]
            unique_keys = sub.drop_duplicates()
            n_unique = len(unique_keys)
        elif chain_mode == "paired":
            sub = df[["cdr3_aa_beta", "cdr3_aa_alpha"]].dropna()
            sub = sub[(sub["cdr3_aa_beta"] != "") & (sub["cdr3_aa_alpha"] != "")]
            unique_keys = sub.drop_duplicates()
            n_unique = len(unique_keys)
        else:
            sub = df[[cdr3_col]].dropna()
            sub = sub[sub[cdr3_col] != ""]
            unique_keys = sub.drop_duplicates()
            n_unique = len(unique_keys)

        actual_method = choose_method(n_unique, method)
        print(f"\n  [{label}] n_unique={n_unique:,} method={actual_method} chain_mode={chain_mode}",
              file=sys.stderr)

        # ── Instantiate Clustering ────────────────────────────────────────────
        cl = Clustering(chain=chain_arg, method=actual_method, n_cpus=n_cpus)

        # ── Call fit() ────────────────────────────────────────────────────────
        try:
            if "vgene" in chain_mode:
                output = fit_vgene(cl, unique_keys, cdr3_col, vgene_col)
            elif chain_mode == "paired":
                beta_series  = unique_keys["cdr3_aa_beta"].reset_index(drop=True)
                alpha_series = unique_keys["cdr3_aa_alpha"].reset_index(drop=True)
                output = fit_paired(cl, beta_series, alpha_series)
            else:
                cdr3_series = unique_keys[cdr3_col].reset_index(drop=True)
                output = fit_novgene(cl, cdr3_series)
        except Exception as e:
            print(f"  [{label}] ERROR: {e}", file=sys.stderr)
            summary.append(dict(mode=label, status=f"FAILED:{e}",
                                n_clustered=0, n_total=len(df),
                                frac_clustered=0.0, n_clusters=0))
            continue

        # ── Normalise clusters_df and merge back ──────────────────────────────
        clusters_df = parse_clusters_df(output.clusters_df)
        merged = merge_back(df, clusters_df, cdr3_col, vgene_col,
                            chain_mode, args.min_cluster_size)
        n_c  = save_mode(merged, clusters_df, output, cdr3_col, args.outdir, label)
        n_cl = int(merged.loc[merged["clustcr_cluster"] != -1,
                               "clustcr_cluster"].nunique())
        summary.append(dict(
            mode           = label,
            status         = "OK",
            method         = actual_method,
            chain_mode     = chain_mode,
            n_total        = len(merged),
            n_clustered    = n_c,
            frac_clustered = round(n_c / max(len(merged), 1), 4),
            n_clusters     = n_cl,
        ))

    # ── Final summary ─────────────────────────────────────────────────────────
    print("\n\n=== Done ===\n", file=sys.stderr)
    sdf = pd.DataFrame(summary)
    print(sdf.to_string(index=False), file=sys.stderr)
    args.outdir.mkdir(parents=True, exist_ok=True)
    sdf.to_csv(args.outdir / "run_summary.tsv", sep="\t", index=False)
    print(f"\n[OK] {args.outdir}/run_summary.tsv", file=sys.stderr)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
