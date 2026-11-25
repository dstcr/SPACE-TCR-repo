#!/usr/bin/env python
"""
SPTCE Stage E — clusTCR clustering (sequence-based, alternative to Stage C / SPACE2-TCR).

Clusters TCRs by CDR3 amino-acid sequence using clusTCR (Valkiers et al. 2021),
following the input convention used in the TCR_Unsupervised_Benchmark study
(Jouannet et al.): paired CDR3a/CDR3b where available, beta-only otherwise.

Must be run with the clustcr_env interpreter and PYTHONNOUSERSITE=1 set
(handled by sbatch/run_E_slurm.sh).

Outputs (mirrors the benchmark's clusTCR output triad):
  <outprefix>_output.csv          one row per input seq with assigned cluster
  <outprefix>_output_edgelist.txt cluster network edges (if available)
  <outprefix>_output_summary.csv  one row per cluster with size + motif
"""
import argparse
import os
import sys
import warnings

# clusTCR's deps pin older numpy/scipy; the installed versions work fine.
# Silence the cosmetic version-mismatch warning seen in the env.
warnings.filterwarnings("ignore", message=r".*NumPy version.*")

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Run clusTCR on a TCR table.")
    p.add_argument("--input", required=True,
                   help="Input TSV/CSV with CDR3 (and optional V gene) columns.")
    p.add_argument("--outprefix", required=True,
                   help="Output path prefix; files get _output.csv etc. appended.")
    p.add_argument("--method", default="mcl", choices=["mcl", "two-step"],
                   help="clusTCR clustering method (default: mcl).")
    p.add_argument("--n-cpus", default="all",
                   help="Number of CPUs, or 'all' (default: all).")
    p.add_argument("--chain", default="beta", choices=["beta", "paired"],
                   help="Cluster on beta CDR3 only, or paired alpha+beta.")
    # Column names — defaults are guesses; override to match your TSV header.
    p.add_argument("--cdr3b-col", default="cdr3b",
                   help="Column holding the CDR3 beta amino-acid sequence.")
    p.add_argument("--vb-col", default="v_b_gene",
                   help="Column holding the TRBV gene (used by MCL for motif/V info).")
    p.add_argument("--cdr3a-col", default="cdr3a",
                   help="Column holding the CDR3 alpha sequence (paired mode).")
    p.add_argument("--va-col", default="v_a_gene",
                   help="Column holding the TRAV gene (paired mode).")
    p.add_argument("--id-col", default=None,
                   help="Optional ID column to carry through to the output.")
    p.add_argument("--sep", default=None,
                   help="Field separator; inferred from extension if omitted.")
    return p.parse_args()


def read_table(path, sep):
    if sep is None:
        sep = "\t" if path.lower().endswith((".tsv", ".tsv.gz", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = [c.strip() for c in df.columns]
    return df


def require_cols(df, cols, path):
    missing = [c for c in cols if c and c not in df.columns]
    if missing:
        sys.exit(
            f"ERROR: column(s) {missing} not found in {path}.\n"
            f"Available columns: {list(df.columns)}\n"
            f"Set the matching --*-col flags (or CLUSTCR_* vars in config.sh)."
        )


def main():
    args = parse_args()
    n_cpus = args.n_cpus
    if n_cpus != "all":
        n_cpus = int(n_cpus)

    df = read_table(args.input, args.sep)

    # Import here so the column-check error (common) doesn't wait on a slow import.
    from clustcr import Clustering

    if args.chain == "beta":
        require_cols(df, [args.cdr3b_col, args.vb_col, args.id_col], args.input)
        work = df.dropna(subset=[args.cdr3b_col]).copy()
        cdr3 = work[args.cdr3b_col].astype(str)
        clustering = Clustering(method=args.method, n_cpus=n_cpus)
        # clusTCR can use V-gene info for MCL; pass it when available.
        try:
            result = clustering.fit(cdr3, v_gene=work[args.vb_col].astype(str))
        except TypeError:
            # Older/newer signatures: fall back to CDR3-only fit.
            result = clustering.fit(cdr3)
    else:  # paired
        require_cols(
            df,
            [args.cdr3a_col, args.cdr3b_col, args.va_col, args.vb_col, args.id_col],
            args.input,
        )
        work = df.dropna(subset=[args.cdr3a_col, args.cdr3b_col]).copy()
        # clusTCR paired API expects a DataFrame with named CDR3a/CDR3b/V columns.
        paired = pd.DataFrame({
            "CDR3_alpha": work[args.cdr3a_col].astype(str),
            "CDR3_beta": work[args.cdr3b_col].astype(str),
            "V_alpha": work[args.va_col].astype(str),
            "V_beta": work[args.vb_col].astype(str),
        })
        clustering = Clustering(method=args.method, n_cpus=n_cpus)
        result = clustering.fit(paired, alpha=0.5)

    clusters_df = result.clusters_df.copy()

    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)) or ".", exist_ok=True)

    # Carry the ID column through by joining on the CDR3 sequence if requested.
    seq_col = "junction_aa" if "junction_aa" in clusters_df.columns else clusters_df.columns[0]
    if args.id_col:
        key = args.cdr3b_col if args.chain == "beta" else args.cdr3b_col
        id_map = work[[key, args.id_col]].drop_duplicates(subset=[key])
        clusters_df = clusters_df.merge(
            id_map, left_on=seq_col, right_on=key, how="left"
        )

    out_csv = f"{args.outprefix}_output.csv"
    clusters_df.to_csv(out_csv, index=False)
    print(f"[clustcr] wrote {out_csv}  ({len(clusters_df)} clustered sequences)")

    # Summary (one row per cluster). summary() returns motif/size per cluster.
    try:
        summary_df = result.summary()
        out_sum = f"{args.outprefix}_output_summary.csv"
        summary_df.to_csv(out_sum, index=False)
        n_clusters = clusters_df["cluster"].nunique() if "cluster" in clusters_df else len(summary_df)
        print(f"[clustcr] wrote {out_sum}  ({n_clusters} clusters)")
    except Exception as e:
        print(f"[clustcr] WARNING: could not write summary: {e}", file=sys.stderr)

    # Edge list (network). clusTCR exposes this differently across versions:
    # a DataFrame, a set of "seq1\tseq2" strings, or a set/list of pairs.
    try:
        edges = getattr(result, "edges", None)
        if edges is not None:
            out_edge = f"{args.outprefix}_output_edgelist.txt"
            if hasattr(edges, "to_csv"):                 # DataFrame
                edges.to_csv(out_edge, sep="\t", index=False, header=False)
            else:                                        # set or list
                with open(out_edge, "w") as fh:
                    for e in edges:
                        if isinstance(e, str):
                            line = e if "\t" in e else e.replace(" ", "\t")
                        else:                            # tuple/list pair
                            line = "\t".join(str(x) for x in e)
                        fh.write(line.rstrip("\n") + "\n")
            print(f"[clustcr] wrote {out_edge}")
    except Exception as e:
        print(f"[clustcr] WARNING: could not write edge list: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
