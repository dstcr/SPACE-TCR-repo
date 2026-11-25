#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

def is_aa(seq):
    return isinstance(seq, str) and seq != "" and all(ch in AA20 for ch in seq)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tra-tsv", required=True)
    ap.add_argument("--trb-tsv", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tra = pd.read_csv(args.tra_tsv, sep="\t", dtype=str)
    trb = pd.read_csv(args.trb_tsv, sep="\t", dtype=str)

    tra = tra[["TCR_name", "TRA_aa"]].copy()
    trb = trb[["TCR_name", "TRB_aa"]].copy()

    tra["TCR_name"] = tra["TCR_name"].astype(str).str.strip()
    trb["TCR_name"] = trb["TCR_name"].astype(str).str.strip()
    tra["TRA_aa"] = tra["TRA_aa"].astype(str).str.replace("*", "", regex=False).str.strip()
    trb["TRB_aa"] = trb["TRB_aa"].astype(str).str.replace("*", "", regex=False).str.strip()

    merged = tra.merge(trb, on="TCR_name", how="inner")

    wrote = 0
    skipped_bad_aa = 0
    skipped_bad_len = 0

    for _, r in merged.iterrows():
        rid = r["TCR_name"]
        a = r["TRA_aa"]
        b = r["TRB_aa"]

        if not is_aa(a) or not is_aa(b):
            skipped_bad_aa += 1
            continue

        # full-length stitched chains can be much longer
        if not (100 <= len(a) <= 500 and 100 <= len(b) <= 500):
            skipped_bad_len += 1
            continue

        out = out_dir / f"{rid}.fasta"
        if out.exists():
            continue

        with open(out, "w") as fh:
            fh.write(">B\n")
            fh.write(f"{b}\n")
            fh.write(">A\n")
            fh.write(f"{a}\n")

        wrote += 1

    print(f"Merged rows: {len(merged)}")
    print(f"Wrote FASTAs: {wrote}")
    print(f"Skipped bad AA rows: {skipped_bad_aa}")
    print(f"Skipped bad length rows: {skipped_bad_len}")
    print(f"Out dir: {out_dir}")

if __name__ == "__main__":
    main()
