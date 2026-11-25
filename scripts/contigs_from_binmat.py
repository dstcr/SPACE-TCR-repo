"""Build a contig-annotation-shaped table from the binarized matrix.

The real vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv was not
available. Its only role in jouannet_10x_labels() is to supply barcode ->
CDR3beta, and the binarized matrix already carries that per cell in
cell_clono_cdr3_aa as "TRA:...;TRA:...;TRB:...". This emits the same columns
the pipeline filters on, so the code path runs unmodified.

The substitution is NOT free: the contig route keeps productive,
high-confidence TRB contigs per contig, while this is the cell's clonotype
consensus. The pipeline's own alignment report against subset_1 is the check.
"""
import csv, sys

import argparse
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("binmat", help="vdj_v1_hs_aggregated_donor1_binarized_matrix.csv")
ap.add_argument("out",    help="where to write the contig-shaped CSV")
ap.add_argument("--donor", default="donor1")
a = ap.parse_args()
src, out = a.binmat, a.out

n_cells = n_rows = n_multi = 0
with open(src, newline="") as f, open(out, "w", newline="") as g:
    r = csv.DictReader(f)
    w = csv.writer(g)
    w.writerow(["barcode", "cdr3", "chain", "is_cell",
                "high_confidence", "productive"])
    for row in r:
        if row.get("donor") != a.donor:
            continue
        n_cells += 1
        clono = row.get("cell_clono_cdr3_aa") or ""
        betas = [p.split(":", 1)[1] for p in clono.split(";")
                 if p.startswith("TRB:") and ":" in p]
        if len(betas) > 1:
            n_multi += 1
        for b in betas:
            if 6 <= len(b) <= 22:
                w.writerow([row["barcode"], b, "TRB", "True", "True", "True"])
                n_rows += 1

print(f"{a.donor} cells      : {n_cells}")
print(f"TRB contig rows   : {n_rows}")
print(f"cells with >1 TRB : {n_multi}")
print(f"wrote {out}")
