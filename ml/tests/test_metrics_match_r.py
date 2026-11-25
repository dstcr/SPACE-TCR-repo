"""Check ml/metrics.py agrees with sptce_metrics.R to floating-point.

A learned model scored by a different equation cannot be compared with
SPACE-TCR. This generates random assignment tables, scores each in both
languages, and fails on any disagreement.

    python3 ml/tests/test_metrics_match_r.py        # from the repository root
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import metrics as M  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

R_DRIVER = r'''
args <- commandArgs(trailingOnly = TRUE)
csv <- args[1]; denom <- as.numeric(args[2]); out <- args[3]
setwd(Sys.getenv("SPTCE_REPO"))
suppressMessages(source("sptce_metrics.R"))
asg <- readr::read_csv(csv, show_col_types = FALSE)
asg$epitope[asg$epitope == "NA"] <- NA
m  <- jm_metrics(asg, denom = denom, unit = "receptor_id")
ct <- cluster_table(asg, unit = "receptor_id")
epis <- sort(unique(asg$epitope[!is.na(asg$epitope)]))
sens <- sapply(epis, function(e) {
  n_in <- sum(asg$epitope == e, na.rm = TRUE)
  jm_sensitivity(asg, e, n_in)
})
res <- list(
  n_clusters        = m$n_clusters,
  n_units_clustered = m$n_units_clustered,
  retention         = m$retention,
  purity            = m$purity,
  pct_clusters_hi   = m$pct_clusters_hi,
  pct_units_hi      = m$pct_units_hi,
  polyspecificity   = m$polyspecificity,
  cluster_sizes     = as.list(setNames(ct$size, ct$cluster_id)),
  cluster_purities  = as.list(setNames(ct$purity, ct$cluster_id)),
  sensitivity       = as.list(setNames(as.numeric(sens), epis))
)
writeLines(jsonlite::toJSON(res, auto_unbox = TRUE, digits = 12, na = "null"), out)
'''


def random_assignments(rng: np.random.Generator, n_clusters: int,
                       n_epitopes: int) -> pd.DataFrame:
    """Random clusters with duplicate receptors, singletons, and unlabelled
    rows — the awkward shapes where the two implementations could diverge."""
    rows = []
    rid = 0
    for c in range(n_clusters):
        size = int(rng.integers(1, 9))
        for _ in range(size):
            # ~15% of the time reuse an earlier receptor, so some receptors sit
            # in two clusters; that is what separates row counts from distinct
            # counts.
            if rows and rng.random() < 0.15:
                r = rows[int(rng.integers(0, len(rows)))][1]
            else:
                r = f"R{rid:04d}"
                rid += 1
            epi = (None if rng.random() < 0.08
                   else f"EPI{int(rng.integers(0, n_epitopes))}")
            rows.append((f"C{c:03d}", r, epi))
    df = pd.DataFrame(rows, columns=["cluster_id", "receptor_id", "epitope"])
    # exact duplicate rows, which dedup must absorb
    if len(df) > 4:
        df = pd.concat([df, df.sample(3, random_state=int(rng.integers(1e6)))],
                       ignore_index=True)
    return df


def run_r(csv: str, denom: int) -> dict:
    with tempfile.NamedTemporaryFile("r", suffix=".R", delete=False) as f:
        drv = f.name
    with open(drv, "w") as fh:
        fh.write(R_DRIVER)
    out = csv + ".json"
    env = dict(os.environ, SPTCE_REPO=REPO)
    p = subprocess.run(["Rscript", drv, csv, str(denom), out],
                       capture_output=True, text=True, env=env)
    if not os.path.exists(out):
        raise RuntimeError(f"R driver failed:\n{p.stdout}\n{p.stderr}")
    with open(out) as fh:
        return json.load(fh)


def main() -> int:
    rng = np.random.default_rng(20261009)
    n_pass = n_fail = 0

    def chk(label, a, b, tol=1e-9):
        nonlocal n_pass, n_fail
        if a is None and b is None:
            same = True
        elif a is None or b is None:
            same = False
        elif isinstance(a, str) or isinstance(b, str):
            same = a == b
        else:
            a, b = float(a), float(b)
            same = (np.isnan(a) and np.isnan(b)) or abs(a - b) < tol
        if same:
            n_pass += 1
        else:
            n_fail += 1
            print(f"  FAIL  {label}: python={a!r} R={b!r}")
        return same

    for trial in range(6):
        n_clusters = int(rng.integers(6, 26))
        n_epitopes = int(rng.integers(2, 7))
        asg = random_assignments(rng, n_clusters, n_epitopes)
        denom = int(asg["receptor_id"].nunique() * 2)

        with tempfile.TemporaryDirectory() as td:
            csv = os.path.join(td, "asg.csv")
            asg.to_csv(csv, index=False, na_rep="NA")
            r = run_r(csv, denom)

            m = M.metrics(asg, denom=denom)
            ct = M.cluster_table(asg)

            print(f"\n-- trial {trial + 1}: {n_clusters} clusters, "
                  f"{len(asg)} rows, {asg['receptor_id'].nunique()} receptors, "
                  f"denom {denom}")
            chk("n_clusters", m.n_clusters, r["n_clusters"])
            chk("n_units_clustered", m.n_units_clustered, r["n_units_clustered"])
            chk("retention", m.retention, r["retention"])
            chk("purity", m.purity, r["purity"])
            chk("pct_clusters_hi", m.pct_clusters_hi, r["pct_clusters_hi"])
            chk("pct_units_hi", m.pct_units_hi, r["pct_units_hi"])
            chk("polyspecificity", m.polyspecificity, r["polyspecificity"])

            # per-cluster agreement, not just the aggregates
            rs = r["cluster_sizes"]
            rp = r["cluster_purities"]
            chk("n clusters kept", len(ct), len(rs))
            py_sizes = dict(zip(ct["cluster_id"], ct["size"]))
            py_pur = dict(zip(ct["cluster_id"], ct["purity"]))
            for cid in sorted(rs):
                chk(f"size[{cid}]", py_sizes.get(cid), rs[cid])
                chk(f"purity[{cid}]", py_pur.get(cid), rp[cid])

            for epi, rv in r["sensitivity"].items():
                n_in = int((asg["epitope"] == epi).sum())
                pv = M.sensitivity(asg, epi, n_in)
                chk(f"sensitivity[{epi}]", pv, rv)
                # and the ceiling invariant must hold in the Python version
                cf = M.clustered_fraction(asg, epi, n_in)
                if not np.isnan(pv) and not np.isnan(cf) and pv > cf + 1e-9:
                    n_fail += 1
                    print(f"  FAIL  sensitivity[{epi}] {pv} exceeds ceiling {cf}")
                else:
                    n_pass += 1

    print(f"\n==== metrics vs R: {n_pass} agreed, {n_fail} disagreed ====")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
