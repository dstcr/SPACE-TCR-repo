"""Jouannet benchmark metrics, in Python.

The Python counterpart of ``sptce_metrics.R``, which follows the reference in
the Jouannet repository (``Functions/purity_functions.R`` and
``Analysis/benchmark_tools_pooled_database_CD8.Rmd``). It lets a learned model
be scored with the same equations as SPACE-TCR; ``tests/test_metrics_match_r.py``
checks the two give the same results.

Definitions, with the reference line:

* a cluster qualifies at ``size >= 2``             (Rmd:229, ``th <- 1``, used as ``size > th``)
* high purity is ``>= 0.9``                        (Rmd:351)
* retention counts distinct receptors              (Rmd:217-226)
* a cluster qualifies for sensitivity when ``size > 3`` and
  ``n_e / size >= 0.67``                           (Rmd:1547)
* the sensitivity numerator counts distinct receptors (Rmd:1575, 1578)

Everything takes a long assignment table with columns ``cluster_id``, a unit
column naming receptor identity, and ``epitope``.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd

MIN_CLUSTER_SIZE = 2
HIGH_PURITY_CUT = 0.9
SENS_MIN_SIZE = 4       # the reference writes `size > 3`
SENS_MIN_FRAC = 0.67    # the reference's literal 0.67, NOT 2/3


def cluster_table(asg: pd.DataFrame, unit: str = "receptor_id",
                  min_size: int = MIN_CLUSTER_SIZE) -> pd.DataFrame:
    """Per-cluster size, dominant-epitope count, label count and purity.

    Deduplicated on (cluster_id, unit) first, so a receptor listed twice in one
    cluster cannot inflate its size — the reference applies ``unique()`` after
    its merge for the same reason.
    """
    for col in ("cluster_id", unit):
        if col not in asg.columns:
            raise KeyError(f"assignments table needs a {col!r} column")
    d = asg[asg[unit].notna()].drop_duplicates(subset=["cluster_id", unit])
    if d.empty:
        return pd.DataFrame(columns=["cluster_id", "size", "n_dom", "n_labels",
                                     "dom_epitope", "purity"])

    sizes = d.groupby("cluster_id", sort=False).size().rename("size")

    lab = d[d["epitope"].notna()]
    if lab.empty:
        out = sizes.to_frame().reset_index()
        out["n_dom"] = 0
        out["n_labels"] = 0
        out["dom_epitope"] = None
    else:
        per = lab.groupby(["cluster_id", "epitope"], sort=False).size().rename("n_epi")
        per = per.reset_index()
        idx = per.groupby("cluster_id", sort=False)["n_epi"].idxmax()
        dom = per.loc[idx, ["cluster_id", "epitope", "n_epi"]].rename(
            columns={"epitope": "dom_epitope", "n_epi": "n_dom"})
        nlab = per.groupby("cluster_id", sort=False)["epitope"].nunique().rename("n_labels")
        out = (sizes.to_frame().reset_index()
               .merge(dom, on="cluster_id", how="left")
               .merge(nlab.reset_index(), on="cluster_id", how="left"))
        out["n_dom"] = out["n_dom"].fillna(0).astype(int)
        out["n_labels"] = out["n_labels"].fillna(0).astype(int)

    out["purity"] = np.where(out["size"] > 0, out["n_dom"] / out["size"], np.nan)
    return out[out["size"] >= min_size].reset_index(drop=True)


def purity(ct: pd.DataFrame) -> float:
    """Global purity: sum(n_dom) / sum(size). Reference: ``purity_function``."""
    if ct.empty:
        return float("nan")
    return float(ct["n_dom"].sum() / ct["size"].sum())


def polyspecificity(ct: pd.DataFrame) -> float:
    """Fraction of non-singleton clusters carrying more than one epitope."""
    if ct.empty:
        return float("nan")
    return float((ct["n_labels"] > 1).mean())


@dataclass
class Metrics:
    n_clusters: int
    n_units_clustered: int
    n_rows_clustered: int
    retention: float
    purity: float
    pct_clusters_hi: float
    pct_units_hi: float
    polyspecificity: float
    unit: str
    denom: int

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def metrics(asg: pd.DataFrame, denom: int, unit: str = "receptor_id",
            min_size: int = MIN_CLUSTER_SIZE,
            hi: float = HIGH_PURITY_CUT) -> Metrics:
    """The headline metric row for one clustering.

    ``retention`` is distinct receptors in non-singleton clusters over
    ``denom``, as in the reference; a receptor in several clusters counts once.
    """
    ct = cluster_table(asg, unit=unit, min_size=min_size)
    if ct.empty:
        return Metrics(0, 0, 0, 0.0, float("nan"), float("nan"),
                       float("nan"), float("nan"), unit, denom)
    keep = asg[asg["cluster_id"].isin(set(ct["cluster_id"])) & asg[unit].notna()]
    n_units = int(keep[unit].nunique())
    hi_ct = ct[ct["purity"] >= hi]
    hi_units = int(keep[keep["cluster_id"].isin(set(hi_ct["cluster_id"]))][unit].nunique())
    return Metrics(
        n_clusters=int(len(ct)),
        n_units_clustered=n_units,
        n_rows_clustered=int(len(keep)),
        retention=float(n_units / denom) if denom else float("nan"),
        purity=purity(ct),
        pct_clusters_hi=float(len(hi_ct) / len(ct)),
        pct_units_hi=float(hi_units / n_units) if n_units else float("nan"),
        polyspecificity=polyspecificity(ct),
        unit=unit,
        denom=int(denom),
    )


def sensitivity(asg: pd.DataFrame, epi: str, n_input_e: int,
                unit: str = "receptor_id") -> float:
    """Epitope-level sensitivity. Reference: Rmd:1547 and Rmd:1575/1578.

    A cluster qualifies for epitope ``e`` when its dominant epitope is ``e``, its
    size exceeds 3, and at least 67% of its members carry ``e``. The numerator is
    the number of DISTINCT receptors labelled ``e`` in those clusters; summing
    per-cluster counts double-counts a receptor in two qualifying clusters and
    can push the result above the clustered fraction.
    """
    if not n_input_e or n_input_e <= 0:
        return float("nan")
    ct = cluster_table(asg, unit=unit)
    if ct.empty:
        return 0.0
    d = asg[asg[unit].notna() & asg["cluster_id"].isin(set(ct["cluster_id"]))]
    d = d.drop_duplicates(subset=["cluster_id", unit])
    n_e = (d[d["epitope"] == epi].groupby("cluster_id", sort=False)
           .size().rename("n_e").reset_index())
    q = ct.merge(n_e, on="cluster_id", how="left")
    q["n_e"] = q["n_e"].fillna(0).astype(int)
    q["frac"] = q["n_e"] / q["size"]
    q = q[(q["dom_epitope"] == epi) & (q["size"] >= SENS_MIN_SIZE)
          & (q["frac"] >= SENS_MIN_FRAC)]
    if q.empty:
        return 0.0
    rec = asg[asg["cluster_id"].isin(set(q["cluster_id"]))
              & (asg["epitope"] == epi) & asg[unit].notna()][unit].nunique()
    return float(rec / n_input_e)


def clustered_fraction(asg: pd.DataFrame, epi: str, n_input_e: int,
                       unit: str = "receptor_id") -> float:
    """Ceiling for ``sensitivity``: the share of this epitope's input that is
    clustered at all. Sensitivity above this is impossible by construction."""
    if not n_input_e:
        return float("nan")
    ct = cluster_table(asg, unit=unit)
    d = asg[asg["cluster_id"].isin(set(ct["cluster_id"]))
            & (asg["epitope"] == epi) & asg[unit].notna()]
    return float(d[unit].nunique() / n_input_e)


def check_metrics(m: Metrics, label: str = "") -> None:
    """Check the metric invariants; raises if any is violated."""
    bad = []
    if m.n_units_clustered > m.denom:
        bad.append("units clustered exceeds the input denominator (modes pooled?)")
    if m.retention > 1 + 1e-9:
        bad.append("retention above 1")
    if not np.isnan(m.purity) and m.purity > 1 + 1e-9:
        bad.append("purity above 1")
    if not np.isnan(m.pct_clusters_hi) and m.pct_clusters_hi > 1 + 1e-9:
        bad.append("high-purity cluster share above 1")
    if bad:
        tag = f" [{label}]" if label else ""
        raise AssertionError("check_metrics failed" + tag + ":\n  - " + "\n  - ".join(bad))


def check_sensitivity(df: pd.DataFrame, label: str = "") -> None:
    """Sensitivity cannot exceed the clustered fraction. This is the invariant
    the published GILGFVFTL value (0.83 against a 0.311 ceiling) broke."""
    need = {"sensitivity", "clustered_frac"}
    if not need.issubset(df.columns):
        raise KeyError(f"need columns {need}")
    bad = df[df["sensitivity"].notna() & df["clustered_frac"].notna()
             & (df["sensitivity"] > df["clustered_frac"] + 1e-9)]
    if len(bad):
        tag = f" [{label}]" if label else ""
        raise AssertionError(
            "check_sensitivity failed" + tag + ": sensitivity exceeds the "
            "clustered fraction for:\n" + bad.to_string(index=False))


def within_cluster_median_ld(asg: pd.DataFrame, seq_col: str = "cdr3b"
                             ) -> pd.DataFrame:
    """Median pairwise Levenshtein distance among each cluster's distinct
    sequences, used for the divergence classes of Table 3."""
    from rapidfuzz.distance import Levenshtein
    d = asg[asg[seq_col].notna()].drop_duplicates(subset=["cluster_id", seq_col])
    rows = []
    for cid, g in d.groupby("cluster_id", sort=False):
        seqs = g[seq_col].tolist()
        if len(seqs) < 2:
            rows.append((cid, len(seqs), 0.0))
            continue
        ds = [Levenshtein.distance(seqs[i], seqs[j])
              for i in range(len(seqs)) for j in range(i + 1, len(seqs))]
        rows.append((cid, len(seqs), float(np.median(ds))))
    return pd.DataFrame(rows, columns=["cluster_id", "n_seq", "median_ld"])


LD_CLASS_BREAKS = [(-np.inf, 1, "Sequence-identical / near-identical"),
                   (1, 3, "Sequence-similar"),
                   (3, 6, "Sequence-diverse"),
                   (6, np.inf, "Sequence-unrelated")]


def ld_class(x: float) -> str:
    """Table 3's divergence classes, by median within-cluster LD."""
    for lo, hi, name in LD_CLASS_BREAKS:
        if lo < x <= hi:
            return name
    return LD_CLASS_BREAKS[-1][2]


def sequence_reachable(cls: str) -> bool:
    """Whether a divergence class is within reach of the benchmarked sequence
    methods. The two that are not are the paper's novelty claim."""
    return cls in ("Sequence-identical / near-identical", "Sequence-similar")
