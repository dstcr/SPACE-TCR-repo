"""Distil SPACE-TCR's structural clustering into a sequence-only model.

WHAT THIS IS FOR
----------------
SPACE-TCR groups receptors by predicted CDR-loop geometry, so every receptor
must first be modelled with TCRBuilder2, which limits throughput (§4 of the
paper).

This trains a siamese encoder on sequence alone, supervised by the structural
clusters SPACE-TCR produced, so that structural co-membership could be
predicted for receptors that were never modelled.

The baseline is Levenshtein distance on CDR3b at matched retention, reported
alongside every result.

WHAT IT MEASURES
----------------
1. pair-level ranking (AUROC, average precision) on held-out pairs;
2. cluster-level metrics after clustering the learned embedding, scored with the
   Jouannet equations used for Tables 1 and 2, via ``ml/metrics.py`` (checked
   against ``sptce_metrics.R`` by ``ml/tests/test_metrics_match_r.py``);
3. recovery restricted to sequence-diverse clusters (median within-cluster
   CDR3b Levenshtein distance >= 4), the 22.7% of clusters the paper describes
   as out of reach of sequence methods. Recovering only near-identical clusters
   would mean the model has learned edit distance.

SPLITS
------
Grouped by epitope, not by receptor: receptors of the same epitope are often
near-clonal, so a receptor-level split would let the model memorise epitopes.
Held-out epitopes test unseen specificities. CDR3b sequences present on both
sides of the split are dropped from the test side and counted.

USAGE
-----
    python3 ml/distill.py --assignments cluster_assignments.tsv \
                          --labels data/jouannet_labelled_input.tsv \
                          --out ml_out

    python3 ml/distill.py --synthetic        # self-test, no data needed

Needs only the cluster assignments and the labelled input; the distance tables
are not required.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metrics as M  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i + 1 for i, a in enumerate(AA)}   # 0 is the pad index
MAX_LEN = 30


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def encode_seq(s: str, max_len: int = MAX_LEN) -> np.ndarray:
    s = (s or "")[:max_len]
    v = np.zeros(max_len, dtype=np.int64)
    for i, ch in enumerate(s):
        v[i] = AA_IDX.get(ch, 0)
    return v


def strip_allele(x) -> str:
    if not isinstance(x, str):
        return "UNK"
    return x.split("*")[0] if "*" in x else x


@dataclass
class Cohort:
    df: pd.DataFrame                 # receptor_id, cdr3b, cdr3a, epitope, v_gene, j_gene, cluster_id, bin_key
    gene_vocab: Dict[str, int]

    @property
    def n(self) -> int:
        return len(self.df)


def _pick(df: pd.DataFrame, cands: Sequence[str]) -> Optional[str]:
    for c in cands:
        if c in df.columns:
            return c
    return None


def load_cohort(assignments: str, labels: str) -> Cohort:
    """Join SPACE-TCR cluster assignments to the labelled sequence table."""
    asg = pd.read_csv(assignments, sep=None, engine="python")
    lab = pd.read_csv(labels, sep=None, engine="python")

    a_id = _pick(asg, ["tcr_id", "id", "receptor_id", "ID"])
    a_cl = _pick(asg, ["cluster_id", "cluster", "cluster_by_rmsd", "clusters"])
    a_bin = _pick(asg, ["bin_key", "cluster_by_length", "bin"])
    if a_id is None or a_cl is None:
        raise SystemExit(f"assignments need an id and a cluster column; got {list(asg.columns)}")

    l_id = _pick(lab, ["id", "tcr_id", "receptor_id"])
    l_b = _pick(lab, ["cdr3_aa_beta", "cdr3b", "CDR3b", "cdr3_beta"])
    l_a = _pick(lab, ["cdr3_aa_alpha", "cdr3a", "CDR3a", "cdr3_alpha"])
    l_e = _pick(lab, ["epitope", "Epitope"])
    l_v = _pick(lab, ["v_b_gene", "v_call_beta", "V_beta", "TRBV"])
    l_j = _pick(lab, ["j_b_gene", "j_call_beta", "J_beta", "TRBJ"])
    if l_id is None or l_b is None or l_e is None:
        raise SystemExit(f"labels need id, CDR3b and epitope; got {list(lab.columns)}")

    L = pd.DataFrame({
        "receptor_id": lab[l_id].astype(str),
        "cdr3b": lab[l_b].astype(str),
        "cdr3a": lab[l_a].astype(str) if l_a else "",
        "epitope": lab[l_e].astype(str),
        "v_gene": lab[l_v].map(strip_allele) if l_v else "UNK",
        "j_gene": lab[l_j].map(strip_allele) if l_j else "UNK",
    })
    A = pd.DataFrame({
        "receptor_id": asg[a_id].astype(str).str.replace(r"\.pdb$", "", regex=True),
        "cluster_id": asg[a_cl].astype(str),
        "bin_key": asg[a_bin].astype(str) if a_bin else "nobin",
    })
    df = A.merge(L, on="receptor_id", how="inner")
    df = df[df["epitope"].notna() & (df["epitope"] != "nan") & (df["epitope"] != "")]
    # keep only non-singleton clusters: a singleton carries no pair supervision
    sizes = df.groupby("cluster_id")["receptor_id"].transform("nunique")
    df = df[sizes >= 2].drop_duplicates(subset=["receptor_id", "cluster_id"])
    df = df.reset_index(drop=True)
    genes = sorted(set(df["v_gene"]) | set(df["j_gene"]) | {"UNK"})
    return Cohort(df=df, gene_vocab={g: i for i, g in enumerate(genes)})


AA_GROUPS = [list("AVLIMFWC"), list("STNQGP"), list("KRH"), list("DE")]


def synthetic_cohort(n_clusters: int = 40, per_cluster: int = 14,
                     length: int = 14, frac_beta_signal: float = 0.35,
                     seed: int = 11) -> Cohort:
    """A fixture that separates "learned structure" from "learned edit distance".

    Modelled on the paper's own section 3.6: SPACE-TCR clusters are coherent in
    CDR3-alpha while remaining diverse in CDR3-beta, which is precisely what a
    CDR3-beta-only method cannot see. So most clusters here carry their signal in
    the **alpha chain and the V/J genes**, with CDR3-beta drawn at random. A
    Levenshtein baseline on CDR3-beta is then at chance *by construction*, not by
    luck, while the signal remains learnable from the paired input.

    Every receptor shares one CDR3 length, so length cannot leak the cluster
    either. A minority (``frac_beta_signal``) additionally carry a beta motif,
    so the evaluation covers the regime edit distance *can* reach rather than
    only the flattering one.

    A harness that cannot beat Levenshtein here is broken. One that can is
    measuring what it claims to.
    """
    rng = np.random.default_rng(seed)
    n_beta = int(round(n_clusters * frac_beta_signal))
    vs = [f"TRBV{i}" for i in range(1, 13)]
    js = [f"TRBJ{i}-1" for i in range(1, 7)]
    rows = []
    for c in range(n_clusters):
        epi = f"EPI{c:03d}"
        beta_signal = c < n_beta
        # alpha motif: three fixed residues at three fixed positions
        a_pos = rng.choice(np.arange(2, 11), size=3, replace=False)
        a_aa = rng.choice(list(AA), size=3)
        v_c, j_c = vs[c % len(vs)], js[c % len(js)]
        b_pos = rng.choice(np.arange(2, 11), size=3, replace=False)
        b_aa = rng.choice(list(AA), size=3)
        for _ in range(per_cluster):
            beta = list(rng.choice(list(AA), size=length))
            if beta_signal:
                for p_, x_ in zip(b_pos, b_aa):
                    beta[p_] = x_
            alpha = list(rng.choice(list(AA), size=length))
            for p_, x_ in zip(a_pos, a_aa):
                alpha[p_] = x_
            rows.append((f"C{c:03d}", "".join(beta), "".join(alpha), epi,
                         v_c, j_c, "beta+alpha" if beta_signal else "alpha-only"))
    df = pd.DataFrame(rows, columns=["cluster_id", "cdr3b", "cdr3a", "epitope",
                                     "v_gene", "j_gene", "kind"])
    df["receptor_id"] = [f"TCR{i:05d}" for i in range(len(df))]
    df["bin_key"] = str(length)            # one bin: length cannot separate
    df = df[["receptor_id", "cluster_id", "bin_key", "cdr3b", "cdr3a",
             "epitope", "v_gene", "j_gene", "kind"]].reset_index(drop=True)
    genes = sorted(set(df["v_gene"]) | set(df["j_gene"]) | {"UNK"})
    return Cohort(df=df, gene_vocab={g: i for i, g in enumerate(genes)})


# ---------------------------------------------------------------------------
# splits
# ---------------------------------------------------------------------------
def epitope_split(df: pd.DataFrame, frac_test: float = 0.3, seed: int = 0
                  ) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Hold out whole epitopes. Also drop test receptors whose CDR3b appears in
    training, and report how many."""
    rng = np.random.default_rng(seed)
    epis = df["epitope"].drop_duplicates().tolist()
    rng.shuffle(epis)
    n_test = max(1, int(round(len(epis) * frac_test)))
    test_epis = set(epis[:n_test])
    tr = df[~df["epitope"].isin(test_epis)].copy()
    te = df[df["epitope"].isin(test_epis)].copy()
    shared = set(tr["cdr3b"]) & set(te["cdr3b"])
    n_dropped = int(te["cdr3b"].isin(shared).sum())
    te = te[~te["cdr3b"].isin(shared)].copy()
    info = {"n_epitopes_total": len(epis), "n_epitopes_test": len(test_epis),
            "test_epitopes": sorted(test_epis),
            "n_train": int(len(tr)), "n_test": int(len(te)),
            "n_test_dropped_shared_cdr3b": n_dropped}
    return tr, te, info


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------
def build_model(n_genes: int, emb_dim: int = 64, seed: int = 0):
    import torch
    import torch.nn as nn

    torch.manual_seed(seed)

    class Encoder(nn.Module):
        """Small by design: with a few thousand receptors a large model would
        memorise the training epitopes."""

        def __init__(self):
            super().__init__()
            self.aa = nn.Embedding(len(AA) + 1, 24, padding_idx=0)
            self.gene = nn.Embedding(n_genes, 8)
            convs = []
            for k in (3, 5, 7):
                convs.append(nn.Sequential(
                    nn.Conv1d(24, 48, kernel_size=k, padding=k // 2),
                    nn.ReLU(),
                ))
            self.convs = nn.ModuleList(convs)
            self.head = nn.Sequential(
                nn.Linear(48 * 3 * 2 + 8 * 2, 128), nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(128, emb_dim),
            )

        def _chain(self, x):
            h = self.aa(x).transpose(1, 2)              # B,24,L
            feats = [c(h).amax(dim=2) for c in self.convs]
            return torch.cat(feats, dim=1)              # B,144

        def forward(self, b, a, v, j):
            hb = self._chain(b)
            ha = self._chain(a)
            g = torch.cat([self.gene(v), self.gene(j)], dim=1)
            z = self.head(torch.cat([hb, ha, g], dim=1))
            return nn.functional.normalize(z, dim=1)     # unit sphere

    return Encoder()


def tensors_for(df: pd.DataFrame, gene_vocab: Dict[str, int]):
    import torch
    b = np.stack([encode_seq(s) for s in df["cdr3b"]])
    a = np.stack([encode_seq(s) for s in df["cdr3a"]])
    v = np.array([gene_vocab.get(g, gene_vocab.get("UNK", 0)) for g in df["v_gene"]])
    j = np.array([gene_vocab.get(g, gene_vocab.get("UNK", 0)) for g in df["j_gene"]])
    return (torch.as_tensor(b), torch.as_tensor(a),
            torch.as_tensor(v), torch.as_tensor(j))


def sample_pairs(df: pd.DataFrame, n_pairs: int, rng: np.random.Generator,
                 hard_frac: float = 0.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Positives share a structural cluster. Negatives are drawn preferentially
    from the same CDR-length bin, because SPACE-TCR only compares receptors
    within a bin; a negative from another bin would be an easier distinction
    than the structural method ever makes."""
    idx_by_cluster: Dict[str, List[int]] = {}
    for i, c in enumerate(df["cluster_id"].to_numpy()):
        idx_by_cluster.setdefault(c, []).append(i)
    idx_by_bin: Dict[str, List[int]] = {}
    for i, b in enumerate(df["bin_key"].to_numpy()):
        idx_by_bin.setdefault(b, []).append(i)
    clusters = [c for c, v in idx_by_cluster.items() if len(v) >= 2]
    if not clusters:
        raise SystemExit("no multi-member clusters to learn from")
    cl_of = df["cluster_id"].to_numpy()

    ii, jj, yy = [], [], []
    half = n_pairs // 2
    for _ in range(half):
        c = clusters[int(rng.integers(len(clusters)))]
        a, b = rng.choice(idx_by_cluster[c], size=2, replace=False)
        ii.append(a); jj.append(b); yy.append(1)
    for _ in range(n_pairs - half):
        if rng.random() < hard_frac:
            bk = list(idx_by_bin)[int(rng.integers(len(idx_by_bin)))]
            pool = idx_by_bin[bk]
            if len(pool) < 2:
                pool = range(len(df))
            a, b = rng.choice(list(pool), size=2, replace=True)
        else:
            a, b = rng.integers(0, len(df), size=2)
        if cl_of[a] == cl_of[b]:
            continue
        ii.append(int(a)); jj.append(int(b)); yy.append(0)
    return np.array(ii), np.array(jj), np.array(yy, dtype=np.float32)


# Minimum median inner-validation pair distance for an epoch to be eligible for
# best-epoch restore. Below this the embedding has collapsed to a point and no
# clustering threshold can separate anything; see the note in train().
MIN_SPREAD = 0.10


def train(cohort: Cohort, tr: pd.DataFrame, epochs: int = 30,
          pairs_per_epoch: int = 20000, margin: float = 0.7,
          lr: float = 1e-3, emb_dim: int = 64, seed: int = 0,
          quiet: bool = False):
    import torch
    import torch.nn as nn

    rng = np.random.default_rng(seed)
    model = build_model(len(cohort.gene_vocab), emb_dim=emb_dim, seed=seed)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    B, A, V, J = tensors_for(tr, cohort.gene_vocab)

    # a small slice of training epitopes held back to watch for overfitting
    inner_epis = tr["epitope"].drop_duplicates().tolist()
    rng.shuffle(inner_epis)
    val_epis = set(inner_epis[:max(1, len(inner_epis) // 8)])
    is_val = tr["epitope"].isin(val_epis).to_numpy()

    hist = []
    best_auc, best_ep, best_state = float("nan"), epochs, None
    best_spread = float("nan")
    for ep in range(epochs):
        model.train()
        ii, jj, yy = sample_pairs(tr, pairs_per_epoch, rng)
        keep = (~is_val[ii]) & (~is_val[jj])
        ii, jj, yy = ii[keep], jj[keep], yy[keep]
        perm = rng.permutation(len(ii))
        ii, jj, yy = ii[perm], jj[perm], yy[perm]
        tot, nb = 0.0, 0
        for s in range(0, len(ii), 256):
            a_i = ii[s:s + 256]; b_i = jj[s:s + 256]
            y = torch.as_tensor(yy[s:s + 256])
            za = model(B[a_i], A[a_i], V[a_i], J[a_i])
            zb = model(B[b_i], A[b_i], V[b_i], J[b_i])
            d = (za - zb).pow(2).sum(1).clamp_min(1e-12).sqrt()
            # contrastive: pull positives to 0, push negatives past the margin
            loss = (y * d.pow(2) + (1 - y) * nn.functional.relu(margin - d).pow(2)).mean()
            opt.zero_grad(); loss.backward(); opt.step()
            tot += loss.detach().item(); nb += 1
        vi, vj, vy = sample_pairs(tr[is_val].reset_index(drop=True),
                                  2000, rng) if is_val.sum() > 10 else (None, None, None)
        vauc = float("nan")
        vspread = float("nan")
        if vi is not None and len(vi):
            sub = tr[is_val].reset_index(drop=True)
            Bv, Av, Vv, Jv = tensors_for(sub, cohort.gene_vocab)
            model.eval()
            with torch.no_grad():
                za = model(Bv[vi], Av[vi], Vv[vi], Jv[vi])
                zb = model(Bv[vj], Av[vj], Vv[vj], Jv[vj])
                d = (za - zb).pow(2).sum(1).sqrt().numpy()
            vauc = auroc(-d, vy)
            # Spread of the embedding: the median inner-validation pair
            # distance. A contrastive encoder can collapse towards a point and
            # still rank pairs well, since the distances keep their order while
            # shrinking; a collapsed embedding then puts everything in one
            # cluster at any threshold. AUROC does not show this, so the spread
            # is tracked separately and used in model selection below.
            vspread = float(np.median(d)) if len(d) else float("nan")
        hist.append({"epoch": ep + 1, "loss": tot / max(nb, 1), "val_auroc": vauc,
                     "val_spread": vspread})

        # Keep the best weights by inner-validation AUROC, among epochs whose
        # embedding has not collapsed. On the 874-receptor cohort AUROC peaks
        # at epoch 1-2, where the maximum pairwise distance on the training
        # receptors is 0.065: every threshold in the evaluation sweep (0.05
        # upward) then merges everything. Training to epoch 40 gives a lower
        # AUROC (0.656 vs 0.684) but a clustering that responds to the
        # threshold (40 clusters at 33.1% retention / 99.1% purity, rising to
        # 79 at 81.3% / 93.1%). MIN_SPREAD excludes collapsed epochs; if none
        # qualifies, the final epoch is kept and a warning is printed.
        ok_spread = vspread == vspread and vspread >= MIN_SPREAD
        if ok_spread and vauc == vauc and (best_auc != best_auc or vauc > best_auc):
            best_auc, best_ep, best_spread = vauc, ep + 1, vspread
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}

        if not quiet and (ep % 5 == 4 or ep == 0):
            print(f"    epoch {ep+1:3d}  loss {tot/max(nb,1):.4f}  "
                  f"inner-val AUROC {vauc:.3f}")

    if best_state is not None and best_ep != epochs:
        model.load_state_dict(best_state)
        if not quiet:
            print(f"    restored epoch {best_ep} (inner-val AUROC {best_auc:.3f}, "
                  f"spread {best_spread:.3f}); last epoch scored "
                  f"{hist[-1]['val_auroc']:.3f}")
    elif best_state is None and not quiet:
        sp = hist[-1].get("val_spread", float("nan")) if hist else float("nan")
        print(f"    WARNING: no epoch cleared the collapse floor "
              f"(spread >= {MIN_SPREAD}); keeping the final epoch "
              f"(spread {sp:.3f}). The embedding may be degenerate, in which "
              f"case every clustering threshold will merge everything.")
    return model, hist


def embed(model, df: pd.DataFrame, gene_vocab) -> np.ndarray:
    import torch
    B, A, V, J = tensors_for(df, gene_vocab)
    model.eval()
    out = []
    with torch.no_grad():
        for s in range(0, len(df), 512):
            out.append(model(B[s:s+512], A[s:s+512], V[s:s+512], J[s:s+512]).numpy())
    return np.vstack(out) if out else np.zeros((0, 1))


# ---------------------------------------------------------------------------
# evaluation
# ---------------------------------------------------------------------------
def auroc(score: np.ndarray, y: np.ndarray) -> float:
    """Rank-based AUROC; ties averaged. No sklearn dependency at call sites."""
    if len(np.unique(y)) < 2:
        return float("nan")
    order = np.argsort(score)
    ranks = np.empty(len(score), dtype=float)
    ranks[order] = np.arange(1, len(score) + 1)
    # average ranks within ties
    s_sorted = score[order]
    i = 0
    while i < len(s_sorted):
        j = i
        while j + 1 < len(s_sorted) and s_sorted[j + 1] == s_sorted[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (ranks[order[i]] + ranks[order[j]]) / 2
        i = j + 1
    n_pos = float((y == 1).sum()); n_neg = float((y == 0).sum())
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(score: np.ndarray, y: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    order = np.argsort(-score)
    y = y[order]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / max(y.sum(), 1))


def cluster_by_threshold(ids: Sequence[str], dist_fn, threshold: float,
                         bins: Optional[Sequence[str]] = None) -> pd.Series:
    """Connected components of the distance-threshold graph — the same procedure
    SPACE-TCR uses, so the comparison is like for like. Restricted within bins
    when given, again mirroring SPACE-TCR."""
    n = len(ids)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    groups: Dict[str, List[int]] = {}
    for i, b in enumerate(bins if bins is not None else ["all"] * n):
        groups.setdefault(b, []).append(i)
    for g in groups.values():
        for a_pos in range(len(g)):
            for b_pos in range(a_pos + 1, len(g)):
                i, j = g[a_pos], g[b_pos]
                if dist_fn(i, j) <= threshold:
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        parent[ri] = rj
    return pd.Series([f"P{find(i)}" for i in range(n)], index=list(ids))


def score_clustering(df: pd.DataFrame, cluster_ids: pd.Series, denom: int,
                     label: str) -> dict:
    asg = pd.DataFrame({
        "cluster_id": cluster_ids.reindex(df["receptor_id"]).to_numpy(),
        "receptor_id": df["receptor_id"].to_numpy(),
        "epitope": df["epitope"].to_numpy(),
        "cdr3b": df["cdr3b"].to_numpy(),
    })
    m = M.metrics(asg, denom=denom)
    M.check_metrics(m, label)
    out = m.as_dict()
    out["label"] = label
    # how divergent are the clusters it forms?
    ld = M.within_cluster_median_ld(asg, "cdr3b")
    if len(ld):
        ct = M.cluster_table(asg)
        hi = ct[ct["purity"] >= M.HIGH_PURITY_CUT]["cluster_id"]
        ldh = ld[ld["cluster_id"].isin(set(hi))]
        if len(ldh):
            cls = ldh["median_ld"].map(M.ld_class)
            unreachable = (~cls.map(M.sequence_reachable)).mean()
            out["pct_hi_clusters_sequence_diverse"] = float(unreachable)
            out["median_ld_of_hi_clusters"] = float(ldh["median_ld"].median())
    return out


def ld_distance_fn(seqs: Sequence[str]):
    from rapidfuzz.distance import Levenshtein
    def f(i, j):
        return Levenshtein.distance(seqs[i], seqs[j])
    return f


def calibrate_threshold(dist_fn, ids, bins, target_retention: float,
                        denom: int, lo: float, hi: float,
                        steps: int = 12) -> Tuple[float, pd.Series]:
    """Pick the threshold whose retention is closest to a target, so two methods
    are compared at matched yield rather than at arbitrary operating points."""
    best = (None, None, 1e9)
    for t in np.linspace(lo, hi, steps):
        cl = cluster_by_threshold(ids, dist_fn, float(t), bins)
        sizes = cl.value_counts()
        kept = cl.isin(sizes[sizes >= 2].index)
        ret = kept.sum() / denom
        gap = abs(ret - target_retention)
        if gap < best[2]:
            best = (float(t), cl, gap)
    return best[0], best[1]


def evaluate(cohort: Cohort, model, tr: pd.DataFrame, te: pd.DataFrame,
             seed: int = 0, quiet: bool = False) -> dict:
    rng = np.random.default_rng(seed + 1)
    res: dict = {}

    # ---- 1. pair-level ranking on held-out epitopes ----------------------
    ii, jj, yy = sample_pairs(te.reset_index(drop=True), 20000, rng)
    Z = embed(model, te, cohort.gene_vocab)
    d_model = np.linalg.norm(Z[ii] - Z[jj], axis=1)
    seqs = te["cdr3b"].tolist()
    from rapidfuzz.distance import Levenshtein
    d_ld = np.array([Levenshtein.distance(seqs[a], seqs[b]) for a, b in zip(ii, jj)],
                    dtype=float)
    res["pair_level"] = {
        "n_pairs": int(len(ii)),
        "n_positive": int(yy.sum()),
        "model_auroc": auroc(-d_model, yy),
        "model_ap": average_precision(-d_model, yy),
        "levenshtein_auroc": auroc(-d_ld, yy),
        "levenshtein_ap": average_precision(-d_ld, yy),
    }

    # the pairs that matter: positives whose CDR3b are far apart
    far = (yy == 1) & (d_ld >= 4)
    near_neg = (yy == 0)
    if far.sum() > 10:
        sel = far | near_neg
        res["pair_level_sequence_diverse_only"] = {
            "n_divergent_positives": int(far.sum()),
            "model_auroc": auroc(-d_model[sel], yy[sel]),
            "levenshtein_auroc": auroc(-d_ld[sel], yy[sel]),
        }

    # ---- 2. cluster-level, swept like Table 1 ----------------------------
    ids = te["receptor_id"].tolist()
    bins = te["bin_key"].tolist()
    denom = len(te)
    truth = pd.Series(te["cluster_id"].to_numpy(), index=ids)
    res["structural_truth"] = score_clustering(te, truth, denom, "SPACE-TCR (truth)")
    truth_purity = res["structural_truth"]["purity"]

    def model_dist(i, j):
        return float(np.linalg.norm(Z[i] - Z[j]))

    # A single "matched retention" point is degenerate when the truth retains
    # everything: matching it forces a threshold loose enough to merge every
    # cluster together. Sweep instead, exactly as Table 1 sweeps RMSD, and let
    # the purity-retention curve speak. The operating point reported as "best"
    # is the loosest threshold still at or above the truth's purity, which is
    # what a user would actually pick.
    ld_fn = ld_distance_fn(seqs)

    def sweep(dist_fn, lo, hi, steps, name):
        rows, best = [], None
        for t in np.linspace(lo, hi, steps):
            cl = cluster_by_threshold(ids, dist_fn, float(t), bins)
            sc = score_clustering(te, cl, denom, name)
            sc["threshold"] = float(t)
            rows.append(sc)
            if sc["n_clusters"] and not np.isnan(sc["purity"]) \
               and sc["purity"] >= truth_purity - 1e-9:
                if best is None or sc["retention"] > best["retention"]:
                    best = sc
        fell_back = best is None
        if fell_back:
            cand = [r for r in rows if r["n_clusters"] and not np.isnan(r["purity"])]
            best = max(cand, key=lambda r: r["purity"]) if cand else rows[0]
        best = dict(best)
        best["reached_truth_purity"] = not fell_back
        return rows, best

    res["model_sweep"], m_best = sweep(model_dist, 0.05, 1.6, 12, "distilled model")
    res["model"] = dict(m_best)
    res["levenshtein_sweep"], l_best = sweep(ld_fn, 1, 8, 8, "Levenshtein CDR3b")
    res["levenshtein"] = dict(l_best)

    cl_model = cluster_by_threshold(ids, model_dist, res["model"]["threshold"], bins)
    cl_ld = cluster_by_threshold(ids, ld_fn, res["levenshtein"]["threshold"], bins)


    # ---- 3. agreement with the structural clusters -----------------------
    def ari(a: Sequence[str], b: Sequence[str]) -> float:
        ct = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy()
        n = ct.sum()
        if n < 2:
            return float("nan")
        comb = lambda x: (x * (x - 1) / 2).sum()
        idx = comb(ct); ra = comb(ct.sum(1)); rb = comb(ct.sum(0))
        exp = ra * rb / (n * (n - 1) / 2); mx = (ra + rb) / 2
        return float((idx - exp) / (mx - exp)) if mx != exp else 0.0

    tr_t = truth.reindex(ids).to_numpy()
    res["agreement_with_structure"] = {
        "model_ari": ari(tr_t, cl_model.reindex(ids).to_numpy()),
        "levenshtein_ari": ari(tr_t, cl_ld.reindex(ids).to_numpy()),
    }
    return res


# ---------------------------------------------------------------------------
def report(res: dict, info: dict) -> None:
    p = res["pair_level"]
    print("\n--- splits ---")
    print(f"  epitopes: {info['n_epitopes_total']} total, "
          f"{info['n_epitopes_test']} held out")
    print(f"  receptors: {info['n_train']} train / {info['n_test']} test"
          f"   ({info['n_test_dropped_shared_cdr3b']} test rows dropped for "
          f"sharing a CDR3b with train)")

    print("\n--- 1. pair-level ranking on held-out epitopes ---")
    print(f"  {'':24s} {'AUROC':>8s} {'AP':>8s}")
    print(f"  {'distilled model':24s} {p['model_auroc']:8.3f} {p['model_ap']:8.3f}")
    print(f"  {'Levenshtein CDR3b':24s} {p['levenshtein_auroc']:8.3f} "
          f"{p['levenshtein_ap']:8.3f}")
    if "pair_level_sequence_diverse_only" in res:
        q = res["pair_level_sequence_diverse_only"]
        print(f"\n  restricted to sequence-divergent positives "
              f"(CDR3b LD >= 4, n={q['n_divergent_positives']}):")
        print(f"  {'distilled model':24s} {q['model_auroc']:8.3f}")
        print(f"  {'Levenshtein CDR3b':24s} {q['levenshtein_auroc']:8.3f}")
        print("  ^ this is the comparison that matters: the paper's claim is")
        print("    that sequence methods cannot reach these.")

    print("\n--- 2. cluster-level, Jouannet metrics (swept, as Table 1) ---")
    print(f"  {'':24s} {'thresh':>7s} {'clusters':>9s} {'retention':>10s} "
          f"{'purity':>8s} {'>=0.9':>7s}")
    r = res["structural_truth"]
    print(f"  {'SPACE-TCR (truth)':24s} {'-':>7s} {r['n_clusters']:9d} "
          f"{r['retention']*100:9.1f}% {r['purity']*100:7.1f}% "
          f"{r['pct_clusters_hi']*100:6.1f}%")
    for key, name in (("model_sweep", "distilled model"),
                      ("levenshtein_sweep", "Levenshtein CDR3b")):
        print(f"  -- {name} --")
        for sc in res[key]:
            if not sc["n_clusters"]:
                continue
            pur = 0.0 if sc["purity"] is None or np.isnan(sc["purity"]) else sc["purity"]
            hic = 0.0 if sc["pct_clusters_hi"] is None or np.isnan(sc["pct_clusters_hi"]) else sc["pct_clusters_hi"]
            print(f"  {'':24s} {sc['threshold']:7.2f} {sc['n_clusters']:9d} "
                  f"{sc['retention']*100:9.1f}% {pur*100:7.1f}% {hic*100:6.1f}%")
    print("\n  operating point: loosest threshold at or above the truth's purity")
    for k in ("model", "levenshtein"):
        sc = res[k]
        if not sc.get("reached_truth_purity", True):
            print(f"  {sc['label']:24s} NEVER reaches the truth's purity at any "
                  f"threshold; showing its best")
        pur = 0.0 if np.isnan(sc["purity"]) else sc["purity"]
        hic = 0.0 if np.isnan(sc["pct_clusters_hi"]) else sc["pct_clusters_hi"]
        print(f"  {sc['label']:24s} {sc['threshold']:7.2f} {sc['n_clusters']:9d} "
              f"{sc['retention']*100:9.1f}% {pur*100:7.1f}% {hic*100:6.1f}%")
    for k in ("structural_truth", "model", "levenshtein"):
        sc = res[k]
        if "pct_hi_clusters_sequence_diverse" in sc:
            print(f"  {sc['label']:24s} sequence-diverse share of high-purity "
                  f"clusters: {sc['pct_hi_clusters_sequence_diverse']*100:.1f}%")


    a = res["agreement_with_structure"]
    print("\n--- 3. agreement with the structural clustering (ARI) ---")
    print(f"  distilled model     {a['model_ari']:.3f}")
    print(f"  Levenshtein CDR3b   {a['levenshtein_ari']:.3f}")

    print("\n--- read this before quoting any of it ---")
    gain = p["model_auroc"] - p["levenshtein_auroc"]
    if np.isnan(gain):
        print("  AUROC undefined — too few positives or negatives.")
    elif gain <= 0.01:
        print("  The model does NOT beat edit distance. On this data the")
        print("  structural signal is not recoverable from sequence, or the")
        print("  cohort is too small to learn it. This is")
        print("  a real negative result about the limits of distillation.")
    else:
        print(f"  The model leads edit distance by {gain:+.3f} AUROC. Check the")
        print("  sequence-divergent block above before claiming it has learned")
        print("  structure rather than a better edit distance.")
        if not res["model"].get("reached_truth_purity", True):
            print()
            print("  BUT: a good pairwise ranker is not the same as a good")
            print("  clustering. The model never matches the structural purity at")
            print("  any threshold here, because connected components chain")
            print("  through borderline pairs. If you take this further, replace")
            print("  the threshold graph with something that resists chaining")
            print("  (average linkage, or a density method) before reporting")
            print("  cluster-level numbers.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assignments", help="cluster_assignments.tsv from Stage C")
    ap.add_argument("--labels", help="the labelled SPACE-TCR input TSV")
    ap.add_argument("--synthetic", action="store_true",
                    help="run on synthetic data; no inputs needed")
    ap.add_argument("--out", default="ml_out")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--emb-dim", type=int, default=64)
    ap.add_argument("--test-frac", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    if a.synthetic:
        print("  [synthetic] cohort whose signal sits in CDR3-alpha and V/J,")
        print("              where a CDR3-beta edit distance is blind to it")
        cohort = synthetic_cohort(seed=a.seed)
    else:
        if not (a.assignments and a.labels):
            ap.error("--assignments and --labels are required unless --synthetic")
        cohort = load_cohort(a.assignments, a.labels)
    print(f"  cohort: {cohort.n} receptors, "
          f"{cohort.df['cluster_id'].nunique()} clusters, "
          f"{cohort.df['epitope'].nunique()} epitopes")

    tr, te, info = epitope_split(cohort.df, a.test_frac, a.seed)
    if len(te) < 20 or te["cluster_id"].nunique() < 3:
        print("  test split too small to evaluate; lower --test-frac or supply "
              "more epitopes")
        return 1
    print(f"  training on {len(tr)} receptors ...")
    model, hist = train(cohort, tr, epochs=a.epochs, emb_dim=a.emb_dim,
                        seed=a.seed, quiet=a.quiet)
    res = evaluate(cohort, model, tr, te, seed=a.seed, quiet=a.quiet)
    report(res, info)

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "distill_results.json"), "w") as fh:
        json.dump({"split": info, "results": res, "history": hist,
                   "config": vars(a)}, fh, indent=2, default=float)
    pd.DataFrame(hist).to_csv(os.path.join(a.out, "distill_history.csv"),
                              index=False)
    print(f"\n  wrote {a.out}/distill_results.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
