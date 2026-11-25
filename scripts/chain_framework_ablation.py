#!/usr/bin/env python3
"""What do the second chain and the framework alignment contribute?

    PYTHONPATH=. python3 scripts/chain_framework_ablation.py \
        --pdb-dir <folder of modelled PDBs, from Zenodo record 23765886> \
        --labels  data/jouannet_labelled_input.tsv \
        --out     results/tables/extensions/chain_framework_ablation.csv

Two choices define the method:

  * PAIRING. SPACE-TCR superposes and measures both chains' CDR loops. The
    sequence comparators are mostly CDR3-beta methods, so clustering on the
    beta loops alone is the structural analogue of what they do; the gap
    between that and the paired setting is what the second chain adds.

  * FRAMEWORK ALIGNMENT. Distances are Ca RMSD over the CDR loops after
    superposing on framework anchors, so they reflect where each loop sits
    relative to the framework as well as its own shape. Superposing on the
    loops themselves removes the relative placement and leaves internal
    conformation only; the gap between the two is what the framework alignment
    contributes. (Raw coordinate RMSD with no superposition is not a useful
    control: it depends on where the predictor placed each model.)

ARMS

  paired_fw      CDR_all   anchored on fwB+fwA     <- the published setting
  beta_fw        CDR_B_all anchored on fwB+fwA
  alpha_fw       CDR_A_all anchored on fwB+fwA
  beta_fwB       CDR_B_all anchored on fwB only    <- a beta-chain-only method
  alpha_fwA      CDR_A_all anchored on fwA only
  paired_noframe CDR_all   anchored on CDR_all     <- internal shape only
  beta_noframe   CDR_B_all anchored on CDR_B_all

Each arm is binned on the CDR lengths it uses, so a beta-only arm is not split
by alpha lengths it does not measure. Clustering is single linkage (connected
components of the thresholded graph), as in the pipeline.

COST: all pairs within bins, one RMSD per pair per arm. The six-loop binning
gives ~40k pairs for the 4,031-receptor cohort; single-chain arms bin on three
loops, so their bins are larger and there are many more pairs. Use --max-bin to
cap bin size, or --epitope to restrict to one specificity.
"""
import argparse, collections, csv, json, os, sys
import numpy as np
from joblib import Parallel, delayed

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from SPACE2TCR.util import reg_def, rmsd, parse_antibodies  # noqa: E402

ARMS = {
    "paired_fw":      ("CDR_all",   "fw_all"),
    "beta_fw":        ("CDR_B_all", "fw_all"),
    "alpha_fw":       ("CDR_A_all", "fw_all"),
    "beta_fwB":       ("CDR_B_all", "fwB_only"),
    "alpha_fwA":      ("CDR_A_all", "fwA_only"),
    "paired_noframe": ("CDR_all",   "CDR_all"),
    "beta_noframe":   ("CDR_B_all", "CDR_B_all"),
}
# Which CDR regions define the length bins for each selection. A beta-only arm
# bins on the beta loops only.
BIN_REGIONS = {"CDR_all":   ["CDRB1", "CDRB2", "CDRB3", "CDRA1", "CDRA2", "CDRA3"],
               "CDR_B_all": ["CDRB1", "CDRB2", "CDRB3"],
               "CDR_A_all": ["CDRA1", "CDRA2", "CDRA3"]}


def regions(name):
    if name == "fwB_only":
        return np.asarray(reg_def["fwB"])
    if name == "fwA_only":
        return np.asarray(reg_def["fwA"])
    return np.concatenate(reg_def[name])


def loop_lengths(numbers, region_names):
    """Observed length of each named CDR region for one model."""
    present = set(int(x) for x in numbers)
    return tuple(sum(1 for r in reg_def[rn] if int(r) in present) for rn in region_names)


def components(n, edges):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    return [find(i) for i in range(n)]


def purity_retention(labels_by_model, groups, n_total):
    """Purity over LABELLED members; retention over ALL receptors.

    On the noise ladder most receptors carry no epitope annotation. They take
    part in the clustering, since contaminating clusters is what background
    does, but they cannot contribute to purity, which is agreement between
    labels. An unlabelled member (label None) counts towards the cluster's size
    for retention and is left out of the purity numerator and denominator.

    A cluster with fewer than two members, or fewer than two LABELLED members,
    contributes nothing to purity; it still counts for retention if it is
    non-singleton overall.

    The two high-purity fractions are those Jouannet et al. report alongside
    pooled purity, and that results/tables/published_method_metrics.csv gives
    for the nine published methods: the share of clusters whose own purity is
    >= 0.9, and the share of clustered receptors in those clusters. A cluster's
    own purity here is dominant / labelled, the same denominator as the pooled
    purity above. On the curated cohort every clustered receptor is labelled,
    so this matches the pipeline (which divides by full cluster size); on the
    noise ladder the two differ, so ladder rows are comparable across arms but
    not with the published methods, which were not scored on a noisy cohort."""
    HI = 0.9
    dom = tot = nclus = 0
    retained = 0
    mono = 0            # clusters carrying exactly one epitope
    n_epi = []          # distinct epitopes per scored cluster
    n_hi = 0            # scored clusters whose own purity clears HI
    sz_scored = sz_hi = 0
    for mem in groups:
        if len(mem) < 2: continue
        retained += len(mem)
        lab = [labels_by_model[m] for m in mem if labels_by_model[m] is not None]
        if len(lab) < 2: continue
        cnt = collections.Counter(lab)
        top = cnt.most_common(1)[0][1]
        dom += top; tot += len(lab); nclus += 1
        n_epi.append(len(cnt))
        if len(cnt) == 1: mono += 1
        sz_scored += len(mem)
        if top / len(lab) >= HI:
            n_hi += 1; sz_hi += len(mem)
    poly = 1 - mono / nclus if nclus else float("nan")
    mean_epi = (sum(n_epi) / len(n_epi)) if n_epi else float("nan")
    pc_hi = n_hi / nclus if nclus else float("nan")
    pt_hi = sz_hi / sz_scored if sz_scored else float("nan")
    return ((dom / tot if tot else float("nan")), retained / n_total, nclus,
            poly, mean_epi, pc_hi, pt_hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--thresholds", default="0.5,0.75,1.0,1.25")
    ap.add_argument("--max-bin", type=int, default=100000,
                    help="skip bins larger than this (pair count is quadratic). "
                         "The default is effectively no cap: a low cap silently "
                         "drops most of a single-chain arm, because binning on "
                         "three loops instead of six makes bins far larger, and "
                         "the arms then stop being comparable.")
    ap.add_argument("--arms", default=None,
                    help="comma-separated subset of arms to run (default: all)")
    ap.add_argument("--bin-on", default=None, choices=list(BIN_REGIONS),
                    help="force every arm to use ONE binning instead of binning "
                         "on the loops it measures. The default (per-arm "
                         "binning) is what a single-chain method would do, "
                         "but it mixes binning into the comparison: binning on three loops "
                         "gives 26 bins for beta and 68 for alpha against 477 "
                         "for the paired set, and a coarser pre-filter puts "
                         "more candidates in each bin, which depresses purity "
                         "for reasons unrelated to the distance. Pass "
                         "CDR_all to hold the binning fixed and vary only the "
                         "distance measure.")
    ap.add_argument("--universe", type=int, default=4292,
                    help="denominator for the cross-method unique-CDR3b unit "
                         "(default 4,292, the paper's universe)")
    ap.add_argument("--subset-ids", default=None,
                    help="file whose first tab-separated column lists the "
                         "receptor ids to include (a noise-ladder rung)")
    ap.add_argument("--epitope", default=None,
                    help="restrict to one epitope, for a fast first look")
    ap.add_argument("--n-jobs", type=int, default=-1)
    a = ap.parse_args()

    # An epitope of "NA" or "" means unlabelled, not absent: those receptors
    # join the clustering and are skipped when purity is computed.
    lab, cdr3b = {}, {}
    for r in csv.DictReader(open(a.labels), delimiter="\t"):
        e = (r.get("epitope") or "").strip()
        lab[r["id"]] = None if e in ("", "NA", "NaN", "None") else e
        # Where the label table carries the beta sequence, the arms are also
        # reported on the paper's cross-method unit: unique CDR3-beta over the
        # 4,292-sequence universe. Without it those columns are blank.
        b = (r.get("cdr3_aa_beta") or r.get("cdr3b") or "").strip()
        if b: cdr3b[r["id"]] = b
    n_unlab = sum(1 for v in lab.values() if v is None)
    if n_unlab:
        print(f"labels: {len(lab)} receptors, {n_unlab} unlabelled "
              f"({100*n_unlab/len(lab):.1f}% noise)", file=sys.stderr)
    if a.subset_ids:
        keep = {ln.split("\t")[0].strip()
                for ln in open(a.subset_ids) if ln.strip()}
        lab = {k: v for k, v in lab.items() if k in keep}
        print(f"restricted to {len(lab)} receptors from {a.subset_ids}",
              file=sys.stderr)
    if a.epitope:
        lab = {k: v for k, v in lab.items() if v == a.epitope}

    files = sorted(os.path.join(a.pdb_dir, f) for f in os.listdir(a.pdb_dir)
                   if f.endswith(".pdb") and f[:-4] in lab)
    if not files:
        sys.exit(f"no .pdb files in {a.pdb_dir} matching the label table.\n"
                 f"If the structures are not there, this is the script that "
                 f"needs them -- see README. Refusing to invent a result.")
    print(f"models: {len(files)}", file=sys.stderr)

    models = parse_antibodies(files, n_jobs=a.n_jobs)
    ids = [os.path.basename(f)[:-4] for f in files]
    ylab = [lab[i] for i in ids]
    n_lab = sum(1 for v in ylab if v is not None)
    print(f"clustering {len(ids)} models, {n_lab} of them labelled",
          file=sys.stderr)
    thresholds = [float(t) for t in a.thresholds.split(",")]

    arms = ARMS
    if a.arms:
        want = [x.strip() for x in a.arms.split(",")]
        bad = [w for w in want if w not in ARMS]
        if bad:
            sys.exit(f"unknown arm(s) {bad}; choices are {list(ARMS)}")
        arms = {k: ARMS[k] for k in want}

    rows = []
    for arm, (sel_name, anc_name) in arms.items():
        sel = regions(sel_name); anc = regions(anc_name)
        bin_name = a.bin_on or sel_name
        bins = collections.defaultdict(list)
        for k, m in enumerate(models):
            bins[loop_lengths(m[0], BIN_REGIONS[bin_name])].append(k)
        skipped = sum(len(v) for v in bins.values() if len(v) > a.max_bin)
        pairs = sum(len(v) * (len(v) - 1) // 2
                    for v in bins.values() if 1 < len(v) <= a.max_bin)
        print(f"{arm:15} bins={len(bins):5} pairs={pairs:,} "
              f"skipped={skipped}", file=sys.stderr)

        def uniq_beta(groups):
            """Unique CDR3-beta sequences sitting in non-singleton clusters."""
            if not cdr3b: return None
            seen = set()
            for mem in groups:
                if len(mem) < 2: continue
                for m in mem:
                    b = cdr3b.get(ids[m])
                    if b: seen.add(b)
            return len(seen)

        def one_bin(idxs):
            """All-pairs within one bin, then components at each threshold."""
            d = {}
            for ii in range(len(idxs)):
                mi = models[idxs[ii]]
                for jj in range(ii + 1, len(idxs)):
                    try:
                        v = float(rmsd(mi, models[idxs[jj]],
                                       selection=sel, anchors=anc))
                    except Exception:
                        continue
                    if np.isfinite(v): d[(ii, jj)] = v
            out = {}
            for t in thresholds:
                lab_i = components(len(idxs), [p for p, v in d.items() if v <= t])
                g = collections.defaultdict(list)
                for pos, root in enumerate(lab_i): g[root].append(idxs[pos])
                out[t] = list(g.values())
            return out

        todo = [v for v in bins.values() if 1 < len(v) <= a.max_bin]
        todo.sort(key=len, reverse=True)          # longest bins first
        results = Parallel(n_jobs=a.n_jobs, backend="loky", verbose=0)(
            delayed(one_bin)(idxs) for idxs in todo)
        per_thr = {t: [] for t in thresholds}
        for res in results:
            for t in thresholds: per_thr[t].extend(res[t])

        for t in thresholds:
            p, r, nc, poly, mepi, pc_hi, pt_hi = purity_retention(
                {k: ylab[k] for k in range(len(ids))}, per_thr[t], len(ids))
            rows.append(dict(arm=arm, selection=sel_name, anchors=anc_name,
                             binned_on=bin_name,
                             threshold=t, purity=round(p, 6),
                             retention=round(r, 6), n_clusters=nc,
                             perc_clusters_gt09=round(pc_hi, 6),
                             perc_tcrs_gt09=round(pt_hi, 6),
                             polyspecific=round(poly, 6),
                             mean_epitopes=round(mepi, 6),
                             n_models=len(ids),
                             n_unique_cdr3b=uniq_beta(per_thr[t]),
                             retention_common=(
                                 round(uniq_beta(per_thr[t]) / a.universe, 6)
                                 if cdr3b else None),
                             universe_n=a.universe if cdr3b else None,
                             epitope=a.epitope or "ALL"))
            print(f"  {arm:15} thr {t:4.2f}  purity {p:.4f}  "
                  f"retention {r:.4f}  clusters {nc}  "
                  f"polyspec {poly:.4f}  %clus>=.9 {pc_hi:.4f}  "
                  f"%tcrs>=.9 {pt_hi:.4f}", file=sys.stderr)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
