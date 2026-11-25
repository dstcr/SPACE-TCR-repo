#!/usr/bin/env python3
"""Is structural cluster purity explained by germline V usage or CDR3 length?

    python3 scripts/cluster_controls.py [n_permutations]

Two alternative explanations of the structural clusters:

  1. V GENE. SPACE-TCR clusters are strongly V-homogeneous although the method
     receives no germline annotation. Receptors sharing a V gene share
     framework and germline-encoded CDR1/CDR2 geometry, so epitope purity could
     reflect the known association between V usage and specificity.

  2. CDR3 LENGTH. Models are binned on CDR lengths before any distance is
     computed, so clusters are length-pure (1 of 2,797 spans two CDR3 length
     pairs). Length correlates with specificity, so purity could come from the
     binning rather than from the geometry within a bin.

THE TEST

The clusters are held fixed and only the epitope labels are permuted, within
strata:

  global          no strata (chance)
  V gene          strata = (V-alpha, V-beta)
  CDR3 length     strata = (len CDR3-alpha, len CDR3-beta)
  V gene + length strata = both

Each null gives the purity the clusters would have if the geometry carried
nothing beyond that stratum; the gap between observed and null is the purity
attributable to geometry beyond it. Because the bins are defined by length,
the length null permutes within nearly the same partition the clustering
starts from, and so tests only the within-bin geometry.

Strata of size one cannot be permuted, so those receptors keep their true
label, which raises each null towards the observed value. Every null is
therefore conservative (V+length most of all), and the script reports how many
receptors are in single-member strata.

SCOPE

data/cohort_structural/cluster_assignments.tsv holds the RMSD assignments at
1.25 A (1,555 receptors in non-singleton clusters). Purity there is well below
the focal 0.75 A figure and is not comparable with Table 1.
"""
import csv, collections, random, statistics, sys, os

N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 200
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_SUM = os.path.join(HERE, "results/tables/extensions/cluster_controls.csv")
OUT_PERM = os.path.join(HERE, "results/tables/extensions/cluster_controls_permutations.csv")

ca = list(csv.DictReader(
    open(os.path.join(HERE, "data/cohort_structural/cluster_assignments.tsv")),
    delimiter="\t"))
meta = {r["id"]: r for r in csv.DictReader(
    open(os.path.join(HERE, "data/jouannet_labelled_input.tsv")), delimiter="\t")}

gene = lambda g: g.split("*")[0]
rows = []
for r in ca:
    m = meta.get(r["tcr_id"])
    if m:
        rows.append(dict(
            cluster=(r["bin_key"], r["cluster_id"]),
            epitope=m["epitope"],
            vpair=(gene(m["v_call_alpha"]), gene(m["v_call_beta"])),
            lpair=(len(m["cdr3_aa_alpha"]), len(m["cdr3_aa_beta"]))))
assert rows, "no receptors joined to labels"

byc = collections.defaultdict(list)
for i, r in enumerate(rows):
    byc[r["cluster"]].append(i)
clusters = [m for m in byc.values() if len(m) > 1]
y = [r["epitope"] for r in rows]


def purity(labels):
    """Eq. S1 over non-singleton clusters."""
    dom = tot = 0
    for m in clusters:
        cnt = collections.Counter(labels[i] for i in m)
        dom += cnt.most_common(1)[0][1]
        tot += len(m)
    return dom / tot


NULLS = {
    "Chance (labels shuffled freely)": lambda r: 0,
    "Matched on V gene": lambda r: r["vpair"],
    "Matched on CDR3 length": lambda r: r["lpair"],
    "Matched on V gene + CDR3 length": lambda r: (r["vpair"], r["lpair"]),
}

rng = random.Random(0)
obs = purity(y)
summary, perms = [], []
print(f"receptors joined         {len(rows)}")
print(f"non-singleton clusters   {len(clusters)} "
      f"({sum(len(m) for m in clusters)} receptors)")
print(f"observed purity          {obs:.4f}\n")
print(f"{'null':34} {'mean':>7} {'sd':>7} {'gain':>8} {'z':>7} {'stuck':>7}")

chance_mean = None
for name, keyf in NULLS.items():
    strata = collections.defaultdict(list)
    for i, r in enumerate(rows):
        strata[keyf(r)].append(i)
    stuck = sum(len(s) for s in strata.values() if len(s) == 1)

    vals = []
    for _ in range(N_PERM):
        p = [None] * len(y)
        for _, idx in strata.items():
            labs = [y[i] for i in idx]
            rng.shuffle(labs)
            for i, l in zip(idx, labs):
                p[i] = l
        vals.append(purity(p))
    m, sd = statistics.mean(vals), statistics.stdev(vals)
    if chance_mean is None:
        chance_mean = m
    z = (obs - m) / sd if sd > 0 else float("inf")
    share = (m - chance_mean) / (obs - chance_mean) if obs != chance_mean else 0.0
    print(f"{name:34} {m:7.4f} {sd:7.4f} {obs-m:+8.4f} {z:7.1f} {stuck:7d}")
    summary.append(dict(null=name, n_strata=len(strata),
                        n_receptors_in_strata_of_one=stuck,
                        null_mean=round(m, 6), null_sd=round(sd, 6),
                        purity_observed=round(obs, 6),
                        gain_over_null=round(obs - m, 6),
                        z_vs_null=round(z, 3),
                        permutations_ge_observed=sum(1 for v in vals if v >= obs),
                        fraction_of_gap_explained=round(share, 6)))
    perms += [dict(null=name, replicate=k, purity=round(v, 6))
              for k, v in enumerate(vals)]

os.makedirs(os.path.dirname(OUT_SUM), exist_ok=True)
with open(OUT_SUM, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0]))
    w.writeheader(); w.writerows(summary)
with open(OUT_PERM, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["null", "replicate", "purity"])
    w.writeheader(); w.writerows(perms)
print(f"\nwrote {OUT_SUM}\nwrote {OUT_PERM}")
