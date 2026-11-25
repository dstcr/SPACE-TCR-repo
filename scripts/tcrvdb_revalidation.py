#!/usr/bin/env python3
"""Re-benchmark SPACE-TCR against functionally validated epitope labels.

    python3 scripts/tcrvdb_revalidation.py /path/to/TCRvdb.csv

BACKGROUND

The purity, sensitivity and polyspecificity numbers in the paper are measured
against epitope labels from IEDB / VDJdb / McPAS, through the Jouannet cohort.
Messemaker et al. (Immunity 59:2775-2788, 2026; doi
10.1016/j.immuni.2026.09.007) synthesised 3,693 VDJdb TCRs and screened them
functionally; about 56% reproduced their annotated reactivity.

Purity (Eq. S1) measures agreement with the label, not with reactivity: ten
receptors all annotated NLV score purity 1.0 whether or not any is NLV-reactive.
The direction of any bias depends on how the label errors are structured.
Random errors lower measured purity (co-specific receptors would carry mixed
labels); errors correlated within a source study raise it (receptors reported
together often come from one donor and one expansion, so they are similar in
sequence and share the same error). This script measures which applies.

SCHEMA (release 2026-10-01, ttd/gil/ylq/nlv/glc)

Reactivity is three binary columns, `pos`, `neg` and `ambiguous`. Keys used:
`cdr3_alpha_aa`, `cdr3_beta_aa`, `epitope_aa`. The script stops if these are
missing.

TWO ANALYSES, BECAUSE OF THE JOIN KEY

  A. FOCAL 0.75 A, unit = unique CDR3-beta. This is the paper's headline
     setting, but results/tables/extensions/cluster_overview_membership.csv
     records membership by CDR3-beta only, so the join uses beta and epitope.
     5.0% of beta+epitope keys in TCRvdb carry conflicting calls (the same
     beta validating with one alpha and not another); these are treated as
     ambiguous rather than resolved by majority.

  B. 1.25 A, unit = paired receptor. data/cohort_structural/
     cluster_assignments.tsv is keyed on tcr_id, so the join can use the full
     alpha+beta+epitope key, where only 3 of 2,029 keys conflict. The threshold
     is looser but the join is exact. If A and B disagree, the join is the
     likelier cause.

TCRvdb is available on request (tcrvdb.nkiai.nl) and its terms do not allow
redistribution, so it is not included here. Pass the path to your copy.
"""
import csv, collections, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COHORT = os.path.join(HERE, "data/jouannet_labelled_input.tsv")
MEMBER = os.path.join(HERE, "results/tables/extensions/cluster_overview_membership.csv")
ASSIGN = os.path.join(HERE, "data/cohort_structural/cluster_assignments.tsv")
OUTDIR = os.path.join(HERE, "results/tables/extensions")

TCRVDB_EPITOPES = ["GILGFVFTL", "NLVPMVATV", "GLCTLVAML", "YLQPRTFLL",
                   "TTDPSFLGRY"]
FOCAL = "SPACE-TCR RMSD (0.75A)"
WIDE = "SPACE-TCR RMSD (1A)"
COMPARATOR = "clusTCR"
REQUIRED = ["pos", "neg", "ambiguous", "cdr3_alpha_aa", "cdr3_beta_aa",
            "epitope_aa"]


def call_of(row):
    if row["pos"] == "1":
        return "validating"
    if row["neg"] == "1":
        return "non-validating"
    return "ambiguous"


def collapse(calls):
    """One call per key. A key seen both validating and non-validating is
    ambiguous -- not a majority vote. Reactivity that depends on the partner
    chain is a real phenomenon, not noise to be averaged away."""
    s = set(calls)
    if s == {"validating"}:
        return "validating"
    if s == {"non-validating"}:
        return "non-validating"
    return "ambiguous"


def purity(groups, label_of, drop=frozenset()):
    """Eq. S1 over non-singleton clusters, after dropping excluded keys."""
    dom = tot = nclus = 0
    for mem in groups.values():
        mem = [k for k in mem if k not in drop and k in label_of]
        if len(mem) < 2:
            continue
        cnt = collections.Counter(label_of[k] for k in mem)
        dom += cnt.most_common(1)[0][1]
        tot += len(mem)
        nclus += 1
    return (dom / tot if tot else float("nan")), tot, nclus


def fisher(a11, a12, a21, a22):
    from scipy.stats import fisher_exact
    return fisher_exact([[a11, a12], [a21, a22]], alternative="greater")


def enrichment(label_name, keys_clustered, keys_all, call):
    a11 = a12 = a21 = a22 = 0
    for k in keys_all:
        c = call.get(k)
        if c not in ("validating", "non-validating"):
            continue
        inside = k in keys_clustered
        if inside and c == "validating":   a11 += 1
        elif inside:                       a12 += 1
        elif c == "validating":            a21 += 1
        else:                              a22 += 1
    if min(a11 + a12, a21 + a22) == 0:
        print(f"    {label_name:26} too few receptors to test")
        return None
    orr, p = fisher(a11, a12, a21, a22)
    r_in, r_out = a11 / (a11 + a12), a21 / (a21 + a22)
    print(f"    {label_name:26} clustered {r_in:.3f} ({a11}/{a11+a12})   "
          f"unclustered {r_out:.3f} ({a21}/{a21+a22})   "
          f"OR {orr:.2f}   p {p:.3g}")
    return dict(method=label_name, clustered_validating=a11,
                clustered_nonvalidating=a12, unclustered_validating=a21,
                unclustered_nonvalidating=a22, rate_clustered=round(r_in, 4),
                rate_unclustered=round(r_out, 4), odds_ratio=round(orr, 4),
                p_one_sided=p)


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: tcrvdb_revalidation.py /path/to/TCRvdb.csv")
    tv = list(csv.DictReader(open(sys.argv[1])))
    missing = [c for c in REQUIRED if c not in tv[0]]
    if missing:
        sys.exit(f"TCRvdb is missing expected columns {missing}.\n"
                 f"Columns present: {list(tv[0])}\nRefusing to guess.")
    print(f"TCRvdb rows {len(tv)}   "
          f"{collections.Counter(call_of(r) for r in tv)}")

    by_ab, by_b = collections.defaultdict(list), collections.defaultdict(list)
    for r in tv:
        c = call_of(r)
        by_ab[(r["cdr3_alpha_aa"], r["cdr3_beta_aa"], r["epitope_aa"])].append(c)
        by_b[(r["cdr3_beta_aa"], r["epitope_aa"])].append(c)
    call_ab = {k: collapse(v) for k, v in by_ab.items()}
    call_b = {k: collapse(v) for k, v in by_b.items()}
    conflicted = sum(1 for k, v in by_b.items()
                     if len(set(v)) > 1 and "ambiguous" not in set(v))
    print(f"beta+epitope keys {len(call_b)} ({conflicted} demoted to ambiguous "
          f"on partner-dependent reactivity);  alpha+beta+epitope keys "
          f"{len(call_ab)}")

    cohort = [r for r in csv.DictReader(open(COHORT), delimiter="\t")
              if r["epitope"] in TCRVDB_EPITOPES]
    print(f"\ncohort receptors in these 5 epitopes: {len(cohort)}")

    rows_out, enr_out, pur_out = [], [], []

    # ================= A. focal 0.75 A, unit = unique CDR3-beta =============
    print("\n" + "=" * 72)
    print("A. FOCAL 0.75 A   unit = unique CDR3-beta   join = beta + epitope")
    print("=" * 72)
    members = list(csv.DictReader(open(MEMBER)))
    clustered = collections.defaultdict(set)
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for m in members:
        k = (m["cdr3b"], m["epitope"])
        clustered[m["method"]].add(k)
        groups[m["method"]][m["cluster_id"]].append(k)

    keys_all = {(r["cdr3_aa_beta"], r["epitope"]) for r in cohort}
    label = {k: k[1] for k in keys_all}
    hit = {k for k in keys_all if k in call_b}
    print(f"unique CDR3b in cohort {len(keys_all)};  matched in TCRvdb "
          f"{len(hit)} ({100*len(hit)/len(keys_all):.1f}%)")

    print(f"\n{'epitope':12} {'cohort':>7} {'matched':>8} {'valid':>7} "
          f"{'nonvalid':>9} {'amb':>5} {'rate':>7}")
    for e in TCRVDB_EPITOPES + ["ALL"]:
        ks = keys_all if e == "ALL" else {k for k in keys_all if k[1] == e}
        h = [k for k in ks if k in call_b]
        v = sum(1 for k in h if call_b[k] == "validating")
        nv = sum(1 for k in h if call_b[k] == "non-validating")
        am = len(h) - v - nv
        rate = v / (v + nv) if v + nv else float("nan")
        print(f"{e:12} {len(ks):7} {len(h):8} {v:7} {nv:9} {am:5} {rate:7.3f}")
        rows_out.append(dict(analysis="focal_0.75A_cdr3b", epitope=e,
                             cohort_keys=len(ks), matched=len(h), validating=v,
                             nonvalidating=nv, ambiguous=am,
                             validation_rate=round(rate, 4)))

    nonval = {k for k in keys_all if call_b.get(k) == "non-validating"}
    ambig = {k for k in keys_all if call_b.get(k) == "ambiguous"}
    unmatched = keys_all - hit
    label_rel = dict(label)
    for k in nonval:
        label_rel[k] = "NON-REACTIVE"

    print("\npurity, restricted to these 5 epitopes:")
    for meth in (FOCAL, WIDE, COMPARATOR):
        g = groups[meth]
        p1 = purity(g, label)
        p2 = purity(g, label, drop=nonval | ambig | unmatched)
        p3 = purity(g, label_rel)
        print(f"  {meth}")
        for nm, (p, n, k) in [("as published", p1),
                              ("validated only", p2),
                              ("non-validating relabelled", p3)]:
            print(f"      {nm:28} {p:.4f}  ({n} CDR3b in {k} clusters)")
            pur_out.append(dict(analysis="focal_0.75A_cdr3b", method=meth,
                                basis=nm.replace(" ", "_"), purity=round(p, 6),
                                n_keys=n, n_clusters=k))

    print("\nare clustered CDR3b enriched for VALIDATING annotations?")
    for meth in (FOCAL, WIDE, COMPARATOR):
        e = enrichment(meth, clustered[meth], keys_all, call_b)
        if e:
            e["analysis"] = "focal_0.75A_cdr3b"
            enr_out.append(e)

    # ================= B. permissive 1.25 A, unit = paired receptor =========
    print("\n" + "=" * 72)
    print("B. PERMISSIVE 1.25 A   unit = paired receptor   "
          "join = alpha + beta + epitope")
    print("=" * 72)
    ca = list(csv.DictReader(open(ASSIGN), delimiter="\t"))
    cid = {r["tcr_id"]: (r["bin_key"], r["cluster_id"]) for r in ca}
    sz = collections.Counter(cid.values())

    keys_b, grp_b, clustered_b = set(), collections.defaultdict(list), set()
    label_b = {}
    for r in cohort:
        k = (r["cdr3_aa_alpha"], r["cdr3_aa_beta"], r["epitope"])
        keys_b.add(k)
        label_b[k] = r["epitope"]
        c = cid.get(r["id"])
        if c is not None and sz[c] > 1:
            clustered_b.add(k)
            grp_b[c].append(k)
    hit_b = {k for k in keys_b if k in call_ab}
    print(f"paired receptors {len(keys_b)};  matched in TCRvdb {len(hit_b)} "
          f"({100*len(hit_b)/len(keys_b):.1f}%);  in non-singleton clusters "
          f"{len(clustered_b)}")

    nv_b = {k for k in keys_b if call_ab.get(k) == "non-validating"}
    am_b = {k for k in keys_b if call_ab.get(k) == "ambiguous"}
    un_b = keys_b - hit_b
    lab_rel_b = dict(label_b)
    for k in nv_b:
        lab_rel_b[k] = "NON-REACTIVE"
    print("\npurity of the 1.25 A clusters:")
    for nm, kw in [("as published", {}),
                   ("validated only", dict(drop=nv_b | am_b | un_b)),
                   ("non-validating relabelled", None)]:
        p, n, k = (purity(grp_b, lab_rel_b) if kw is None
                   else purity(grp_b, label_b, **kw))
        print(f"    {nm:28} {p:.4f}  ({n} receptors in {k} clusters)")
        pur_out.append(dict(analysis="permissive_1.25A_paired",
                            method="SPACE-TCR RMSD (1.25A)",
                            basis=nm.replace(" ", "_"), purity=round(p, 6),
                            n_keys=n, n_clusters=k))

    print("\nare clustered receptors enriched for VALIDATING annotations?")
    e = enrichment("SPACE-TCR RMSD (1.25A)", clustered_b, keys_b, call_ab)
    if e:
        e["analysis"] = "permissive_1.25A_paired"
        enr_out.append(e)

    # ================= C. is the enrichment a confound? ====================
    # Two things could produce it without structural information:
    #   (i)  VDJdb's own confidence score, if clustering tracked it;
    #   (ii) sequence "publicness": receptors with many near-identical
    #        relatives are both easier to cluster and more likely to be real.
    # Both are tested by stratifying and repeating the test within each
    # stratum. ClusTCR goes through the same procedure for comparison.
    print("\n" + "=" * 72)
    print("C. CONTROLS  (stratified; pooled ORs are confounded by stratum)")
    print("=" * 72)
    strat_out = []

    vdj = collections.defaultdict(set)
    for r in tv:
        vdj[(r["cdr3_beta_aa"], r["epitope_aa"])].add(r["vdjdb_score"])
    vdj = {k: max(v) for k, v in vdj.items()}

    from rapidfuzz.distance import Levenshtein
    from rapidfuzz import process
    uniq = sorted({k[0] for k in keys_all})
    M = process.cdist(uniq, uniq, scorer=Levenshtein.distance, workers=-1)
    deg = dict(zip(uniq, (M <= 2).sum(1) - 1))

    def bucket(d):
        return "0" if d == 0 else ("1-2" if d <= 2 else ("3-9" if d <= 9 else "10+"))

    for cname, keyf, order in [
            ("VDJdb confidence", lambda k: vdj.get(k),
             ["0.0", "1.0", "2.0", "3.0"]),
            ("LD<=2 sequence degree", lambda k: bucket(deg[k[0]]),
             ["0", "1-2", "3-9", "10+"])]:
        print(f"\n  stratified on {cname}")
        for meth in (FOCAL, COMPARATOR):
            print(f"    {meth}")
            for s in order:
                a = b = c = d = 0
                for k in keys_all:
                    if keyf(k) != s:
                        continue
                    cl = call_b.get(k)
                    if cl not in ("validating", "non-validating"):
                        continue
                    ins = k in clustered[meth]
                    if ins and cl == "validating":  a += 1
                    elif ins:                       b += 1
                    elif cl == "validating":        c += 1
                    else:                           d += 1
                if min(a + b, c + d) == 0 or (a + b) < 5:
                    print(f"       {s:5} too few (clustered {a+b})")
                    continue
                orr, p = fisher(a, b, c, d)
                print(f"       {s:5} clustered {a/(a+b):.3f} ({a}/{a+b})  "
                      f"unclustered {c/(c+d):.3f} ({c}/{c+d})  "
                      f"OR {orr:.2f}  p {p:.3g}")
                strat_out.append(dict(control=cname, method=meth, stratum=s,
                                      clustered_validating=a,
                                      clustered_nonvalidating=b,
                                      unclustered_validating=c,
                                      unclustered_nonvalidating=d,
                                      rate_clustered=round(a/(a+b), 4),
                                      rate_unclustered=round(c/(c+d), 4),
                                      odds_ratio=round(orr, 4), p_one_sided=p))

    os.makedirs(OUTDIR, exist_ok=True)
    for name, data in [("tcrvdb_validation_by_epitope.csv", rows_out),
                       ("tcrvdb_purity_recomputed.csv", pur_out),
                       ("tcrvdb_cluster_enrichment.csv", enr_out),
                       ("tcrvdb_enrichment_stratified.csv", strat_out)]:
        if not data:
            continue
        with open(os.path.join(OUTDIR, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader(); w.writerows(data)
    print(f"\nwrote 4 tables to {OUTDIR}")


if __name__ == "__main__":
    main()
