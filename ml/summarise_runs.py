#!/usr/bin/env python3
"""Summarise a set of distill.py runs, with a paired test across seeds.

    python3 ml/summarise_runs.py ml/results/*.json
    python3 ml/summarise_runs.py ml/results/{125es,seed1,seed2,seed3,seed4}.json

A single run of distill.py on a cohort this small tells you very little: across
five seeds on the same 1.25 A supervision the overall AUROC ranged 0.684 to
0.887, purely from which 53 epitopes landed in the held-out set. This reduces a
set of runs to the question that matters -- does the model beat edit distance,
and by how much, with what spread -- and runs a paired t-test over the seeds.

Pass only runs that ARE replicates of each other. Runs at different supervision
thresholds are not: they are nested clusterings of one cohort, so pooling them
as if they were independent overstates the evidence.

No dependencies beyond the standard library, so it runs wherever distill.py did.
"""
import json
import math
import statistics as st
import sys
from pathlib import Path


def _betacf(a, b, x, maxit=200, eps=3e-14, fpmin=1e-300):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = fpmin if abs(d) < fpmin else d
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d; d = fpmin if abs(d) < fpmin else d
        c = 1.0 + aa / c; c = fpmin if abs(c) < fpmin else c
        d = 1.0 / d; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d; d = fpmin if abs(d) < fpmin else d
        c = 1.0 + aa / c; c = fpmin if abs(c) < fpmin else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def _betai(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lb = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
          + a * math.log(x) + b * math.log1p(-x))
    bt = math.exp(lb)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def paired_t(diffs):
    """Two-sided paired t-test against zero. Returns (mean, sd, t, p, lo, hi)."""
    n = len(diffs)
    if n < 2:
        return (diffs[0] if n else float("nan"),) + (float("nan"),) * 5
    m, s = st.mean(diffs), st.stdev(diffs)
    if s == 0:
        return m, 0.0, float("inf"), 0.0, m, m
    se = s / math.sqrt(n)
    t = m / se
    df = n - 1
    p = _betai(df / 2.0, 0.5, df / (df + t * t))
    # two-sided 95% critical values for small df; falls back to the normal
    crit = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
            7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228}.get(df, 1.96)
    return m, s, t, p, m - crit * se, m + crit * se


def load(path):
    d = json.loads(Path(path).read_text())
    r, s = d["results"], d["split"]
    hist = [h for h in d.get("history", []) if h["val_auroc"] == h["val_auroc"]]
    best_ep = max(hist, key=lambda h: h["val_auroc"])["epoch"] if hist else None
    return {
        "name": Path(path).stem,
        "seed": d.get("config", {}).get("seed"),
        "best_epoch": best_ep,
        "n_epochs": len(d.get("history", [])),
        "n_train": s["n_train"], "n_test": s["n_test"],
        "n_epitopes_test": s["n_epitopes_test"],
        "model": r["pair_level"]["model_auroc"],
        "ld": r["pair_level"]["levenshtein_auroc"],
        "model_div": r["pair_level_sequence_diverse_only"]["model_auroc"],
        "ld_div": r["pair_level_sequence_diverse_only"]["levenshtein_auroc"],
        "model_ari": r["agreement_with_structure"]["model_ari"],
        "ld_ari": r["agreement_with_structure"]["levenshtein_ari"],
        "reached_purity": r["model"].get("reached_truth_purity"),
    }


def main(paths):
    runs = [load(p) for p in paths]
    if not runs:
        raise SystemExit("no runs given")

    print()
    print(f"{'run':<12}{'seed':>5}{'bestEp':>7}{'/ep':>5}{'test':>6}{'epi':>5} | "
          f"{'model':>6}{'LD':>7}{'diff':>7} | {'mDIV':>6}{'ldDIV':>7}{'diff':>7} | "
          f"{'mARI':>6}{'lARI':>6}  pur?")
    print("-" * 104)
    for r in runs:
        print(f"{r['name']:<12}{str(r['seed']):>5}{str(r['best_epoch']):>7}"
              f"{r['n_epochs']:>5}{r['n_test']:>6}{r['n_epitopes_test']:>5} | "
              f"{r['model']:>6.3f}{r['ld']:>7.3f}{r['model']-r['ld']:>+7.3f} | "
              f"{r['model_div']:>6.3f}{r['ld_div']:>7.3f}"
              f"{r['model_div']-r['ld_div']:>+7.3f} | "
              f"{r['model_ari']:>6.3f}{r['ld_ari']:>6.3f}  {r['reached_purity']}")

    if len(runs) < 2:
        return
    print()
    print("paired across the runs given (only valid if they are replicates):")
    print()
    for key, label in (("div", "divergent-pair advantage"),
                       ("all", "overall advantage")):
        if key == "div":
            d = [r["model_div"] - r["ld_div"] for r in runs]
        else:
            d = [r["model"] - r["ld"] for r in runs]
        m, s, t, p, lo, hi = paired_t(d)
        wins = sum(1 for x in d if x > 0)
        sig = "" if p < 0.05 else "   NOT SIGNIFICANT"
        print(f"  {label:<26} {m:+.3f} +/- {s:.3f}   "
              f"t({len(d)-1})={t:.2f}  p={p:.3f}  "
              f"95% CI [{lo:+.3f}, {hi:+.3f}]  wins {wins}/{len(d)}{sig}")

    print()
    if not any(r["reached_purity"] for r in runs):
        print("  No run reached the structural clustering's purity at any threshold.")
        print("  Cluster-level numbers are not usable until the threshold graph is")
        print("  replaced with something chaining-resistant.")
    eps = [r["best_epoch"] for r in runs if r["best_epoch"]]
    if eps and max(eps) <= 3:
        print(f"  Best epoch was <= {max(eps)} in every run. An encoder whose best")
        print("  weights are that young has not learned a metric; fix the selection")
        print("  (bigger inner-validation set, patience schedule) before judging it.")
    print()


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        print(__doc__)
        raise SystemExit(0)
    main(args)
