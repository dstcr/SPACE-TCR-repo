#!/usr/bin/env python3
"""Does the beta chain's predicted shape depend on which alpha it was modelled with?

    python3 scripts/dummy_alpha_control.py seqs   --pdb-dir DIR --out seqs.tsv
    python3 scripts/dummy_alpha_control.py verify --seqs seqs.tsv --pdb-dir DIR -n 10
    python3 scripts/dummy_alpha_control.py bins   --pdb-dir DIR --seqs seqs.tsv \
            --target 400 --out-ids ids.txt
    python3 scripts/dummy_alpha_control.py model  --seqs seqs.tsv --ids ids.txt \
            --mode dummy --out-dir DIR

TCRBuilder2 predicts the alpha-beta Fv jointly, so every beta coordinate in the
cohort was generated with that receptor's own alpha chain present. The "beta
only" arm of the chain ablation is therefore a beta-shaped readout of a paired
model, and some apparent beta-chain signal could be alpha information passing
through the predictor. This control measures that.

THE CONTROL

Every receptor is re-modelled as {A: one fixed alpha, B: its own beta}. The
alpha partner is then constant, so differences between receptors' beta
coordinates come from the beta sequence alone. Those betas are clustered and
compared with the betas from the real paired models:

  curves coincide  ->  the effect is negligible; the beta arm measures beta
                       geometry
  dummy is worse   ->  the gap measures the alpha contribution

TCRBuilder2 cannot model one chain alone (predict({'B': seq}) raises
KeyError 'A'), so a constant partner is the closest approximation: the
predictor still models an interface, but the same one every time.

WHOLE BINS

Clustering compares receptors only within a CDR-length bin, so taking whole
bins reproduces exactly the clusters the full run gives for those bins, which a
random subsample would not. At ~19 s per structure on two cores the full 4,031
receptors are impractical, so the control uses whole bins.

ARM

beta_fwB: beta CDR loops anchored on the beta framework only. It uses no alpha
coordinates, so the dummy alpha can affect it only through the predictor.
"""
import argparse, collections, csv, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
       "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
       "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
       "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def chain_seqs(path):
    """One-letter sequence per chain, in residue order, from CA records."""
    out, seen = collections.defaultdict(list), collections.defaultdict(set)
    for l in open(path):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            ch, resi = l[21], l[22:27]
            if resi in seen[ch]:
                continue
            seen[ch].add(resi)
            out[ch].append(AA3.get(l[17:20].strip(), "X"))
    return {c: "".join(v) for c, v in out.items()}


def cmd_seqs(a):
    files = sorted(f for f in os.listdir(a.pdb_dir) if f.endswith(".pdb"))
    n_bad = 0
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "alpha", "beta"])
        for f in files:
            s = chain_seqs(os.path.join(a.pdb_dir, f))
            if "A" not in s or "B" not in s or "X" in s["A"] + s["B"]:
                n_bad += 1
                continue
            w.writerow([f[:-4], s["A"], s["B"]])
    print(f"wrote {a.out}  ({len(files) - n_bad} of {len(files)} models; "
          f"{n_bad} skipped for a missing chain or an unknown residue)",
          file=sys.stderr)


def _builder():
    from ImmuneBuilder import TCRBuilder2
    return TCRBuilder2()


def cmd_verify(a):
    """Re-model N receptors with their REAL alpha and compare with the PDB on
    disk. This separates two explanations for any later difference: a genuine
    alpha effect, or simply a different predictor version from the one that
    built the cohort. If these come back near zero the cohort PDBs and today's
    ImmuneBuilder agree, and the dummy arm needs no version-matched twin."""
    import numpy as np
    from SPACE2TCR.util import parse_antibodies, reg_def
    rows = list(csv.DictReader(open(a.seqs), delimiter="\t"))[:a.n]
    tb = _builder()
    os.makedirs(a.tmp, exist_ok=True)
    sel = np.concatenate(reg_def["CDR_B_all"])
    print(f"{'id':12}{'CDRB Ca RMSD vs cohort PDB':>28}", file=sys.stderr)
    vals = []
    for r in rows:
        p = os.path.join(a.tmp, r["id"] + ".pdb")
        if not os.path.exists(p):
            tb.predict({"A": r["alpha"], "B": r["beta"]}).save(p)
        old = os.path.join(a.pdb_dir, r["id"] + ".pdb")
        m_new, m_old = parse_antibodies([p, old], n_jobs=1)
        # Same model, same numbering: compare the beta CDR CA coordinates
        # directly, with no superposition, after matching on IMGT number.
        def coords(m):
            num, xyz = m[0], m[1]
            keep = {int(x) for x in sel}
            return {int(n): xyz[i] for i, n in enumerate(num) if int(n) in keep}
        ca_n, ca_o = coords(m_new), coords(m_old)
        common = sorted(set(ca_n) & set(ca_o))
        if len(common) < 10:
            print(f"{r['id']:12}{'too few shared residues':>28}", file=sys.stderr)
            continue
        d = np.array([ca_n[k] - ca_o[k] for k in common])
        v = float(np.sqrt((d ** 2).sum(axis=1).mean()))
        vals.append(v)
        print(f"{r['id']:12}{v:>28.3f}", file=sys.stderr)
    if vals:
        import statistics
        print(f"\nmedian {statistics.median(vals):.3f} A over {len(vals)} "
              f"receptors (no superposition; coordinates compared as built)",
              file=sys.stderr)


def cmd_compare(a):
    """Per-receptor beta-CDR Ca RMSD between two directories of models.

    Run it twice to get the two numbers that matter, on the same receptors:

      NOISE FLOOR   real-alpha re-model  vs  the cohort PDB
                    how much the predictor moves when nothing changed
      ALPHA EFFECT  dummy-alpha re-model vs  real-alpha re-model
                    how much the beta moves when only the partner changed

    If the alpha effect sits inside the noise floor, the beta coordinates do
    not depend on which alpha they were modelled with, and there is no
    alpha leak to find. If it is clearly larger, the
    leak is real and this is its size, in the same units as the 0.75 A
    clustering threshold.

    Coordinates are compared as built, matched on IMGT number, with no
    superposition: both models come from the same predictor on the same beta
    sequence, so a common frame is already theirs and superposing would hide
    exactly the rigid-body shift that an interface change would cause."""
    import numpy as np, statistics
    from SPACE2TCR.util import parse_antibodies, reg_def
    sel = {int(x) for x in np.concatenate(reg_def["CDR_B_all"])}
    ids = sorted(f[:-4] for f in os.listdir(a.dir_a) if f.endswith(".pdb"))
    ids = [i for i in ids if os.path.exists(os.path.join(a.dir_b, i + ".pdb"))]
    if not ids:
        sys.exit("no receptors present in both directories")
    rows, vals = [], []
    for i in ids:
        ma, mb = parse_antibodies([os.path.join(a.dir_a, i + ".pdb"),
                                   os.path.join(a.dir_b, i + ".pdb")], n_jobs=1)
        def ca(m):
            return {int(n): m[1][k] for k, n in enumerate(m[0]) if int(n) in sel}
        A, B = ca(ma), ca(mb)
        common = sorted(set(A) & set(B))
        if len(common) < 10:
            continue
        d = np.array([A[k] - B[k] for k in common])
        v = float(np.sqrt((d ** 2).sum(axis=1).mean()))
        rows.append({"id": i, "n_residues": len(common), "rmsd": round(v, 4)})
        vals.append(v)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "n_residues", "rmsd"])
        w.writeheader(); w.writerows(rows)
    q = statistics.quantiles(vals, n=4) if len(vals) > 3 else [float("nan")] * 3
    print(f"{a.label}: n={len(vals)}  median {statistics.median(vals):.3f} A  "
          f"q25 {q[0]:.3f}  q75 {q[2]:.3f}  max {max(vals):.3f}", file=sys.stderr)
    print(f"  fraction over 0.75 A: {sum(1 for v in vals if v > 0.75)/len(vals):.3f}",
          file=sys.stderr)
    print(f"wrote {a.out}", file=sys.stderr)


def cmd_bins(a):
    """Pick whole beta-length bins, largest first, up to --target receptors."""
    from SPACE2TCR.util import parse_antibodies, reg_def
    keep = {r["id"] for r in csv.DictReader(open(a.seqs), delimiter="\t")}
    files = sorted(os.path.join(a.pdb_dir, f) for f in os.listdir(a.pdb_dir)
                   if f.endswith(".pdb") and f[:-4] in keep)
    models = parse_antibodies(files, n_jobs=a.n_jobs)
    ids = [os.path.basename(f)[:-4] for f in files]
    bins = collections.defaultdict(list)
    for k, m in enumerate(models):
        present = {int(x) for x in m[0]}
        key = tuple(sum(1 for r in reg_def[rn] if int(r) in present)
                    for rn in ("CDRB1", "CDRB2", "CDRB3"))
        bins[key].append(ids[k])
    chosen, total = [], 0
    for key, mem in sorted(bins.items(), key=lambda kv: -len(kv[1])):
        if total >= a.target:
            break
        chosen.append((key, mem)); total += len(mem)
    with open(a.out_ids, "w") as fh:
        for key, mem in chosen:
            for i in mem:
                fh.write(i + "\n")
    print(f"{len(bins)} beta bins over {len(ids)} receptors; took the "
          f"{len(chosen)} largest = {total} receptors", file=sys.stderr)
    for key, mem in chosen:
        print(f"   CDRB lengths {key}  n={len(mem)}", file=sys.stderr)
    print(f"wrote {a.out_ids}", file=sys.stderr)


def cmd_model(a):
    rows = {r["id"]: r for r in csv.DictReader(open(a.seqs), delimiter="\t")}
    ids = [ln.strip() for ln in open(a.ids) if ln.strip()]
    missing = [i for i in ids if i not in rows]
    if missing:
        sys.exit(f"{len(missing)} ids have no sequence, e.g. {missing[:3]}")
    # The fixed partner. A real alpha V-domain from the cohort, so the
    # predictor stays in distribution; which one is arbitrary and is recorded
    # in the output directory so the run is reproducible.
    dummy = rows[a.dummy_id]["alpha"] if a.dummy_id else rows[ids[0]]["alpha"]
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "_PROVENANCE.txt"), "w") as fh:
        fh.write(f"mode={a.mode}\nn={len(ids)}\n"
                 f"dummy_alpha_from={a.dummy_id or ids[0]}\n"
                 f"dummy_alpha={dummy}\n")
    tb = _builder()
    t0, done = time.time(), 0
    for i in ids:
        p = os.path.join(a.out_dir, i + ".pdb")
        if os.path.exists(p):
            done += 1; continue
        alpha = dummy if a.mode == "dummy" else rows[i]["alpha"]
        try:
            tb.predict({"A": alpha, "B": rows[i]["beta"]}).save(p)
        except Exception as e:
            print(f"  FAIL {i}: {type(e).__name__} {e}", file=sys.stderr)
            continue
        done += 1
        if done % 10 == 0:
            el = time.time() - t0
            print(f"  {done}/{len(ids)}  {el/60:.1f} min elapsed, "
                  f"~{(len(ids)-done)*el/max(done,1)/60:.0f} min left",
                  file=sys.stderr, flush=True)
    print(f"wrote {done} models to {a.out_dir}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("seqs");   p.set_defaults(fn=cmd_seqs)
    p.add_argument("--pdb-dir", required=True); p.add_argument("--out", required=True)

    p = sub.add_parser("verify"); p.set_defaults(fn=cmd_verify)
    p.add_argument("--seqs", required=True); p.add_argument("--pdb-dir", required=True)
    p.add_argument("-n", type=int, default=10)
    p.add_argument("--tmp", default="/tmp/verify_remodel")

    p = sub.add_parser("compare"); p.set_defaults(fn=cmd_compare)
    p.add_argument("--dir-a", required=True); p.add_argument("--dir-b", required=True)
    p.add_argument("--out", required=True); p.add_argument("--label", default="compare")

    p = sub.add_parser("bins");   p.set_defaults(fn=cmd_bins)
    p.add_argument("--pdb-dir", required=True); p.add_argument("--seqs", required=True)
    p.add_argument("--target", type=int, default=400)
    p.add_argument("--out-ids", required=True); p.add_argument("--n-jobs", type=int, default=-1)

    p = sub.add_parser("model");  p.set_defaults(fn=cmd_model)
    p.add_argument("--seqs", required=True); p.add_argument("--ids", required=True)
    p.add_argument("--mode", choices=["dummy", "real"], required=True)
    p.add_argument("--dummy-id", default=None)
    p.add_argument("--out-dir", required=True)

    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
