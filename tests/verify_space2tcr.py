"""Verify the vendored SPACE2TCR against the interface the pipeline needs,
and compute a real RMSD end to end.

Run from the repo root:  PYTHONPATH=. python3 -I verify_space2tcr.py
"""
import sys, os, inspect
import numpy as np

REPO = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
sys.path.insert(0, REPO)

fails = []
def ok(label, cond, extra=None):
    txt = ""
    if extra is not None:
        txt = str(extra).replace("\n", " ")
        txt = ("  " + (txt[:90] + "..." if len(txt) > 90 else txt))
    print(("  PASS  " if bool(cond) else "  FAIL  ") + label + txt)
    if not bool(cond):
        fails.append(label)

print("== 1. the six names the pipeline imports ==")
from SPACE2TCR.util import (rmsd, reg_def, reg_def_fw_all, reg_def_CDR_all,
                            cluster_antibodies_by_CDR_length, same_length_only)
ok("all six import from SPACE2TCR.util", True)
ok("rmsd is callable", callable(rmsd), "(" + ", ".join(inspect.signature(rmsd).parameters) + ")")
ok("cluster_antibodies_by_CDR_length is callable",
   callable(cluster_antibodies_by_CDR_length))
ok("same_length_only is an ndarray of 6 ones (NOT a callable)",
   isinstance(same_length_only, np.ndarray) and same_length_only.tolist() == [1]*6,
   same_length_only)

print("\n== 2. reg_def: the TCR region definitions ==")
keys = sorted(reg_def.keys())
ok("reg_def has the keys the scripts use",
   all(k in reg_def for k in ["CDR_all", "CDR_A_all", "CDR_B_all", "fw_all"]), keys)
ok("CDR_all == CDR_B_all + CDR_A_all (the make_bins fallback holds)",
   len(reg_def["CDR_all"]) == len(reg_def["CDR_B_all"]) + len(reg_def["CDR_A_all"])
   and all(np.array_equal(a, b) for a, b in
           zip(reg_def["CDR_all"], reg_def["CDR_B_all"] + reg_def["CDR_A_all"])))
ok("CDR_all is 6 parts (the padded selection make_bins builds)",
   len(reg_def["CDR_all"]) == 6)
ok("reg_def_CDR_all flattens to all CDR residues",
   reg_def_CDR_all.ndim == 1 and len(reg_def_CDR_all) == sum(len(x) for x in reg_def["CDR_all"]),
   f"{len(reg_def_CDR_all)} residues")
ok("reg_def_fw_all flattens to the framework anchors",
   reg_def_fw_all.ndim == 1, f"{len(reg_def_fw_all)} residues")
print(f"      CDRB3 = {reg_def['CDRB3'].min()}-{reg_def['CDRB3'].max()}  "
      f"CDRA3 = {reg_def['CDRA3'].min()}-{reg_def['CDRA3'].max()}")
ok("alpha is offset by +128 from beta (paired-chain numbering)",
   (reg_def["CDRA3"] - reg_def["CDRB3"]).tolist() == [128]*13)
ok("make_bins --selection-key choices are non-empty", len(keys) > 0, len(keys))

print("\n== 3. a real RMSD through the real code path ==")
# Build two minimal paired-chain PDBs: chain B numbered 1..128, chain A 1..128.
# get_antibody() offsets chain A by +128 and discards any residue numbered >128.
def pdb_text(coords_b, coords_a):
    lines = []
    n = 1
    for ch, coords in (("B", coords_b), ("A", coords_a)):
        for i, (x, y, z) in enumerate(coords, start=1):
            lines.append(
                f"ATOM  {n:5d}  CA  ALA {ch}{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
            n += 1
    return "\n".join(lines) + "\n"

rng = np.random.default_rng(0)
# 128 residues per chain so every CDR and framework position exists.
base_b = rng.normal(0, 10, (128, 3))
base_a = rng.normal(0, 10, (128, 3))

from SPACE2TCR.util import get_antibody, parse_antibody, get_CDR_lengths

ab1 = get_antibody(pdb_text(base_b, base_a))
ok("get_antibody parses CA atoms", ab1[0].size > 0, f"{ab1[0].size} CA atoms")
ok("numbering spans both chains with the offset",
   ab1[0].max() > 128 and ab1[0].min() >= 1, f"{ab1[0].min()}-{ab1[0].max()}")

# Identical structure -> RMSD 0.
r_same = rmsd(ab1, ab1, selection=reg_def_CDR_all, anchors=reg_def_fw_all)
ok("RMSD of a structure against itself is 0", abs(r_same) < 1e-8, r_same)

# Rigid-body transform of the whole structure -> RMSD still ~0 after alignment
# on the framework anchors. This is the property the method depends on: the
# distance must reflect CDR geometry, not model orientation.
from SPACE2TCR.util import random_rot
R = random_rot()
ab2 = ((ab1[0]).copy(), (ab1[1] @ R.T) + np.array([13.0, -7.0, 21.0]))
r_rot = rmsd(ab1, ab2, selection=reg_def_CDR_all, anchors=reg_def_fw_all)
ok("RMSD is invariant to rigid-body rotation+translation", abs(r_rot) < 1e-6, r_rot)

# Perturb only the CDR loops -> RMSD must rise.
pert_b = base_b.copy()
cdrb3 = reg_def["CDRB3"]
pert_b[cdrb3 - 1] += 2.0        # residue number n is row n-1 within the chain
ab3 = get_antibody(pdb_text(pert_b, base_a))
r_pert = rmsd(ab1, ab3, selection=reg_def_CDR_all, anchors=reg_def_fw_all)
ok("perturbing CDR3B raises the RMSD", r_pert > 0.1, f"{r_pert:.4f}")

# Perturb only framework residues far from any CDR -> CDR RMSD should stay small.
fw_only = np.setdiff1d(reg_def["fwB"], reg_def_CDR_all)
pert_fw = base_b.copy()
pert_fw[fw_only[:10] - 1] += 0.05   # tiny, so the anchor fit barely moves
ab4 = get_antibody(pdb_text(pert_fw, base_a))
r_fw = rmsd(ab1, ab4, selection=reg_def_CDR_all, anchors=reg_def_fw_all)
ok("a small framework-only perturbation keeps CDR RMSD below the CDR one",
   r_fw < r_pert, f"{r_fw:.4f} < {r_pert:.4f}")

print("\n== 4. binning, as make_bins.py calls it ==")
# make_bins pads the selection to 6 parts and passes clustering="bins".
sel = [np.asarray(x).ravel().astype(int) for x in reg_def["CDR_all"]]
abs_ = [ab1, ab3, ab4]
ids = ["t1", "t3", "t4"]
cdr_clusters, cdr_ids = cluster_antibodies_by_CDR_length(
    abs_, ids, selection=sel, clustering="bins", tolerance=same_length_only)
ok("binning returns a dict of bin key -> ids", isinstance(cdr_ids, dict), list(cdr_ids.keys()))
ok("every input id is placed in some bin",
   sorted(i for v in cdr_ids.values() for i in v) == sorted(ids))
ok("equal-length structures land in one bin (same_length_only)",
   len(cdr_ids) == 1, f"{len(cdr_ids)} bin(s)")

print("\n== 5. the pipeline's own Stage C scripts import cleanly ==")
import importlib.util
for mod in ("make_bins", "cluster_one_bin"):
    p = os.path.join(REPO, "scripts", mod + ".py")
    try:
        spec = importlib.util.spec_from_file_location(mod, p)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        ok(f"scripts/{mod}.py imports", True)
    except Exception as e:
        ok(f"scripts/{mod}.py imports", False, f"{type(e).__name__}: {e}")

print()
if fails:
    print(f"==== {len(fails)} FAILURE(S): {fails} ====")
    sys.exit(1)
print("==== all checks passed — the repository is self-contained ====")
