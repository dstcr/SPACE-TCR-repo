#!/usr/bin/env python3
"""
extract_dihedrals.py  —  Stage D (optional, runs after Stage C / make_bins)

Extracts raw CDR3 backbone dihedral angles (phi, psi) for every TCR in the
PDB directory and writes them to a single TSV.  This replaces the composite
distance mode: instead of using dihedrals as a clustering criterion, you get
them as a per-TCR feature table that can be used downstream however you like
(e.g. principal components, distance matrices, visualisation, annotation of
cluster outputs).

Output columns
--------------
tcr_id          : PDB stem (e.g. tcr1)
chain           : A (alpha) or B (beta)
resi            : IMGT residue number (105-117 by default)
phi_rad         : phi angle in radians  (-pi to pi, 0.0 if atoms missing)
psi_rad         : psi angle in radians  (-pi to pi, 0.0 if atoms missing)
phi_deg         : phi in degrees
psi_deg         : psi in degrees
atoms_present   : comma-separated list of backbone atoms found at this residue

One row per (tcr_id, chain, residue).  A TCR with CDR3 residues 105-117 on
both chains will produce up to 26 rows (13 residues x 2 chains).

Usage
-----
python extract_dihedrals.py \
    --pdb-dir  /path/to/pdb \
    --out-tsv  /path/to/dihedrals.tsv.gz \
    --cdr3-start 105 \
    --cdr3-end   117 \
    --workers  8

Then in R:
    dihedrals <- read_tsv("dihedrals.tsv.gz")
"""

from __future__ import annotations

import argparse
import gzip
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np


# ── Backbone parsing ──────────────────────────────────────────────────────────

def parse_backbone(pdb_path: Path,
                   chains: Tuple[str, ...] = ("A", "B"),
                   max_resi: int = 128,
                   ) -> Dict[Tuple[str, int], Dict[str, np.ndarray]]:
    """
    Parse N, CA, C backbone atoms from a PDB file.
    Returns dict: (chain, resi) -> {'N': array(3,), 'CA': array(3,), 'C': array(3,)}
    Missing atoms simply absent from inner dict.
    """
    atoms: Dict[Tuple[str, int], Dict[str, np.ndarray]] = {}
    try:
        text = pdb_path.read_text(errors="ignore")
    except OSError:
        return atoms

    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name not in ("N", "CA", "C"):
            continue
        chain = line[21:22].strip()
        if chain not in chains:
            continue
        try:
            resi = int(line[22:26])
            if resi <= 0 or resi > max_resi:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue

        key = (chain, resi)
        if key not in atoms:
            atoms[key] = {}
        atoms[key][atom_name] = np.array([x, y, z], dtype=np.float64)

    return atoms


# ── Dihedral calculation ──────────────────────────────────────────────────────

def _dihedral(p0: np.ndarray, p1: np.ndarray,
               p2: np.ndarray, p3: np.ndarray) -> float:
    """
    Dihedral angle (radians, -pi..pi) for four 3D points using atan2.
    Returns 0.0 if any bond vector is degenerate.
    """
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    n1 = np.linalg.norm(b1)
    if n1 < 1e-8:
        return 0.0
    b1h = b1 / n1

    v = b0 - np.dot(b0, b1h) * b1h
    w = b2 - np.dot(b2, b1h) * b1h

    x = np.dot(v, w)
    y = np.dot(np.cross(b1h, v), w)
    return float(np.arctan2(y, x))


# ── Per-TCR extraction ────────────────────────────────────────────────────────

def extract_dihedrals_for_pdb(
        pdb_path: Path,
        cdr3_residues: List[int],
) -> List[dict]:
    """
    Extract phi/psi for each (chain, residue) in cdr3_residues.
    Returns a list of row dicts; empty list if PDB cannot be parsed.
    """
    backbone = parse_backbone(pdb_path)
    if not backbone:
        return []

    tcr_id = pdb_path.stem
    rows: List[dict] = []

    for chain in ("A", "B"):
        for k, resi in enumerate(cdr3_residues):
            curr = backbone.get((chain, resi), {})
            prev = backbone.get((chain, cdr3_residues[k - 1]), {}) if k > 0 else {}
            nxt  = backbone.get((chain, cdr3_residues[k + 1]), {}) \
                   if k < len(cdr3_residues) - 1 else {}

            atoms_present = sorted(curr.keys())

            # phi = C(i-1) — N(i) — CA(i) — C(i)
            if "C" in prev and "N" in curr and "CA" in curr and "C" in curr:
                phi = _dihedral(prev["C"], curr["N"], curr["CA"], curr["C"])
            else:
                phi = 0.0

            # psi = N(i) — CA(i) — C(i) — N(i+1)
            if "N" in curr and "CA" in curr and "C" in curr and "N" in nxt:
                psi = _dihedral(curr["N"], curr["CA"], curr["C"], nxt["N"])
            else:
                psi = 0.0

            rows.append({
                "tcr_id":        tcr_id,
                "chain":         chain,
                "resi":          resi,
                "phi_rad":       round(phi, 6),
                "psi_rad":       round(psi, 6),
                "phi_deg":       round(math.degrees(phi), 4),
                "psi_deg":       round(math.degrees(psi), 4),
                "atoms_present": ",".join(atoms_present) if atoms_present else "",
            })

    return rows


# ── Worker (used by ProcessPoolExecutor) ──────────────────────────────────────

def _worker(args: Tuple[Path, List[int]]) -> List[dict]:
    pdb_path, cdr3_residues = args
    return extract_dihedrals_for_pdb(pdb_path, cdr3_residues)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Extract raw CDR3 backbone dihedral angles (phi, psi) for every PDB "
            "in a directory and write a flat TSV.  One row per (TCR, chain, residue)."
        )
    )
    ap.add_argument(
        "--pdb-dir", required=True,
        help="Directory containing PDB files produced by Stage B (TCRBuilder2).",
    )
    ap.add_argument(
        "--out-tsv", required=True,
        help="Output path for the dihedral table (.tsv or .tsv.gz).",
    )
    ap.add_argument(
        "--cdr3-start", type=int, default=105,
        help="First IMGT residue number of the CDR3 loop (default: 105).",
    )
    ap.add_argument(
        "--cdr3-end", type=int, default=117,
        help="Last IMGT residue number of the CDR3 loop inclusive (default: 117).",
    )
    ap.add_argument(
        "--workers", type=int, default=4,
        help="Number of parallel processes (default: 4).",
    )
    ap.add_argument(
        "--progress-every", type=int, default=1000,
        help="Log a progress line every N PDBs (default: 1000).",
    )
    args = ap.parse_args()

    pdb_dir = Path(args.pdb_dir).expanduser().resolve()
    out_tsv = Path(args.out_tsv).expanduser().resolve()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    cdr3_residues = list(range(args.cdr3_start, args.cdr3_end + 1))

    pdbs = sorted(pdb_dir.glob("*.pdb"))
    if not pdbs:
        raise SystemExit(f"[FATAL] No *.pdb files found in {pdb_dir}")

    print(f"[INFO] PDBs found: {len(pdbs):,}", flush=True)
    print(f"[INFO] CDR3 residues: {cdr3_residues[0]}–{cdr3_residues[-1]} "
          f"({len(cdr3_residues)} positions, {len(cdr3_residues)*2} angles per chain)",
          flush=True)
    print(f"[INFO] Workers: {args.workers}", flush=True)
    print(f"[INFO] Output: {out_tsv}", flush=True)

    header = "tcr_id\tchain\tresi\tphi_rad\tpsi_rad\tphi_deg\tpsi_deg\tatoms_present\n"

    open_fn = gzip.open if str(out_tsv).endswith(".gz") else open
    t0      = time.time()
    done    = 0
    written = 0
    errors  = 0

    with open_fn(out_tsv, "wt") as fh:
        fh.write(header)

        work_items = [(p, cdr3_residues) for p in pdbs]

        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for rows in pool.map(_worker, work_items, chunksize=50):
                done += 1
                if not rows:
                    errors += 1
                else:
                    for r in rows:
                        fh.write(
                            f"{r['tcr_id']}\t{r['chain']}\t{r['resi']}\t"
                            f"{r['phi_rad']}\t{r['psi_rad']}\t"
                            f"{r['phi_deg']}\t{r['psi_deg']}\t"
                            f"{r['atoms_present']}\n"
                        )
                        written += 1

                if args.progress_every > 0 and (
                    done % args.progress_every == 0 or done == len(pdbs)
                ):
                    elapsed = time.time() - t0
                    rate    = done / elapsed if elapsed > 0 else float("nan")
                    print(
                        f"[INFO] {done:,}/{len(pdbs):,} PDBs "
                        f"({done/len(pdbs):.1%})  "
                        f"rows_written={written:,}  errors={errors:,}  "
                        f"elapsed={elapsed:,.0f}s  rate={rate:,.0f}/s",
                        flush=True,
                    )

    print(f"\n[DONE] rows written: {written:,}  errors (empty/bad PDB): {errors:,}")
    print(f"[DONE] output: {out_tsv}")
    print(
        f"\nTo load in R:\n"
        f"  library(readr)\n"
        f"  dihedrals <- read_tsv('{out_tsv}')\n"
        f"\nColumns: tcr_id, chain (A/B), resi, phi_rad, psi_rad, phi_deg, psi_deg, atoms_present"
    )


if __name__ == "__main__":
    main()
