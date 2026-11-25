#!/usr/bin/env python3
"""
cluster_one_bin.py

Cluster a single predefined length bin from bins.jsonl and, if requested,
write a full pairwise distance table for every attempted comparison in the bin.

Behaviour
---------
- Reads one bin record from bins.jsonl
- Parses PDBs from the filesystem (--pdb-base-dir) or from a tar archive
- Uses exact connected-components clustering when n <= max_exact_bin
- Uses greedy representative clustering when n > max_exact_bin
- With --write-distance-table, writes a gzipped TSV of every attempted pairwise
  comparison in the bin (status + RMSD / dihedral / composite)

Distance modes
--------------
--distance-mode rmsd        (default) as SPACE2: cluster on RMSD alone
--distance-mode dihedral    cluster on mean CDR3 backbone dihedral distance alone
--distance-mode composite   weighted combination of normalised RMSD + dihedral distance

Composite weights and scaling:
  --w-rmsd          weight on normalised RMSD          (default 0.5)
  --w-dihedral      weight on normalised dihedral dist  (default 0.5)
  --rmsd-scale      Angstroms used to normalise RMSD    (default 3.0)

A complete pairwise distance table needs exact mode, so --write-distance-table
with a bin larger than --max-exact-bin stops with an error.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
import tarfile
import time
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from SPACE2TCR.util import rmsd, reg_def, reg_def_fw_all, reg_def_CDR_all

# scripts/ is not a package, so make the sibling module importable whether this
# file is run as a script, via sbatch, or from the repository root.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dihedral_util import dihedral_distance as _shared_dihedral_distance


# ============================================================
# Logging
# ============================================================

def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _log(msg: str) -> None:
    print(f"[{_now()}] {msg}", flush=True)


# ============================================================
# PDB parsing — Cα (for RMSD) and backbone N/CA/C (for dihedrals)
# ============================================================

def _parse_pdb_lines_ca(
    lines: Iterable[str],
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Parse CA atoms for chains A/B, resSeq 1..128.
    Returns (numbers, coords) or None if parsing fails / no usable atoms found.
    This is the original parser used by the SPACE2 rmsd() function.
    """
    numbers: List[int] = []
    coords: List[Tuple[float, float, float]] = []

    for line in lines:
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue

        chain = line[21].strip()
        if chain not in ("A", "B"):
            continue

        try:
            resi = int(line[22:26])
            if resi <= 0 or resi > 128:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            parts = line.split()
            try:
                chain2 = parts[4]
                if chain2 not in ("A", "B"):
                    continue
                resi = int(parts[5])
                if resi <= 0 or resi > 128:
                    continue
                x = float(parts[6])
                y = float(parts[7])
                z = float(parts[8])
                chain = chain2
            except (ValueError, IndexError):
                return None

        idx = resi + (128 if chain == "A" else 0)
        numbers.append(idx)
        coords.append((x, y, z))

    if not numbers:
        return None

    return np.asarray(numbers, dtype=int), np.asarray(coords, dtype=np.float32)


def _parse_pdb_backbone(
    lines: Iterable[str],
    chains: Tuple[str, ...] = ("A", "B"),
    max_resi: int = 128,
) -> Dict[Tuple[str, int], Dict[str, np.ndarray]]:
    """
    Parse N, CA, C backbone atoms from PDB lines.

    Returns a dict keyed by (chain, resi) -> {'N': array, 'CA': array, 'C': array}.
    Used for dihedral angle extraction.  Separate from the CA-only parser so
    that the original SPACE2 rmsd() call is unaffected.
    """
    atoms: Dict[Tuple[str, int], Dict[str, np.ndarray]] = {}

    for line in lines:
        if not line.startswith("ATOM"):
            continue

        atom_name = line[12:16].strip()
        if atom_name not in ("N", "CA", "C"):
            continue

        chain = line[21].strip()
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
        atoms[key][atom_name] = np.array([x, y, z], dtype=np.float32)

    return atoms


# ============================================================
# Dihedral angle computation
# ============================================================

# CDR3 residue numbers in IMGT numbering (positions 105-117).
# Chain A = alpha, chain B = beta.
_CDR3_IMGT = list(range(105, 118))


def _dihedral_4pts(
    p0: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
) -> float:
    """
    Dihedral angle (radians, range -pi..pi) defined by four 3D points.
    Uses the Gram-Schmidt / atan2 formula, which is numerically stable.
    """
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    norm_b1 = np.linalg.norm(b1)
    if norm_b1 < 1e-8:
        return 0.0
    b1_hat = b1 / norm_b1

    # Project b0 and b2 onto plane perpendicular to b1
    v = b0 - np.dot(b0, b1_hat) * b1_hat
    w = b2 - np.dot(b2, b1_hat) * b1_hat

    x = np.dot(v, w)
    y = np.dot(np.cross(b1_hat, v), w)

    return float(np.arctan2(y, x))


def _extract_chain_dihedrals(
    backbone: Dict[Tuple[str, int], Dict[str, np.ndarray]],
    chain: str,
    residues: List[int],
) -> List[float]:
    """
    Extract phi and psi angles for a list of residue numbers on one chain.

    phi(i) = C(i-1) -- N(i) -- CA(i) -- C(i)
    psi(i) = N(i)   -- CA(i) -- C(i) -- N(i+1)

    Missing atoms are filled with 0.0 (neutral / no-information).
    Returns a flat list [phi_r0, psi_r0, phi_r1, psi_r1, ...].
    """
    angles: List[float] = []

    for k, resi in enumerate(residues):
        curr = backbone.get((chain, resi), {})
        prev = backbone.get((chain, residues[k - 1]), {}) if k > 0 else {}
        nxt  = backbone.get((chain, residues[k + 1]), {}) if k < len(residues) - 1 else {}

        # phi
        if "C" in prev and "N" in curr and "CA" in curr and "C" in curr:
            phi = _dihedral_4pts(prev["C"], curr["N"], curr["CA"], curr["C"])
        else:
            phi = 0.0

        # psi
        if "N" in curr and "CA" in curr and "C" in curr and "N" in nxt:
            psi = _dihedral_4pts(curr["N"], curr["CA"], curr["C"], nxt["N"])
        else:
            psi = 0.0

        angles.extend([phi, psi])

    return angles


def _extract_cdr3_dihedrals(
    backbone: Dict[Tuple[str, int], Dict[str, np.ndarray]],
    cdr3_residues: List[int] = _CDR3_IMGT,
) -> Tuple[List[float], List[float]]:
    """
    Extract CDR3 backbone dihedral angles for both chains.

    Returns:
        dihedrals_alpha: flat [phi, psi] list for chain A CDR3
        dihedrals_beta:  flat [phi, psi] list for chain B CDR3
    """
    da = _extract_chain_dihedrals(backbone, "A", cdr3_residues)
    db = _extract_chain_dihedrals(backbone, "B", cdr3_residues)
    return da, db


def _circular_dist(a: float, b: float) -> float:
    """
    Circular distance between two angles in radians.
    Result is always in [0, pi].
    """
    diff = abs(a - b) % (2.0 * np.pi)
    return float(min(diff, 2.0 * np.pi - diff))


# Single-sourced from scripts/dihedral_util.py, so there is one Python copy;
# tests/test_dihedral_distance.R checks it against the R definition.
def _dihedral_vector_distance(
    angles_i: List[float],
    angles_j: List[float],
) -> float:
    """Mean circular distance between two dihedral-angle vectors, in RADIANS.

    Length mismatches pad the shorter vector with pi (maximum disagreement).
    Returns a value in [0, pi].
    """
    return float(_shared_dihedral_distance(angles_i, angles_j))


def _composite_dist(
    rmsd_val: float,
    dd_alpha: float,
    dd_beta: float,
    w_rmsd: float,
    w_dihedral: float,
    rmsd_scale: float,
) -> float:
    """
    Weighted sum of normalised RMSD and mean normalised CDR3 dihedral distance.

    Both components are clipped to [0, 1] before weighting so the composite
    is also in [0, 1] (assuming w_rmsd + w_dihedral == 1).

    rmsd_scale  — Angstrom value mapped to 1.0 (e.g. 3.0)
    dihedral is normalised by pi (its theoretical maximum).
    """
    rmsd_norm     = min(rmsd_val / rmsd_scale, 1.0)
    dihedral_mean = (dd_alpha + dd_beta) / 2.0
    dihedral_norm = dihedral_mean / np.pi   # pi is the maximum possible

    return w_rmsd * rmsd_norm + w_dihedral * dihedral_norm


# ============================================================
# Helpers
# ============================================================

def _clean_error_message(msg: object, max_len: int = 300) -> str:
    s = str(msg).replace("\t", " ").replace("\n", " ").strip()
    if len(s) > max_len:
        s = s[: max_len - 3] + "..."
    return s


def _tar_id(tar_path: str, member_name: str) -> str:
    return f"tar:{str(Path(tar_path).expanduser().resolve())}::{member_name}"


def _iter_tar_member_lines(tf: tarfile.TarFile, member_name: str) -> Iterator[str]:
    try:
        fh = tf.extractfile(member_name)
    except Exception:
        fh = None
    if fh is None:
        return
    for raw in fh:
        yield raw.decode("utf-8", errors="ignore")


def _resolve_selection(selection_key: str) -> np.ndarray:
    if selection_key == "CDR_all":
        return np.asarray(reg_def_CDR_all).ravel().astype(int)
    x = reg_def[selection_key]
    if isinstance(x, dict):
        x = x[sorted(x.keys())[0]]
    if isinstance(x, (list, tuple)):
        return np.concatenate([np.asarray(xx).ravel() for xx in x]).astype(int)
    return np.asarray(x).ravel().astype(int)


def _resolve_anchors(anchors_key: str) -> np.ndarray:
    if anchors_key == "fw_all":
        return np.asarray(reg_def_fw_all).ravel().astype(int)
    x = reg_def[anchors_key]
    if isinstance(x, dict):
        x = x[sorted(x.keys())[0]]
    if isinstance(x, (list, tuple)):
        return np.concatenate([np.asarray(xx).ravel() for xx in x]).astype(int)
    return np.asarray(x).ravel().astype(int)


def _pair_chunks(n: int, chunk_pairs: int) -> Iterator[List[Tuple[int, int]]]:
    buf: List[Tuple[int, int]] = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            buf.append((i, j))
            if len(buf) >= chunk_pairs:
                yield buf
                buf = []
    if buf:
        yield buf


# ============================================================
# Worker — runs in subprocess via multiprocessing fork
# ============================================================

_G_CLUSTER          = None   # List of (numbers, coords) CA models
_G_BACKBONES        = None   # List of backbone dicts (N/CA/C)
_G_SELECTION        = None
_G_ANCHORS          = None
_G_CUTOFF           = None
_G_PRUNE            = None
_G_WRITE_DIST_ROWS  = None
_G_DISTANCE_MODE    = None   # "rmsd" | "dihedral" | "composite"
_G_W_RMSD           = None
_G_W_DIHEDRAL       = None
_G_RMSD_SCALE       = None
_G_CDR3_RESIDUES    = None


def _init_worker(
    cluster,
    backbones,
    selection,
    anchors,
    cutoff,
    prune,
    write_distance_rows,
    distance_mode,
    w_rmsd,
    w_dihedral,
    rmsd_scale,
    cdr3_residues,
) -> None:
    global _G_CLUSTER, _G_BACKBONES, _G_SELECTION, _G_ANCHORS
    global _G_CUTOFF, _G_PRUNE, _G_WRITE_DIST_ROWS, _G_DISTANCE_MODE
    global _G_W_RMSD, _G_W_DIHEDRAL, _G_RMSD_SCALE, _G_CDR3_RESIDUES
    _G_CLUSTER         = cluster
    _G_BACKBONES       = backbones
    _G_SELECTION       = selection
    _G_ANCHORS         = anchors
    _G_CUTOFF          = cutoff
    _G_PRUNE           = prune
    _G_WRITE_DIST_ROWS = write_distance_rows
    _G_DISTANCE_MODE   = distance_mode
    _G_W_RMSD          = w_rmsd
    _G_W_DIHEDRAL      = w_dihedral
    _G_RMSD_SCALE      = rmsd_scale
    _G_CDR3_RESIDUES   = cdr3_residues


def _process_pair_chunk(
    pairs: List[Tuple[int, int]],
) -> Tuple[
    List[Tuple[int, int, str, Optional[float], Optional[float], Optional[float], Optional[float], str]],
    List[Tuple[int, int]],
]:
    """
    Process a chunk of (i, j) index pairs.

    For each pair:
      - compute RMSD (always needed for rmsd/composite modes)
      - compute CDR3 dihedral distance (needed for dihedral/composite modes)
      - derive the clustering distance according to the selected mode
      - record a distance table row if requested
      - record a graph edge if distance <= threshold

    Distance table row format:
        (i, j, status, rmsd, dihedral_alpha, dihedral_beta, clustering_dist, error)
    """
    distance_rows = []
    edges: List[Tuple[int, int]] = []

    mode     = _G_DISTANCE_MODE
    w_rmsd   = _G_W_RMSD
    w_dih    = _G_W_DIHEDRAL
    r_scale  = _G_RMSD_SCALE
    residues = _G_CDR3_RESIDUES

    edge_thresh = float(_G_CUTOFF)
    if _G_PRUNE is not None:
        edge_thresh = min(edge_thresh, float(_G_PRUNE))

    for i, j in pairs:

        # ── RMSD ──────────────────────────────────────────────────────────────
        rmsd_val: Optional[float] = None

        if mode in ("rmsd", "composite"):
            try:
                r = rmsd(
                    _G_CLUSTER[i], _G_CLUSTER[j],
                    selection=_G_SELECTION,
                    anchors=_G_ANCHORS,
                )
                if not np.isfinite(r):
                    raise ValueError("non-finite RMSD")
                rmsd_val = float(r)
            except Exception as exc:
                if _G_WRITE_DIST_ROWS:
                    distance_rows.append(
                        (i, j, "error", None, None, None, None,
                         _clean_error_message(f"{type(exc).__name__}: {exc}"))
                    )
                continue

        # ── Dihedral distance ──────────────────────────────────────────────────
        dd_alpha: Optional[float] = None
        dd_beta:  Optional[float] = None

        if mode in ("dihedral", "composite"):
            try:
                da_i, db_i = _extract_cdr3_dihedrals(_G_BACKBONES[i], residues)
                da_j, db_j = _extract_cdr3_dihedrals(_G_BACKBONES[j], residues)
                dd_alpha = _dihedral_vector_distance(da_i, da_j)
                dd_beta  = _dihedral_vector_distance(db_i, db_j)
            except Exception as exc:
                if _G_WRITE_DIST_ROWS:
                    distance_rows.append(
                        (i, j, "error", rmsd_val, None, None, None,
                         _clean_error_message(f"dihedral: {type(exc).__name__}: {exc}"))
                    )
                continue

        # ── Clustering distance ────────────────────────────────────────────────
        if mode == "rmsd":
            clustering_dist = rmsd_val

        elif mode == "dihedral":
            # Mean of the alpha and beta dihedral distances, in radians, the
            # same unit as --cutoff and as CFG$focal_dihedral in the pipeline.
            # (The published results use --distance-mode rmsd; the gated
            # RMSD+dih clusters are built in R from the dihedral table.)
            clustering_dist = (dd_alpha + dd_beta) / 2.0

        else:  # composite
            clustering_dist = _composite_dist(
                rmsd_val, dd_alpha, dd_beta, w_rmsd, w_dih, r_scale
            )

        # ── Record row ────────────────────────────────────────────────────────
        if _G_WRITE_DIST_ROWS:
            distance_rows.append(
                (i, j, "ok", rmsd_val, dd_alpha, dd_beta, clustering_dist, "")
            )

        # ── Edge ──────────────────────────────────────────────────────────────
        if clustering_dist <= edge_thresh:
            edges.append((i, j))

    return distance_rows, edges


# ============================================================
# Exact connected-components clustering
# ============================================================

def _write_empty_distance_table(distance_table_path: Path) -> None:
    distance_table_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(distance_table_path, "wt") as fh:
        fh.write("i\tj\tid_i\tid_j\tstatus\trmsd\tdihedral_alpha\tdihedral_beta\tclustering_dist\terror\n")


def _connected_components_threshold(
    cluster_models: List,
    backbones: List,
    ids: List[str],
    selection: np.ndarray,
    anchors: np.ndarray,
    cutoff: float,
    prune: float,
    n_jobs: int,
    chunk_pairs: int,
    progress_every_chunks: int,
    distance_mode: str,
    w_rmsd: float,
    w_dihedral: float,
    rmsd_scale: float,
    cdr3_residues: List[int],
    distance_table_path: Optional[Path] = None,
) -> Tuple[np.ndarray, Dict[int, str]]:
    """
    Exact all-pairs clustering via connected components.
    Optionally writes a full pairwise distance table.
    """
    n = len(ids)

    if distance_table_path is not None:
        _write_empty_distance_table(distance_table_path)

    if n == 1:
        return np.array([0], dtype=int), {0: ids[0]}

    import multiprocessing as mp
    ctx = mp.get_context("fork")

    total_pairs  = n * (n - 1) // 2
    total_chunks = max(1, math.ceil(total_pairs / max(1, chunk_pairs)))
    done_chunks  = 0
    t0           = time.time()

    rows: List[int] = []
    cols: List[int] = []

    distance_fh = None
    try:
        if distance_table_path is not None:
            distance_fh = gzip.open(distance_table_path, "wt")
            distance_fh.write(
                "i\tj\tid_i\tid_j\tstatus\trmsd\t"
                "dihedral_alpha\tdihedral_beta\tclustering_dist\terror\n"
            )

        with ctx.Pool(
            processes=max(1, int(n_jobs)),
            initializer=_init_worker,
            initargs=(
                cluster_models,
                backbones,
                selection,
                anchors,
                cutoff,
                prune,
                distance_table_path is not None,
                distance_mode,
                w_rmsd,
                w_dihedral,
                rmsd_scale,
                cdr3_residues,
            ),
            maxtasksperchild=50,
        ) as pool:
            for distance_rows, edges in pool.imap_unordered(
                _process_pair_chunk,
                _pair_chunks(n, chunk_pairs),
                chunksize=1,
            ):
                done_chunks += 1

                if distance_fh is not None and distance_rows:
                    for row in distance_rows:
                        (ri, rj, status,
                         r_rmsd, r_dda, r_ddb, r_cd, error) = row

                        rmsd_str = "" if r_rmsd is None else f"{r_rmsd:.6f}"
                        dda_str  = "" if r_dda  is None else f"{r_dda:.6f}"
                        ddb_str  = "" if r_ddb  is None else f"{r_ddb:.6f}"
                        cd_str   = "" if r_cd   is None else f"{r_cd:.6f}"

                        distance_fh.write(
                            f"{ri}\t{rj}\t{ids[ri]}\t{ids[rj]}\t{status}\t"
                            f"{rmsd_str}\t{dda_str}\t{ddb_str}\t{cd_str}\t{error}\n"
                        )

                for ei, ej in edges:
                    rows.extend([ei, ej])
                    cols.extend([ej, ei])

                if progress_every_chunks > 0 and (
                    done_chunks % progress_every_chunks == 0
                    or done_chunks == total_chunks
                ):
                    elapsed = time.time() - t0
                    _log(
                        f"[INFO] exact: n={n:,} chunks {done_chunks}/{total_chunks} "
                        f"({done_chunks / total_chunks:.1%}) edges={len(rows) // 2:,} "
                        f"elapsed={elapsed:,.1f}s"
                    )
    finally:
        if distance_fh is not None:
            distance_fh.close()

    A = coo_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)),
        shape=(n, n),
    ).tocsr()
    A.setdiag(1)

    n_comp, labels = connected_components(A, directed=False, return_labels=True)

    rep: Dict[int, str] = {}
    for idx, lab in enumerate(labels):
        lab = int(lab)
        if lab not in rep:
            rep[lab] = ids[idx]

    _log(
        f"[INFO] exact: n={n:,} connected components={n_comp:,} "
        f"(cutoff={cutoff:.3f} mode={distance_mode})"
    )
    if distance_table_path is not None:
        _log(f"[INFO] wrote pairwise distance table: {distance_table_path}")

    return labels.astype(int), rep


# ============================================================
# Greedy representative clustering (large bins)
# ============================================================

def _greedy_representative_clustering(
    cluster_models: List,
    backbones: List,
    ids: List[str],
    selection: np.ndarray,
    anchors: np.ndarray,
    cutoff: float,
    distance_mode: str,
    w_rmsd: float,
    w_dihedral: float,
    rmsd_scale: float,
    cdr3_residues: List[int],
    progress_every: int = 5000,
) -> Tuple[np.ndarray, Dict[int, str]]:
    """
    Greedy representative clustering for bins too large for exact all-pairs.
    Respects the same distance mode as the exact path.
    """
    n = len(ids)
    reps_idx: List[int] = []
    labels = np.full(n, -1, dtype=int)

    # Each model's CDR3 dihedrals are extracted once, not once per
    # representative compared.
    dih_cache: Dict[int, Optional[Tuple[List[float], List[float]]]] = {}
    if distance_mode in ("dihedral", "composite"):
        for k in range(n):
            try:
                dih_cache[k] = _extract_cdr3_dihedrals(backbones[k], cdr3_residues)
            except Exception as exc:
                dih_cache[k] = None
                _log(f"[WARN] greedy: dihedral extraction failed for {ids[k]}: "
                     f"{type(exc).__name__}: {exc}")

    n_rmsd_fail = 0

    for i in range(n):
        assigned = False

        for c, r in enumerate(reps_idx):

            # RMSD component
            rmsd_val: Optional[float] = None
            if distance_mode in ("rmsd", "composite"):
                try:
                    rv = rmsd(
                        cluster_models[i], cluster_models[r],
                        selection=selection, anchors=anchors,
                    )
                    if np.isfinite(rv):
                        rmsd_val = float(rv)
                except Exception as exc:
                    # A failed RMSD is counted and reported, as the exact path
                    # records an "error" row.
                    n_rmsd_fail += 1
                    if n_rmsd_fail <= 20:
                        _log(f"[WARN] greedy: RMSD failed {ids[i]} vs {ids[r]}: "
                             f"{type(exc).__name__}: {exc}")

                if rmsd_val is None:
                    continue

            # Dihedral component
            dd_alpha: Optional[float] = None
            dd_beta:  Optional[float] = None
            if distance_mode in ("dihedral", "composite"):
                ci, cr = dih_cache.get(i), dih_cache.get(r)
                if ci is None or cr is None:
                    continue
                dd_alpha = _dihedral_vector_distance(ci[0], cr[0])
                dd_beta  = _dihedral_vector_distance(ci[1], cr[1])

            # Clustering distance
            if distance_mode == "rmsd":
                d = rmsd_val
            elif distance_mode == "dihedral":
                d = (dd_alpha + dd_beta) / 2.0   # radians; see _process_pair_chunk
            else:
                d = _composite_dist(rmsd_val, dd_alpha, dd_beta, w_rmsd, w_dihedral, rmsd_scale)

            if d <= cutoff:
                labels[i] = c
                assigned = True
                break

        if not assigned:
            reps_idx.append(i)
            labels[i] = len(reps_idx) - 1

        if progress_every > 0 and (i + 1) % progress_every == 0:
            _log(f"[WARN] greedy progress: {i + 1:,}/{n:,} reps={len(reps_idx):,}")

    if n_rmsd_fail:
        _log(f"[WARN] greedy: {n_rmsd_fail:,} RMSD computations failed and were "
             f"skipped; affected models may have been split into new clusters")

    rep = {c: ids[r] for c, r in enumerate(reps_idx)}
    return labels, rep


# ============================================================
# Bin loading
# ============================================================

def _load_bin_record(bins_jsonl: str, bin_index: int) -> dict:
    rec = None
    with open(bins_jsonl, "r") as fh:
        for line in fh:
            obj = json.loads(line)
            if int(obj["bin_index"]) == int(bin_index):
                rec = obj
                break
    if rec is None:
        raise SystemExit(f"[FATAL] bin_index {bin_index} not found in {bins_jsonl}")
    return rec


# ============================================================
# Main
# ============================================================

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Cluster one predefined length bin with optional dihedral angle support."
    )

    # ── Required ──────────────────────────────────────────────────────────────
    ap.add_argument("--bins-jsonl",  required=True)
    ap.add_argument("--bin-index",   type=int, required=True)
    ap.add_argument("--out-dir",     required=True)
    ap.add_argument("--cutoff",      type=float, required=True)
    ap.add_argument("--prune",       type=float, required=True)

    # ── PDB source ────────────────────────────────────────────────────────────
    ap.add_argument(
        "--pdb-base-dir",
        type=str,
        default=None,
        help=(
            "Filesystem mode: base dir joined with member paths from bins.jsonl. "
            "If set, tar fallback is not used."
        ),
    )

    # ── SPACE2 selection ──────────────────────────────────────────────────────
    ap.add_argument("--selection-key", default="CDR_all")
    ap.add_argument("--anchors-key",   default="fw_all")

    # ── Distance mode ─────────────────────────────────────────────────────────
    ap.add_argument(
        "--distance-mode",
        choices=["rmsd", "dihedral", "composite"],
        default="rmsd",
        help=(
            "rmsd      — original SPACE2 behaviour (default)\n"
            "dihedral  — cluster on CDR3 backbone dihedral distance only\n"
            "composite — weighted combination of normalised RMSD + dihedral"
        ),
    )
    ap.add_argument(
        "--w-rmsd",
        type=float,
        default=0.5,
        help="Weight on normalised RMSD in composite mode (default 0.5).",
    )
    ap.add_argument(
        "--w-dihedral",
        type=float,
        default=0.5,
        help="Weight on normalised dihedral distance in composite mode (default 0.5).",
    )
    ap.add_argument(
        "--rmsd-scale",
        type=float,
        default=3.0,
        help=(
            "RMSD value (Angstroms) mapped to 1.0 in composite normalisation. "
            "Should be approximately the maximum expected RMSD in your dataset "
            "(default 3.0)."
        ),
    )
    ap.add_argument(
        "--cdr3-imgt-start",
        type=int,
        default=105,
        help="First IMGT residue number of the CDR3 loop (default 105).",
    )
    ap.add_argument(
        "--cdr3-imgt-end",
        type=int,
        default=117,
        help="Last IMGT residue number of the CDR3 loop inclusive (default 117).",
    )

    # ── Performance ───────────────────────────────────────────────────────────
    ap.add_argument("--n-jobs",                 type=int, default=1)
    ap.add_argument("--max-exact-bin",          type=int, default=10000)
    ap.add_argument("--chunk-pairs",            type=int, default=25000)
    ap.add_argument("--progress-every-chunks",  type=int, default=25)
    ap.add_argument("--parse-progress-every",   type=int, default=1000)
    ap.add_argument("--greedy-seed",            type=int, default=42)

    # ── Distance table output ─────────────────────────────────────────────────
    ap.add_argument(
        "--write-distance-table",
        action="store_true",
        help=(
            "Write a gzipped TSV of all attempted pairwise comparisons. "
            "Requires exact mode (bin size <= max-exact-bin)."
        ),
    )
    ap.add_argument(
        "--distance-dir",
        default=None,
        help="Directory for per-bin pairwise tables. Default: <out-dir>/distance_tables",
    )

    args = ap.parse_args()

    # ── Validate weights ──────────────────────────────────────────────────────
    if args.distance_mode == "composite":
        if abs(args.w_rmsd + args.w_dihedral - 1.0) > 1e-6:
            raise SystemExit(
                f"[FATAL] --w-rmsd ({args.w_rmsd}) + --w-dihedral ({args.w_dihedral}) "
                f"must sum to 1.0"
            )

    cdr3_residues = list(range(args.cdr3_imgt_start, args.cdr3_imgt_end + 1))

    selection = _resolve_selection(args.selection_key)
    anchors   = _resolve_anchors(args.anchors_key)

    rec       = _load_bin_record(args.bins_jsonl, args.bin_index)
    tar_path  = rec.get("tar_path", None)
    members   = rec["members"]
    key       = tuple(rec["key"])
    n_declared = int(rec["n"])

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"bin_{args.bin_index:05d}.csv"

    distance_table_path: Optional[Path] = None
    if args.write_distance_table:
        distance_dir = (
            Path(args.distance_dir) if args.distance_dir
            else (out_dir / "distance_tables")
        )
        distance_dir.mkdir(parents=True, exist_ok=True)
        distance_table_path = distance_dir / f"bin_{args.bin_index:05d}.pairwise.tsv.gz"

    first_member = members[0] if members else "NA"
    _log(
        f"[INFO] bin={args.bin_index} size={n_declared:,} key={key} "
        f"mode={args.distance_mode} members[0]={first_member}"
    )
    if args.distance_mode != "rmsd":
        _log(
            f"[INFO] dihedral params: CDR3 IMGT {args.cdr3_imgt_start}-{args.cdr3_imgt_end} "
            f"w_rmsd={args.w_rmsd} w_dihedral={args.w_dihedral} rmsd_scale={args.rmsd_scale}"
        )

    # ── Parse members ─────────────────────────────────────────────────────────
    ids: List[str] = []
    models: List[Tuple[np.ndarray, np.ndarray]] = []
    backbones: List[Dict] = []
    bad: List[str] = []

    t0 = time.time()

    if args.pdb_base_dir:
        base  = Path(args.pdb_base_dir).expanduser().resolve()
        total = len(members)

        for i, mem in enumerate(members, start=1):
            pdb_path = base / mem
            pdb_id   = str(pdb_path)

            try:
                lines = pdb_path.read_text(errors="ignore").splitlines(keepends=True)
            except Exception:
                bad.append(pdb_id)
                continue

            model = _parse_pdb_lines_ca(lines)
            if model is None:
                bad.append(pdb_id)
                continue

            backbone = (
                _parse_pdb_backbone(lines)
                if args.distance_mode in ("dihedral", "composite")
                else {}
            )

            ids.append(pdb_id)
            models.append(model)
            backbones.append(backbone)

            if args.parse_progress_every > 0 and (
                i % args.parse_progress_every == 0 or i == total
            ):
                _log(
                    f"[INFO] parse: {i:,}/{total:,} "
                    f"kept={len(ids):,} bad={len(bad):,}"
                )

    else:
        if not tar_path:
            raise SystemExit(
                "[FATAL] No tar_path in bins.jsonl and --pdb-base-dir not provided."
            )

        total = len(members)
        with tarfile.open(tar_path, "r:*") as tf:
            for i, mem in enumerate(members, start=1):
                pdb_id = _tar_id(tar_path, mem)
                lines  = list(_iter_tar_member_lines(tf, mem))
                model  = _parse_pdb_lines_ca(lines)

                if model is None:
                    bad.append(pdb_id)
                    continue

                backbone = (
                    _parse_pdb_backbone(lines)
                    if args.distance_mode in ("dihedral", "composite")
                    else {}
                )

                ids.append(pdb_id)
                models.append(model)
                backbones.append(backbone)

                if args.parse_progress_every > 0 and (
                    i % args.parse_progress_every == 0 or i == total
                ):
                    _log(
                        f"[INFO] parse(tar): {i:,}/{total:,} "
                        f"kept={len(ids):,} bad={len(bad):,}"
                    )

    _log(
        f"[INFO] parsed bin: kept={len(ids):,} bad={len(bad):,} "
        f"elapsed={time.time() - t0:,.1f}s"
    )

    # ── Empty bin ─────────────────────────────────────────────────────────────
    if len(ids) == 0:
        pd.DataFrame(
            [],
            columns=[
                "ID", "cluster_by_length", "cluster_by_rmsd",
                "bin_index", "clustering_mode", "distance_mode",
            ],
        ).to_csv(out_csv, index=False)

        if distance_table_path is not None:
            _write_empty_distance_table(distance_table_path)
            _log(f"[INFO] wrote empty distance table: {distance_table_path}")

        _log(f"[INFO] wrote empty {out_csv}")
        return

    # ── Guard: distance table requires exact mode ──────────────────────────────
    if args.write_distance_table and len(ids) > int(args.max_exact_bin):
        raise SystemExit(
            f"[FATAL] --write-distance-table requested for bin {args.bin_index}, "
            f"but parsed bin size is {len(ids):,} and max_exact_bin={args.max_exact_bin:,}. "
            f"A full pairwise table requires exact all-pairs mode. "
            f"Increase --max-exact-bin or run this bin separately."
        )

    # ── Cluster ───────────────────────────────────────────────────────────────
    rng = np.random.default_rng(int(args.greedy_seed))

    if len(ids) <= int(args.max_exact_bin):
        labels, rep = _connected_components_threshold(
            cluster_models=models,
            backbones=backbones,
            ids=ids,
            selection=selection,
            anchors=anchors,
            cutoff=float(args.cutoff),
            prune=float(args.prune),
            n_jobs=int(args.n_jobs),
            chunk_pairs=max(1000, int(args.chunk_pairs)),
            progress_every_chunks=max(0, int(args.progress_every_chunks)),
            distance_mode=args.distance_mode,
            w_rmsd=args.w_rmsd,
            w_dihedral=args.w_dihedral,
            rmsd_scale=args.rmsd_scale,
            cdr3_residues=cdr3_residues,
            distance_table_path=distance_table_path,
        )
        mode_str = "exact"

    else:
        _log(
            f"[WARN] bin exceeds max_exact_bin={args.max_exact_bin:,}; "
            f"using greedy (shuffled)."
        )

        perm       = rng.permutation(len(ids))
        ids_shuf   = [ids[i]       for i in perm]
        models_shuf= [models[i]    for i in perm]
        backs_shuf = [backbones[i] for i in perm]

        labels_shuf, rep = _greedy_representative_clustering(
            cluster_models=models_shuf,
            backbones=backs_shuf,
            ids=ids_shuf,
            selection=selection,
            anchors=anchors,
            cutoff=float(args.cutoff),
            distance_mode=args.distance_mode,
            w_rmsd=args.w_rmsd,
            w_dihedral=args.w_dihedral,
            rmsd_scale=args.rmsd_scale,
            cdr3_residues=cdr3_residues,
            progress_every=5000,
        )

        labels = np.empty_like(labels_shuf)
        labels[perm] = labels_shuf
        mode_str = "greedy"

    # ── Write output CSV ───────────────────────────────────────────────────────
    key_str = "_".join(map(str, key))
    rows = [
        {
            "ID":               ids[i],
            "cluster_by_length":key_str,
            "cluster_by_rmsd":  rep[int(labels[i])],
            "bin_index":        int(args.bin_index),
            "clustering_mode":  mode_str,
            "distance_mode":    args.distance_mode,
        }
        for i, lab in enumerate(labels)
    ]

    pd.DataFrame(rows).to_csv(out_csv, index=False)
    _log(
        f"[INFO] wrote {out_csv} rows={len(rows):,} "
        f"mode={mode_str} distance={args.distance_mode}"
    )


if __name__ == "__main__":
    main()
