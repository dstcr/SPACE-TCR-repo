# scripts/dihedral_util.py

from __future__ import annotations
import numpy as np
from typing import List, Optional, Tuple

# ── Backbone atom parsing ──────────────────────────────────────────────────────

def parse_pdb_backbone(
    lines,
    chains=("A", "B"),
    max_resi=128
) -> dict:
    """
    Parse N, CA, C backbone atoms from PDB lines.
    Returns dict keyed by (chain, resi) -> {'N': array, 'CA': array, 'C': array}
    """
    atoms = {}
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

# ── Dihedral calculation ───────────────────────────────────────────────────────

def _dihedral(p0, p1, p2, p3) -> float:
    """
    Compute dihedral angle (radians, -pi to pi) defined by four points.
    Uses the standard cross-product / atan2 formula.
    """
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1_norm = b1 / (np.linalg.norm(b1) + 1e-10)

    # Gram-Schmidt: project b0 and b2 onto plane perpendicular to b1
    v = b0 - np.dot(b0, b1_norm) * b1_norm
    w = b2 - np.dot(b2, b1_norm) * b1_norm

    x = np.dot(v, w)
    y = np.dot(np.cross(b1_norm, v), w)

    return float(np.arctan2(y, x))


def extract_cdr3_dihedrals(
    backbone: dict,
    cdr3_residues_alpha: List[int],
    cdr3_residues_beta:  List[int],
) -> Tuple[List[float], List[float]]:
    """
    Extract phi/psi dihedral angles for CDR3alpha and CDR3beta residues.

    backbone: output of parse_pdb_backbone — dict keyed by (chain, resi)
    cdr3_residues_alpha: list of residue numbers (IMGT) for CDR3alpha
    cdr3_residues_beta:  list of residue numbers (IMGT) for CDR3beta

    Returns:
        dihedrals_alpha: flat list [phi_105, psi_105, phi_106, psi_106, ...]
        dihedrals_beta:  same for beta chain
    """
    def _chain_dihedrals(chain: str, residues: List[int]) -> List[float]:
        angles = []
        for i, resi in enumerate(residues):
            curr = backbone.get((chain, resi), {})
            prev = backbone.get((chain, residues[i-1]), {}) if i > 0 else {}
            nxt  = backbone.get((chain, residues[i+1]), {}) if i < len(residues)-1 else {}

            # phi = C(i-1) - N(i) - CA(i) - C(i)
            if all(k in prev for k in ("C",)) and \
               all(k in curr for k in ("N", "CA", "C")):
                phi = _dihedral(prev["C"], curr["N"], curr["CA"], curr["C"])
            else:
                phi = 0.0  # missing atoms: use 0 as neutral padding

            # psi = N(i) - CA(i) - C(i) - N(i+1)
            if all(k in curr for k in ("N", "CA", "C")) and \
               all(k in nxt  for k in ("N",)):
                psi = _dihedral(curr["N"], curr["CA"], curr["C"], nxt["N"])
            else:
                psi = 0.0

            angles.extend([phi, psi])
        return angles

    da = _chain_dihedrals("A", cdr3_residues_alpha)
    db = _chain_dihedrals("B", cdr3_residues_beta)
    return da, db

# ── Angular distance ───────────────────────────────────────────────────────────

def _circular_dist(a: float, b: float) -> float:
    """Angular distance between two angles in radians, result in [0, pi]."""
    diff = abs(a - b) % (2 * np.pi)
    return float(min(diff, 2 * np.pi - diff))


def dihedral_distance(
    angles_i: List[float],
    angles_j: List[float],
) -> float:
    """
    Mean circular distance between two dihedral angle vectors.
    Handles length mismatch by padding shorter vector with pi (max distance).
    Returns value in [0, pi].
    """
    if not angles_i and not angles_j:
        return 0.0
    if not angles_i or not angles_j:
        return np.pi  # one loop completely missing

    max_len = max(len(angles_i), len(angles_j))

    # Pad shorter with pi — signals maximum disagreement for missing residues
    ai = list(angles_i) + [np.pi] * (max_len - len(angles_i))
    aj = list(angles_j) + [np.pi] * (max_len - len(angles_j))

    dists = [_circular_dist(a, b) for a, b in zip(ai, aj)]
    return float(np.mean(dists))


def composite_distance(
    rmsd_val: float,
    dihedral_dist_alpha: float,
    dihedral_dist_beta: float,
    w_rmsd: float = 0.5,
    w_dihedral: float = 0.5,
    rmsd_scale: float = 3.0,      # expected max RMSD in Angstroms
    dihedral_scale: float = np.pi, # max possible dihedral distance
) -> float:
    """
    Weighted composite of normalised RMSD and mean dihedral distance.
    Both components normalised to [0, 1] before weighting.
    """
    rmsd_norm     = min(rmsd_val / rmsd_scale, 1.0)
    dihedral_mean = (dihedral_dist_alpha + dihedral_dist_beta) / 2.0
    dihedral_norm = dihedral_mean / dihedral_scale

    return w_rmsd * rmsd_norm + w_dihedral * dihedral_norm
