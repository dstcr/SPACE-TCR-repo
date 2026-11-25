#!/usr/bin/env python3
"""
make_bins.py — Stage3-equivalent length binning, writing bins.jsonl + member_to_bin.tsv

This version supports BOTH:
1. filesystem mode via --pdb-dir
2. tar mode via --pdb-tar

Recommended for your reuse pipeline:
    --pdb-dir /path/to/pdb/TCRPredDab

Outputs:
- bins.jsonl
- member_to_bin.tsv

For filesystem mode, "members" in bins.jsonl are stored relative to --pdb-dir
(e.g. tcr1.pdb, tcr2.pdb, ...).
"""

from __future__ import annotations

import argparse
import json
import tarfile
import time
from pathlib import Path
from typing import Iterable, List, Optional, Set, Tuple, Union

import numpy as np

from SPACE2TCR.util import (
    cluster_antibodies_by_CDR_length,
    reg_def,
    same_length_only,
)


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _log(msg: str) -> None:
    print(f"[{_now()}] {msg}", flush=True)


def _tar_id(tar_path: str, member_name: str) -> str:
    return f"tar:{str(Path(tar_path).expanduser().resolve())}::{member_name}"


def _fs_id(base_dir: Path, file_path: Path) -> str:
    return str(file_path.resolve())


def _read_bad_members(path: str) -> Set[str]:
    """
    Accept either:
      - member/path/in/tar.pdb
      - tar:/abs/file.tar::member/path/in/tar.pdb
      - bare filename for filesystem mode
      - full path for filesystem mode

    Store strings as-is after normalisation where possible.
    """
    p = Path(path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"bad-pdbs-in not found: {path}")
    bad: Set[str] = set()
    for line in p.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("tar:") and "::" in s:
            bad.add(s.split("::", 1)[1])
        else:
            bad.add(s)
    return bad


def _sel_for_length_padded(sel_key: str) -> List[np.ndarray]:
    """
    Same intent as stage3_space2_tcr3.py: 6-part selection list [B1,B2,B3,A1,A2,A3].
    """
    empty3 = [np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)]

    if sel_key == "CDR_A_all":
        return empty3 + [np.asarray(x).ravel().astype(int) for x in reg_def["CDR_A_all"]]

    if sel_key == "CDR_B_all":
        return [np.asarray(x).ravel().astype(int) for x in reg_def["CDR_B_all"]] + empty3

    both = reg_def.get("CDR_all", reg_def["CDR_B_all"] + reg_def["CDR_A_all"])
    return [np.asarray(x).ravel().astype(int) for x in both]


def _parse_pdb_ca_numbers_coords_bytes(raw_iter: Iterable[bytes]) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Fast CA-only parser (Stage3-compatible output):
    Returns (numbers, coords) where
      numbers: int indices resi (B) or resi+128 (A), resi in 1..128
      coords: float32 array (n,3)
    """
    numbers: List[int] = []
    coords: List[Tuple[float, float, float]] = []

    for raw in raw_iter:
        if len(raw) < 54:
            continue
        if raw[0:4] != b"ATOM":
            continue
        if raw[12:16].strip() != b"CA":
            continue

        chain = raw[21:22]
        if chain not in (b"A", b"B"):
            continue

        try:
            resi = int(raw[22:26])
            if resi <= 0 or resi > 128:
                continue
            x = float(raw[30:38])
            y = float(raw[38:46])
            z = float(raw[46:54])
        except Exception:
            try:
                parts = raw.decode("utf-8", errors="ignore").split()
                chain2 = parts[4]
                if chain2 not in ("A", "B"):
                    continue
                resi = int(parts[5])
                if resi <= 0 or resi > 128:
                    continue
                x = float(parts[6])
                y = float(parts[7])
                z = float(parts[8])
                chain = chain2.encode()
            except Exception:
                return None

        idx = resi + (128 if chain == b"A" else 0)
        numbers.append(idx)
        coords.append((x, y, z))

    if not numbers:
        return None

    return np.asarray(numbers, dtype=int), np.asarray(coords, dtype=np.float32)


def _key_to_int_list(k: Union[str, Tuple, List, np.ndarray]) -> List[int]:
    if isinstance(k, str):
        return [int(x) for x in k.split("_")]
    return [int(x) for x in list(k)]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build SPACE2 bins using Stage3 binning logic from either filesystem PDBs or a tar archive."
    )

    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdb-dir", default=None, help="Directory containing PDB files.")
    src.add_argument("--pdb-tar", default=None, help="Tar archive containing PDB files.")

    ap.add_argument("--tar-member-prefix", default=None)
    ap.add_argument("--tar-member-suffix", default=".pdb")

    ap.add_argument("--selection-key", default="CDR_all", choices=list(reg_def.keys()))
    ap.add_argument("--bad-pdbs-in", default=None)
    ap.add_argument("--out-manifest", required=True)
    ap.add_argument("--out-index", required=True)
    ap.add_argument("--progress-every", type=int, default=50000)

    args = ap.parse_args()

    prefix = args.tar_member_prefix
    suffix = args.tar_member_suffix

    _log("[INFO] make_bins starting")
    if args.pdb_dir:
        _log(f"[INFO] source mode = filesystem")
        _log(f"[INFO] pdb_dir={args.pdb_dir}")
    else:
        _log(f"[INFO] source mode = tar")
        _log(f"[INFO] pdb_tar={args.pdb_tar}")
        _log(f"[INFO] prefix={prefix!r} suffix={suffix!r}")

    _log(f"[INFO] selection-key={args.selection_key}")

    bad_members: Set[str] = set()
    if args.bad_pdbs_in:
        bad_members = _read_bad_members(args.bad_pdbs_in)
        _log(f"[INFO] loaded bad members: n={len(bad_members):,}")

    antibodies: List[Tuple[np.ndarray, np.ndarray]] = []
    ids: List[str] = []
    members_for_manifest: List[str] = []

    seen = 0
    kept = 0
    malformed = 0
    skipped_bad = 0
    skipped_filter = 0
    t0 = time.time()

    # ------------------------------------------------------------------
    # Filesystem mode
    # ------------------------------------------------------------------
    if args.pdb_dir:
        pdb_dir = Path(args.pdb_dir).expanduser().resolve()
        if not pdb_dir.exists():
            raise SystemExit(f"[FATAL] pdb-dir does not exist: {pdb_dir}")

        files = sorted(pdb_dir.glob("*.pdb"))
        if not files:
            raise SystemExit(f"[FATAL] no *.pdb files found in: {pdb_dir}")

        total = len(files)

        for i, fp in enumerate(files, start=1):
            seen += 1
            member_name = fp.name

            if member_name in bad_members or str(fp) in bad_members:
                skipped_bad += 1
                continue

            try:
                with open(fp, "rb") as f:
                    model = _parse_pdb_ca_numbers_coords_bytes(f)
            except Exception:
                model = None

            if model is None:
                malformed += 1
                continue

            antibodies.append(model)
            ids.append(_fs_id(pdb_dir, fp))
            members_for_manifest.append(member_name)
            kept += 1

            if args.progress_every and (kept % int(args.progress_every) == 0 or i == total):
                elapsed = time.time() - t0
                rate = kept / elapsed if elapsed > 0 else float("nan")
                _log(
                    f"[INFO] fs parse progress: kept={kept:,} seen={seen:,} rate={rate:,.1f}/s "
                    f"malformed={malformed:,} skipped_bad={skipped_bad:,}"
                )

    # ------------------------------------------------------------------
    # Tar mode
    # ------------------------------------------------------------------
    else:
        tar_path = str(Path(args.pdb_tar).expanduser().resolve())

        with tarfile.open(tar_path, "r:*") as tf:
            for m in tf:
                if not m.isfile():
                    continue

                name = m.name
                seen += 1

                if prefix is not None and not name.startswith(prefix):
                    skipped_filter += 1
                    continue
                if suffix and not name.endswith(suffix):
                    skipped_filter += 1
                    continue
                if name in bad_members:
                    skipped_bad += 1
                    continue

                f = tf.extractfile(m)
                if f is None:
                    malformed += 1
                    continue

                model = _parse_pdb_ca_numbers_coords_bytes(f)
                if model is None:
                    malformed += 1
                    continue

                antibodies.append(model)
                ids.append(_tar_id(tar_path, name))
                members_for_manifest.append(name)
                kept += 1

                if args.progress_every and (kept % int(args.progress_every) == 0):
                    elapsed = time.time() - t0
                    rate = kept / elapsed if elapsed > 0 else float("nan")
                    _log(
                        f"[INFO] tar parse progress: kept={kept:,} seen={seen:,} rate={rate:,.1f}/s "
                        f"malformed={malformed:,} skipped_bad={skipped_bad:,}"
                    )

    _log(
        f"[INFO] parsing done: kept={kept:,} malformed={malformed:,} "
        f"skipped_bad={skipped_bad:,} elapsed={time.time()-t0:,.1f}s"
    )

    if kept < 2:
        raise SystemExit("[FATAL] Too few models kept to bin.")

    sel_for_len = _sel_for_length_padded(args.selection_key)

    _log("[INFO] running cluster_antibodies_by_CDR_length")
    cdr_clusters, cdr_ids = cluster_antibodies_by_CDR_length(
        antibodies,
        ids,
        selection=sel_for_len,
        clustering="bins",
        tolerance=same_length_only,
    )

    keys = sorted(cdr_ids, key=lambda kk: len(cdr_ids[kk]), reverse=True)
    _log(f"[INFO] Length bins: {len(keys):,}")
    if keys:
        k0 = _key_to_int_list(keys[0])
        _log(f"[INFO] Largest bin size: {len(cdr_ids[keys[0]]):,} key={'_'.join(map(str, k0))}")

    out_manifest = Path(args.out_manifest)
    out_index = Path(args.out_index)
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_index.parent.mkdir(parents=True, exist_ok=True)

    tmp_manifest = out_manifest.with_suffix(out_manifest.suffix + ".tmp")
    tmp_index = out_index.with_suffix(out_index.suffix + ".tmp")

    # Map id -> member name stored in manifest/index
    id_to_member = dict(zip(ids, members_for_manifest))

    _log(f"[INFO] writing manifest: {out_manifest}")
    _log(f"[INFO] writing index:    {out_index}")

    with tmp_manifest.open("w") as mf, tmp_index.open("w") as ix:
        ix.write("member\tbin_index\tkey\n")

        for bin_i, k in enumerate(keys):
            bin_ids = cdr_ids[k]
            members = [id_to_member[tid] for tid in bin_ids]

            k_ints = _key_to_int_list(k)
            key_str = "_".join(map(str, k_ints))

            rec = {
                "bin_index": int(bin_i),
                "key": k_ints,
                "n": int(len(members)),
                "members": members,
            }
            mf.write(json.dumps(rec) + "\n")

            for mem in members:
                ix.write(f"{mem}\t{bin_i}\t{key_str}\n")

    tmp_manifest.replace(out_manifest)
    tmp_index.replace(out_index)

    _log("[INFO] make_bins finished OK")


if __name__ == "__main__":
    main()
