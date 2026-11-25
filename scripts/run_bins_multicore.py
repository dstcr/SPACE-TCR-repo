#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple


def load_bin_indices(bins_jsonl: str, start: int = 0, stop: Optional[int] = None) -> List[int]:
    idx: List[int] = []
    with open(bins_jsonl, "r") as fh:
        for line in fh:
            rec = json.loads(line)
            bin_index = int(rec["bin_index"])
            if bin_index < start:
                continue
            if stop is not None and bin_index >= stop:
                continue
            idx.append(bin_index)
    return idx


def run_one(
    py: str,
    cluster_one_bin: str,
    bins_jsonl: str,
    pdb_base_dir: str,
    out_dir: str,
    cutoff: float,
    prune: float,
    selection_key: str,
    anchors_key: str,
    n_jobs_per_bin: int,
    max_exact_bin: int,
    chunk_pairs: int,
    progress_every_chunks: int,
    parse_progress_every: int,
    greedy_seed: int,
    write_distance_table: bool,
    distance_dir: Optional[str],
    # ── new distance-mode arguments ──────────────────────────
    distance_mode: str,
    w_rmsd: float,
    w_dihedral: float,
    rmsd_scale: float,
    cdr3_imgt_start: int,
    cdr3_imgt_end: int,
    # ─────────────────────────────────────────────────────────
    bin_index: int,
) -> Tuple[int, int]:
    cmd = [
        py,
        "-u",
        cluster_one_bin,
        "--bins-jsonl",      bins_jsonl,
        "--bin-index",       str(bin_index),
        "--pdb-base-dir",    pdb_base_dir,
        "--out-dir",         out_dir,
        "--cutoff",          str(cutoff),
        "--prune",           str(prune),
        "--selection-key",   selection_key,
        "--anchors-key",     anchors_key,
        "--n-jobs",          str(n_jobs_per_bin),
        "--max-exact-bin",   str(max_exact_bin),
        "--chunk-pairs",     str(chunk_pairs),
        "--progress-every-chunks", str(progress_every_chunks),
        "--parse-progress-every",  str(parse_progress_every),
        "--greedy-seed",     str(greedy_seed),
        # ── distance mode ────────────────────────────────────
        "--distance-mode",   distance_mode,
        "--w-rmsd",          str(w_rmsd),
        "--w-dihedral",      str(w_dihedral),
        "--rmsd-scale",      str(rmsd_scale),
        "--cdr3-imgt-start", str(cdr3_imgt_start),
        "--cdr3-imgt-end",   str(cdr3_imgt_end),
    ]

    if write_distance_table:
        cmd.append("--write-distance-table")
        if distance_dir:
            cmd.extend(["--distance-dir", distance_dir])

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        fail_log = Path(out_dir) / f"bin_{bin_index:05d}.FAIL.log"
        fail_log.write_text(
            "CMD:\n"
            + " ".join(cmd)
            + "\n\nSTDOUT:\n"
            + proc.stdout
            + "\n\nSTDERR:\n"
            + proc.stderr
            + "\n"
        )

    return bin_index, proc.returncode


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run predefined-bin SPACE2 clustering in parallel within one job."
    )

    # ── Required ──────────────────────────────────────────────────────────────
    ap.add_argument("--bins-jsonl",       required=True)
    ap.add_argument("--pdb-base-dir",     required=True)
    ap.add_argument("--out-dir",          required=True)
    ap.add_argument("--cluster-one-bin",  required=True)
    ap.add_argument("--python",           default=sys.executable)
    ap.add_argument("--cutoff",           type=float, required=True)
    ap.add_argument("--prune",            type=float, required=True)

    # ── SPACE2 selection ──────────────────────────────────────────────────────
    ap.add_argument("--selection-key",    default="CDR_all")
    ap.add_argument("--anchors-key",      default="fw_all")

    # ── Distance mode (forwarded to cluster_one_bin.py) ───────────────────────
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
        help="RMSD value (Å) mapped to 1.0 in composite normalisation (default 3.0).",
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
    ap.add_argument("--n-jobs-per-bin",          type=int, default=1)
    ap.add_argument("--max-exact-bin",           type=int, default=5000)
    ap.add_argument("--chunk-pairs",             type=int, default=100000)
    ap.add_argument("--progress-every-chunks",   type=int, default=25)
    ap.add_argument("--parse-progress-every",    type=int, default=2000)
    ap.add_argument("--greedy-seed",             type=int, default=42)
    ap.add_argument("--start",                   type=int, default=0)
    ap.add_argument("--stop",                    type=int, default=None)
    ap.add_argument("--workers",                 type=int, default=None)

    # ── Distance table ────────────────────────────────────────────────────────
    ap.add_argument("--write-distance-table",    action="store_true")
    ap.add_argument("--distance-dir",            default=None)

    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    slurm_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    workers = args.workers if args.workers is not None else min(16, slurm_cpus)

    if args.distance_dir is None and args.write_distance_table:
        distance_dir = str(out_dir / "distance_tables")
    else:
        distance_dir = args.distance_dir

    bins = load_bin_indices(args.bins_jsonl, start=args.start, stop=args.stop)

    print(f"[INFO] bins loaded: {len(bins)} (range {args.start}..{args.stop})", flush=True)
    print(
        f"[INFO] workers={workers} slurm_cpus={slurm_cpus} "
        f"n_jobs_per_bin={args.n_jobs_per_bin} max_exact_bin={args.max_exact_bin}",
        flush=True,
    )
    print(
        f"[INFO] distance_mode={args.distance_mode} "
        f"w_rmsd={args.w_rmsd} w_dihedral={args.w_dihedral} "
        f"rmsd_scale={args.rmsd_scale}",
        flush=True,
    )
    if args.write_distance_table:
        print(f"[INFO] write-distance-table=ON distance-dir={distance_dir}", flush=True)

    ok = 0
    fail = 0

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [
            ex.submit(
                run_one,
                args.python,
                args.cluster_one_bin,
                args.bins_jsonl,
                args.pdb_base_dir,
                str(out_dir),
                args.cutoff,
                args.prune,
                args.selection_key,
                args.anchors_key,
                args.n_jobs_per_bin,
                args.max_exact_bin,
                args.chunk_pairs,
                args.progress_every_chunks,
                args.parse_progress_every,
                args.greedy_seed,
                args.write_distance_table,
                distance_dir,
                # ── new args ──────────────────────────────────
                args.distance_mode,
                args.w_rmsd,
                args.w_dihedral,
                args.rmsd_scale,
                args.cdr3_imgt_start,
                args.cdr3_imgt_end,
                # ─────────────────────────────────────────────
                bin_index,
            )
            for bin_index in bins
        ]

        for fut in as_completed(futures):
            _bin_index, rc = fut.result()
            if rc == 0:
                ok += 1
            else:
                fail += 1

            done = ok + fail
            if done % 50 == 0 or done == len(bins):
                print(f"[INFO] done {done}/{len(bins)} ok={ok} fail={fail}", flush=True)

    print(f"[INFO] finished: ok={ok} fail={fail}", flush=True)


if __name__ == "__main__":
    main()
