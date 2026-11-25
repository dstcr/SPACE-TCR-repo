#!/usr/bin/env python3
"""
Stage 1 — CSV → (optional) Thimble/STITCHr → FASTA for TCRBuilder2

What it does:
- Reads every .csv or .csv.gz in --in-dir (delimiter/BOM tolerant).
- Normalises headers case-insensitively and accepts many aliases.
- If full variable-domain AA are present (alpha/beta), uses them directly.
- Otherwise reconstructs V-domain AA via Thimble/STITCHr from TRAV/TRAJ/CDR3a and TRBV/TRBJ/CDR3b.
- Writes one FASTA per clonotype (order: >B then >A) with resume-safe behaviour:
  * does not overwrite; duplicates get _d1, _d2, … suffixes, consistent across reruns.
  * skips any FASTA already on disk.
- Emits a per-file summary and a global summary CSV with counts and reasons.

CLI example
-----------
python stage1_stitchr.py \
  --in-dir       /path/input_csv \
  --stitchr-dir  /path/stitchr \
  --fasta-aa-dir /path/fasta_aa \
  --work         /path/stitchr/work \
  --species      human \
  --reuse-thimble
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


# ─────────────────────────────── Utilities ─────────────────────────────── #

def read_csv_any(path: Path) -> pd.DataFrame:
    """Read CSV/CSV.GZ with delimiter/BOM tolerance; force string dtype."""
    if path.suffix == ".gz":
        with gzip.open(path, mode="rt", encoding="utf-8", errors="replace") as fh:
            return pd.read_csv(fh, sep=None, engine="python", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str, encoding_errors="replace")


def normalise_headers(df: pd.DataFrame) -> pd.DataFrame:
    """Lower-case, strip spaces/quotes from headers."""
    # Drop surrounding quotes some exporters add.
    def _clean(c: str) -> str:
        c = c.strip().strip('"').strip("'")
        c = re.sub(r"\s+", "", c)
        return c.lower()
    return df.rename(columns={c: _clean(c) for c in df.columns})


def first_present(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[pd.Series]:
    """Return the first present column (case-insensitive); else None."""
    lc = {c.lower(): c for c in df.columns}
    for cand in candidates:
        cc = cand.lower()
        if cc in lc:
            return df[lc[cc]]
    return None


def series_norm(s: Optional[pd.Series], n: int, name: str) -> pd.Series:
    """Coerce to string series; NA/""/"NA"/"N/A" → NA."""
    if s is None:
        return pd.Series([pd.NA] * n, name=name, dtype="string")
    s = s.astype("string")
    return s.mask(s.isna() | s.str.strip().eq("") | s.str.upper().isin({"NA", "N/A"}))


AA20 = set("ACDEFGHIKLMNPQRSTVWY")  # we reject non-standard letters for modelling


def is_aa(seq: str) -> bool:
    return all(ch in AA20 for ch in seq)


def sane_cdr3(seq: str) -> bool:
    seq = seq.replace("*", "")
    return (5 <= len(seq) <= 30) and is_aa(seq)


def next_unique_path(base: Path) -> Path:
    """Return a unique path: base or base_dN if base exists (stable across reruns)."""
    if not base.exists():
        return base
    stem, ext = base.stem, base.suffix
    i = 1
    while True:
        p = base.with_name(f"{stem}_d{i}{ext}")
        if not p.exists():
            return p
        i += 1


def init_dup_counters_from_disk(fasta_dir: Path) -> Dict[str, int]:
    """Scan existing FASTAs to initialise duplicate counters per base."""
    counters: Dict[str, int] = defaultdict(int)
    pat = re.compile(r"^(?P<base>.+?)(?:_d(?P<d>\d+))?$")
    for p in fasta_dir.glob("*.fasta"):
        m = pat.match(p.stem)
        if not m:
            continue
        base = m.group("base")
        d = int(m.group("d")) if m.group("d") else 0
        counters[base] = max(counters[base], d)
    return counters


# ────────────────────── Header/alias handling ────────────────────── #

ALIASES = {
    # identifiers
    "id":      ("barcode", "cell_id", "contig_id", "clone_id", "row", "row_id", "id", "name"),
    # alpha chain
    "trav":    ("trav", "v_call_alpha", "v_call_a", "tra_v", "trav_gene", "tcrv_alpha"),
    "traj":    ("traj", "j_call_alpha", "j_call_a", "tra_j", "traj_gene", "tcrj_alpha"),
    "cdr3a":   ("cdr3alpha", "cdr3a", "cdr3_aa_alpha", "tra_cdr3", "tra_cdr3aa", "cdr3a_aa"),
    # beta chain
    "trbv":    ("trbv", "v_call_beta", "v_call_b", "trb_v", "trbv_gene", "tcrv_beta"),
    "trbj":    ("trbj", "j_call_beta", "j_call_b", "trb_j", "trbj_gene", "tcrj_beta"),
    "cdr3b":   ("cdr3beta", "cdr3b", "cdr3_aa_beta", "trb_cdr3", "trb_cdr3aa", "cdr3b_aa"),
    # full V-domain AA (if present)
    "alpha_aa": ("sequence_alignment_aa_alpha", "alpha_aa", "aa_alpha", "vdomain_alpha", "v_alpha_seq", "v_aa_alpha"),
    "beta_aa":  ("sequence_alignment_aa_beta",  "beta_aa",  "aa_beta",  "vdomain_beta",  "v_beta_seq",  "v_aa_beta"),
}


def coerce_table(df_raw: pd.DataFrame, fname: str) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Return a canonical table with columns:
      row_id, trav, traj, cdr3a, trbv, trbj, cdr3b, alpha_aa, beta_aa
    and a dict of drop reasons.
    """
    df = normalise_headers(df_raw)
    n = len(df)
    drops = Counter()

    def pull(key: str) -> pd.Series:
        return series_norm(first_present(df, ALIASES[key]), n, key)

    row_id = pull("id")
    if row_id.isna().any():
        # deterministic fallback if missing
        row_id = row_id.fillna(pd.Series([f"{Path(fname).stem}_{i:06d}" for i in range(n)], dtype="string"))

    out = pd.DataFrame({
        "row_id":  row_id,
        "trav":    pull("trav"),
        "traj":    pull("traj"),
        "cdr3a":   pull("cdr3a"),
        "trbv":    pull("trbv"),
        "trbj":    pull("trbj"),
        "cdr3b":   pull("cdr3b"),
        "alpha_aa": series_norm(first_present(df, ALIASES["alpha_aa"]), n, "alpha_aa"),
        "beta_aa":  series_norm(first_present(df, ALIASES["beta_aa"]),  n, "beta_aa"),
    })

    # Clean CDR3s and AA; record reasons for drops
    keep_mask = []
    for _, r in out.iterrows():
        ok = True
        # V/J must exist
        for nm in ("trav", "traj", "trbv", "trbj"):
            if pd.isna(r[nm]):
                drops[f"missing_{nm}"] += 1
                ok = False
                break
        if not ok:
            keep_mask.append(False)
            continue
        # CDR3s must be sane
        a = (r["cdr3a"] or "").replace("*", "")
        b = (r["cdr3b"] or "").replace("*", "")
        if not (sane_cdr3(a) and sane_cdr3(b)):
            drops["bad_cdr3"] += 1
            keep_mask.append(False)
            continue
        keep_mask.append(True)

    kept = out[pd.Series(keep_mask, index=out.index)].reset_index(drop=True)
    return kept, dict(drops)


# ─────────────────────────── Thimble I/O ─────────────────────────── #

def write_thimble_inputs(df: pd.DataFrame, tra_tsv: Path, trb_tsv: Path) -> None:
    """
    Build chain-specific Thimble TSVs with required headers.
    Both files are written from the same table; Thimble selects by -r a|b.
    """
    header = [
        "TCR_name", "TRAV", "TRAJ", "TRA_CDR3",
        "TRBV", "TRBJ", "TRB_CDR3",
        "TRAC", "TRBC", "TRA_leader", "TRB_leader",
        "Linker", "Link_order",
        "TRA_5_prime_seq", "TRA_3_prime_seq", "TRB_5_prime_seq", "TRB_3_prime_seq",
    ]
    tra_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(tra_tsv, "w", newline="") as fa, open(trb_tsv, "w", newline="") as fb:
        wA, wB = csv.writer(fa, delimiter="\t"), csv.writer(fb, delimiter="\t")
        wA.writerow(header); wB.writerow(header)
        for _, r in df.iterrows():
            rid = str(r["row_id"])
            wA.writerow([rid, r["trav"], r["traj"], r["cdr3a"], "", "", "", "", "", "", "", "", "", "", "", "", ""])
            wB.writerow([rid, "", "", "", r["trbv"], r["trbj"], r["cdr3b"], "", "", "", "", "", "", "", "", "", ""])


def run_thimble(in_path: Path, out_path: Path, receptor_flag: str, species: str) -> None:
    """
    Call thimble with the correct CLI used by current builds.
    -in / -o / -r a|b / -s human|mouse
    """
    cmd = [
        "thimble", "-in", str(in_path), "-o", str(out_path),
        "-r", receptor_flag, "-s", species.lower()
    ]
    subprocess.run(cmd, check=True)


def parse_thimble_output(path: Path, chain: str = 'alpha') -> pd.DataFrame:
    """Return a tidy DataFrame with columns: row_id, aa (variable-domain)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = normalise_headers(df)
    # Identify name and AA columns across builds
    name = first_present(df, ("tcr_name", "row_id", "id", "name"))
    if chain == "beta":
        aa = first_present(df, ("trb_aa", "aa", "aaseq", "proteinseq"))
    else:
        aa = first_present(df, ("tra_aa", "aa", "aaseq", "proteinseq"))
    if name is None or aa is None:
        raise ValueError(f"Thimble output missing required columns: {path}")
    out = pd.DataFrame({
        "row_id": name.astype("string").str.strip(),
        "aa": aa.astype("string").str.replace("*", "", regex=False).str.strip()
    })
    out = out.dropna()
    # Keep only plausible AA strings
    out = out[out["aa"].map(lambda s: isinstance(s, str) and len(s) >= 70 and len(s) <= 600 and is_aa(s))]
    return out.reset_index(drop=True)


# ───────────────────────── FASTA emission ───────────────────────── #

def write_fastas(df: pd.DataFrame, fasta_dir: Path,
                 alpha_map: Dict[str, str], beta_map: Dict[str, str],
                 summary_counts: Counter) -> None:
    """
    For each row: pick β/α AA (prefer AA maps; no maps → skip), then write FASTA:
       >B\n<beta>\n>A\n<alpha>\n
    Resume-safe; stable suffixing.
    """
    fasta_dir.mkdir(parents=True, exist_ok=True)
    dup = init_dup_counters_from_disk(fasta_dir)

    for _, r in df.iterrows():
        rid = str(r["row_id"])
        a = r.get("alpha_aa")
        b = r.get("beta_aa")
        # Prefer thimble maps when present; else direct full AA if available.
        a_seq = (alpha_map.get(rid) or (a if pd.notna(a) and is_aa(a) else None))
        b_seq = (beta_map.get(rid)  or (b if pd.notna(b) and is_aa(b) else None))
        if not (a_seq and b_seq):
            summary_counts["skip_missing_aa"] += 1
            continue

        # Base name and suffix management
        base = rid
        if base not in dup:
            dup[base] = 0
            dest = fasta_dir / f"{base}.fasta"
            if dest.exists():  # already present → advance suffix
                dup[base] = 1
                dest = fasta_dir / f"{base}_d1.fasta"
        else:
            dup[base] += 1
            dest = fasta_dir / f"{base}_d{dup[base]}.fasta"

        if dest.exists():
            continue  # we do not overwrite
        with open(dest, "w") as fh:
            fh.write(">B\n"); fh.write(f"{b_seq}\n")
            fh.write(">A\n"); fh.write(f"{a_seq}\n")
        summary_counts["written_fastas"] += 1


# ─────────────────────────────── Main ─────────────────────────────── #

def main() -> None:
    ap = argparse.ArgumentParser(description="CSV → Thimble → FASTA(>B,>A) for TCRBuilder2 (robust).")
    ap.add_argument("--in-dir",        required=True, help="Directory with input .csv / .csv.gz files.")
    ap.add_argument("--stitchr-dir",   required=True, help="Directory for Thimble/STITCHr outputs (TSV).")
    ap.add_argument("--fasta-aa-dir",  required=True, help="Output directory for FASTA files.")
    ap.add_argument("--work",          default="work", help="Working directory (TSVs/logs).")
    ap.add_argument("--species",       default="human", choices=["human", "mouse"], help="Species for Thimble.")
    ap.add_argument("--reuse-thimble", action="store_true", help="Reuse existing Thimble outputs if present.")
    args = ap.parse_args()

    in_dir      = Path(args.in_dir)
    stitchr_dir = Path(args.stitchr_dir); stitchr_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir   = Path(args.fasta_aa_dir); fasta_dir.mkdir(parents=True, exist_ok=True)
    work_dir    = Path(args.work); (work_dir / "logs").mkdir(parents=True, exist_ok=True)

    # Discover CSVs
    csvs: List[Path] = []
    for p in sorted(in_dir.iterdir()):
        if p.is_file() and (p.suffix == ".csv" or p.suffixes[-2:] == [".csv", ".gz"]):
            csvs.append(p)
    if not csvs:
        sys.exit(f"[stage1] No CSV/CSV.GZ files found in {in_dir}")

    global_summary_rows: List[dict] = []
    kept_frames: List[pd.DataFrame] = []

    # Parse and coerce each file independently, to report per-file stats.
    for csv_path in csvs:
        summary = {"file": csv_path.name}
        try:
            raw = read_csv_any(csv_path)
        except Exception as e:
            summary.update({"n_rows": 0, "n_kept": 0, "note": f"read_error:{e}"})
            global_summary_rows.append(summary); continue

        n0 = len(raw)
        table, drop_reasons = coerce_table(raw, csv_path.name)
        nkept = len(table)
        summary.update({"n_rows": n0, "n_kept": nkept})
        summary.update({f"drop_{k}": v for k, v in drop_reasons.items()})
        global_summary_rows.append(summary)

        if not table.empty:
            kept_frames.append(table)

    if not kept_frames:
        # Write a summary even on failure, recording the reason.
        pd.DataFrame(global_summary_rows).to_csv(work_dir / "stage1_summary.csv", index=False)
        sys.exit("[stage1] No usable clonotypes after coercion. See summary for reasons.")

    df_all = pd.concat(kept_frames, ignore_index=True)

    # Separate rows that already provide full V-domain AA
    have_alpha = df_all["alpha_aa"].notna() & df_all["alpha_aa"].map(lambda s: isinstance(s, str) and is_aa(s))
    have_beta  = df_all["beta_aa"].notna()  & df_all["beta_aa"].map(lambda s: isinstance(s, str) and is_aa(s))
    df_direct  = df_all[have_alpha & have_beta].copy()
    df_need    = df_all[~(have_alpha & have_beta)].copy()

    # Build/Run Thimble if needed
    tra_in = work_dir / "thimble_TRA_input.tsv"
    trb_in = work_dir / "thimble_TRB_input.tsv"
    tra_out = stitchr_dir / "thimble_TRA_output.tsv"
    trb_out = stitchr_dir / "thimble_TRB_output.tsv"

    alpha_map: Dict[str, str] = {}
    beta_map:  Dict[str, str] = {}

    if not df_need.empty:
        write_thimble_inputs(df_need, tra_in, trb_in)
        if args.reuse_thimble and tra_out.exists() and trb_out.exists():
            print("[stage1] Reusing existing Thimble outputs.")
        else:
            print("[stage1] Running Thimble (TRA)…")
            run_thimble(tra_in, tra_out, receptor_flag="a", species=args.species)
            print("[stage1] Running Thimble (TRB)…")
            run_thimble(trb_in, trb_out, receptor_flag="b", species=args.species)

        try:
            A = parse_thimble_output(tra_out, chain="alpha")
            B = parse_thimble_output(trb_out, chain="beta")
            alpha_map = {str(r["row_id"]): str(r["aa"]) for _, r in A.iterrows()}
            beta_map  = {str(r["row_id"]): str(r["aa"]) for _, r in B.iterrows()}
        except Exception as e:
            print(f"[stage1] WARNING: failed to parse Thimble outputs ({e}); "
                  f"only rows with direct AA will be written.", file=sys.stderr)

    # Emit FASTAs (resume-safe)
    counts = Counter()
    # Rows with direct AA first (so users see immediate output even if Thimble fails)
    write_fastas(df_direct, fasta_dir, alpha_map={}, beta_map={}, summary_counts=counts)
    # Rows reconstructed by Thimble
    if alpha_map or beta_map:
        write_fastas(df_need, fasta_dir, alpha_map=alpha_map, beta_map=beta_map, summary_counts=counts)

    # Final per-file summary + global totals
    summary_df = pd.DataFrame(global_summary_rows)
    summary_df.to_csv(work_dir / "stage1_summary.csv", index=False)

    totals = {
        "total_clonotypes_kept": len(df_all),
        "with_direct_full_AA": len(df_direct),
        "submitted_to_thimble": len(df_need),
        "written_fastas": counts.get("written_fastas", 0),
        "skipped_missing_aa": counts.get("skip_missing_aa", 0),
        "fasta_dir": str(fasta_dir),
        "thimble_out_TRA": str(tra_out),
        "thimble_out_TRB": str(trb_out),
        "summary_csv": str(work_dir / "stage1_summary.csv"),
    }
    print("[stage1] Summary:")
    for k, v in totals.items():
        print(f"  - {k}: {v}")

    # Exit code reflects whether anything usable was produced
    if totals["written_fastas"] == 0:
        sys.exit("[stage1] No FASTAs were written. Check stage1_summary.csv and Thimble logs.")

if __name__ == "__main__":
    main()

