#!/usr/bin/env python3
"""
Stage 2 — run TCRBuilder2 on FASTAs produced by stage1_stitchr.py

Input FASTAs are expected to contain:
>B
<beta sequence>
>A
<alpha sequence>

Outputs:
- one PDB per FASTA into --pdb-dir
- one failure log per Slurm array task
"""

import os
import shutil
from pathlib import Path
from typing import Dict, Tuple

from Bio import SeqIO
from ImmuneBuilder import TCRBuilder2

AA_STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta_AB(fasta_path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        hdr = rec.id.strip()
        seq = str(rec.seq).strip()
        if hdr == "A" or hdr.endswith("|A"):
            seqs["A"] = seq
        elif hdr == "B" or hdr.endswith("|B"):
            seqs["B"] = seq
        else:
            # tolerate simple headers containing A/B somewhere
            parts = hdr.split("|")
            for p in reversed(parts):
                if p in ("A", "B"):
                    seqs[p] = seq
                    break
    return seqs


def is_valid_chain(seq: str, min_len: int = 100, max_len: int = 500) -> Tuple[bool, str]:
    if seq is None:
        return False, "seq_none"
    s = seq.strip().upper()
    if not s:
        return False, "seq_empty"
    bad = set(s) - AA_STANDARD
    if bad:
        return False, f"bad_chars={''.join(sorted(bad))}"
    if len(s) < min_len:
        return False, f"too_short={len(s)}"
    if len(s) > max_len:
        return False, f"too_long={len(s)}"
    return True, "ok"


def pdb_ok(p: Path) -> bool:
    if not p.exists():
        return False
    try:
        if p.stat().st_size < 2000:
            return False
        with p.open() as fh:
            for _ in range(400):
                line = fh.readline()
                if not line:
                    break
                if line.startswith("ATOM"):
                    return True
        return False
    except Exception:
        return False


def main():
    fasta_dir = Path(os.environ.get("FASTA_DIR", "fasta_aa")).resolve()
    pdb_dir = Path(os.environ.get("PDB_DIR", "pdb")).resolve()
    pdb_dir.mkdir(parents=True, exist_ok=True)

    if shutil.which("hmmscan") is None:
        raise SystemExit("hmmscan not found in PATH. Load HMMER module in the Slurm script.")

    fastas = sorted(fasta_dir.glob("*.fasta"))
    if not fastas:
        raise SystemExit(f"No FASTAs found in {fasta_dir}")

    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    task_count = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", "1"))

    my_fastas = [p for i, p in enumerate(fastas) if (i % task_count) == task_id]
    if not my_fastas:
        print(f"Task {task_id}: nothing to do.")
        return

    predictor = TCRBuilder2()

    failed_path = pdb_dir / f"tcrbuilder_failed_{task_id}.tsv"
    wrote = 0
    skipped_ok = 0
    removed_bad = 0
    failed = 0
    invalid = 0
    missing_chain = 0

    with failed_path.open("w") as bad:
        bad.write("fasta_name\treason\tdetail\tlenA\tlenB\tfasta\n")

        for fasta_path in my_fastas:
            uid = fasta_path.stem
            pdb_path = pdb_dir / f"{uid}.pdb"

            if pdb_ok(pdb_path):
                skipped_ok += 1
                continue

            if pdb_path.exists() and not pdb_ok(pdb_path):
                removed_bad += 1
                try:
                    pdb_path.unlink()
                except Exception:
                    pass

            seqs = read_fasta_AB(fasta_path)
            if "A" not in seqs or "B" not in seqs:
                missing_chain += 1
                bad.write(
                    f"{uid}\tmissing_chain\tA_or_B_missing\t"
                    f"{len(seqs.get('A',''))}\t{len(seqs.get('B',''))}\t{fasta_path}\n"
                )
                continue

            okA, whyA = is_valid_chain(seqs["A"])
            okB, whyB = is_valid_chain(seqs["B"])
            if not okA or not okB:
                invalid += 1
                bad.write(
                    f"{uid}\tinvalid_sequence\tA:{whyA};B:{whyB}\t"
                    f"{len(seqs['A'])}\t{len(seqs['B'])}\t{fasta_path}\n"
                )
                continue

            try:
                tcr = predictor.predict({"A": seqs["A"], "B": seqs["B"]})
                tcr.save(str(pdb_path))
                wrote += 1
            except AssertionError as e:
                failed += 1
                msg = str(e).replace("\t", " ").replace("\n", " ")
                bad.write(
                    f"{uid}\tassertion_error\t{msg}\t"
                    f"{len(seqs['A'])}\t{len(seqs['B'])}\t{fasta_path}\n"
                )
            except Exception as e:
                failed += 1
                msg = f"{type(e).__name__}: {e}".replace("\t", " ").replace("\n", " ")
                bad.write(
                    f"{uid}\texception\t{msg}\t"
                    f"{len(seqs['A'])}\t{len(seqs['B'])}\t{fasta_path}\n"
                )

    print(
        f"Task {task_id}: assigned={len(my_fastas)} "
        f"wrote={wrote} skipped_ok={skipped_ok} removed_bad_pdb={removed_bad} "
        f"missing_chain={missing_chain} invalid={invalid} failed={failed} "
        f"failed_log={failed_path.name}"
    )


if __name__ == "__main__":
    main()
