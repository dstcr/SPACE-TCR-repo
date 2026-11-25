#!/bin/bash
#SBATCH --job-name=sptce_C
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --partition=long
#SBATCH --output=logs/stage_C_%j.out
#SBATCH --error=logs/stage_C_%j.err

set -euo pipefail
: "${SPTCE_CONFIG:=${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
source "$SPTCE_CONFIG"
echo "[config] $SPTCE_CONFIG"

module purge
module load Anaconda3

# Use the copy of SPACE2TCR in this repository; $DATA/TCR/SPACE2TCR is used
# only if that copy is missing.
if [[ -f "$SPTCE_ROOT/SPACE2TCR/util.py" ]]; then
  export PYTHONPATH="$SPTCE_ROOT:${PYTHONPATH:-}"
  echo "[Stage C] SPACE2TCR: in-repo copy under $SPTCE_ROOT"
else
  export PYTHONPATH="$DATA/TCR/SPACE2TCR:${PYTHONPATH:-}"
  echo "[Stage C] SPACE2TCR: falling back to $DATA/TCR/SPACE2TCR"
fi
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PDB_DIR="$WORK_DIR/pdb"
CLUSTER_DIR="$WORK_DIR/clusters"
BIN_RESULTS_DIR="$CLUSTER_DIR/bin_results"
DIST_DIR="$CLUSTER_DIR/distance_tables"
BINS_JSONL="$CLUSTER_DIR/bins.jsonl"
BINS_INDEX="$CLUSTER_DIR/member_to_bin.tsv"

mkdir -p "$CLUSTER_DIR" "$BIN_RESULTS_DIR" "$DIST_DIR"

# Fail here, with a clear message, rather than part-way through binning.
"$PY_SPACE2" -c 'import SPACE2TCR.util' || {
  echo "[Stage C] FATAL: cannot import SPACE2TCR. Run from the repo root with PYTHONPATH=. -- see INSTALL.md." >&2
  exit 1
}

# ── Step 1: CDR3-length binning ───────────────────────────────────────────────
echo "[Stage C] Step 1/3: binning..."
"$PY_SPACE2" -u "$SPTCE_ROOT/scripts/make_bins.py" \
  --pdb-dir       "$PDB_DIR" \
  --selection-key CDR_all \
  --out-manifest  "$BINS_JSONL" \
  --out-index     "$BINS_INDEX" \
  --progress-every 50000

N_BINS=$(wc -l < "$BINS_JSONL")
echo "[Stage C] Bins: $N_BINS"

# ── Step 2: RMSD clustering + pairwise distance tables ───────────────────────
# --write-distance-table writes a gzipped TSV per bin with every pairwise RMSD.
# These are the raw distance values you can load in R for downstream analysis.
echo "[Stage C] Step 2/3: clustering + distance tables..."
"$PY_SPACE2" -u "$SPTCE_ROOT/scripts/run_bins_multicore.py" \
  --python           "$PY_SPACE2" \
  --cluster-one-bin  "$SPTCE_ROOT/scripts/cluster_one_bin.py" \
  --bins-jsonl       "$BINS_JSONL" \
  --pdb-base-dir     "$PDB_DIR" \
  --out-dir          "$BIN_RESULTS_DIR" \
  --cutoff           "$CUTOFF" \
  --prune            "$CUTOFF" \
  --selection-key    CDR_all \
  --anchors-key      fw_all \
  --distance-mode    rmsd \
  --max-exact-bin    10000 \
  --chunk-pairs      100000 \
  --workers          "$SLURM_CPUS_PER_TASK" \
  --write-distance-table \
  --distance-dir     "$DIST_DIR"

echo "[Stage C] Distance tables written to: $DIST_DIR"
echo "[Stage C] $(ls $DIST_DIR/*.tsv.gz 2>/dev/null | wc -l) distance table files"

# ── Step 3: merge bin CSVs into single cluster assignments TSV ───────────────
echo "[Stage C] Step 3/3: merging..."
export BIN_RESULTS_DIR
export CLUSTER_ASSIGNMENTS="$CLUSTER_DIR/cluster_assignments.tsv"

"$PY_SPACE2" - << 'PYEOF'
import os, glob, pandas as pd
from pathlib import Path
bin_dir  = os.environ["BIN_RESULTS_DIR"]
out_file = os.environ["CLUSTER_ASSIGNMENTS"]
csvs = sorted(glob.glob(f"{bin_dir}/bin_*.csv"))
if not csvs:
    raise SystemExit(f"No bin CSVs in {bin_dir}")
merged = pd.concat([pd.read_csv(f) for f in csvs], ignore_index=True)
merged = merged.rename(columns={
    "ID": "tcr_id",
    "cluster_by_rmsd": "cluster_id",
    "cluster_by_length": "bin_key",
})
merged["tcr_id"] = merged["tcr_id"].apply(lambda x: Path(x).stem)
merged["cluster_size"] = merged.groupby("cluster_id")["tcr_id"].transform("count")
merged.to_csv(out_file, sep="\t", index=False)
clustered = (merged["cluster_size"] > 1).sum()
print(f"TCRs: {len(merged):,}  Clusters: {merged['cluster_id'].nunique():,}  Clustered: {clustered:,} ({clustered/len(merged):.1%})")
PYEOF

echo "[Stage C] Done."
echo "[Stage C] Cluster assignments: $CLUSTER_DIR/cluster_assignments.tsv"
echo "[Stage C] Distance tables:     $DIST_DIR/"
