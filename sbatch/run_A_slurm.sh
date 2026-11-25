#!/bin/bash
#SBATCH --job-name=sptce_A
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=logs/stage_A_%j.out
#SBATCH --error=logs/stage_A_%j.err

set -euo pipefail
: "${SPTCE_CONFIG:=${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
source "$SPTCE_CONFIG"
echo "[config] $SPTCE_CONFIG"

module purge
module load Anaconda3
source activate "$STITCHR_ENV"

INPUT_DIR="$WORK_DIR/stage_A_work/input_csv"
mkdir -p "$INPUT_DIR" "$WORK_DIR/fasta" "$WORK_DIR/stitchr" "$WORK_DIR/stage_A_work"

cp "$SPTCE_ROOT/$INPUT_CSV" "$INPUT_DIR/input.csv"

python "$SPTCE_ROOT/scripts/stage1_stitchr.py" \
  --in-dir       "$INPUT_DIR" \
  --stitchr-dir  "$WORK_DIR/stitchr" \
  --fasta-aa-dir "$WORK_DIR/fasta" \
  --work         "$WORK_DIR/stage_A_work" \
  --species      "$SPECIES" \
  --reuse-thimble

echo "[Stage A] Done. FASTAs written: $(ls $WORK_DIR/fasta/*.fasta 2>/dev/null | wc -l)"
