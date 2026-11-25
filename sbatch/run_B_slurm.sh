#!/bin/bash
#SBATCH --job-name=sptce_B
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-99
#SBATCH --mem=56G
#SBATCH --time=48:00:00
#SBATCH --partition=long
#SBATCH --output=logs/stage_B_%A_%a.out
#SBATCH --error=logs/stage_B_%A_%a.err

set -euo pipefail
: "${SPTCE_CONFIG:=${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
source "$SPTCE_CONFIG"
echo "[config] $SPTCE_CONFIG"

module purge
module load Anaconda3
module load HMMER/3.3.2-gompi-2022a

export FASTA_DIR="$WORK_DIR/fasta"
export PDB_DIR="$WORK_DIR/pdb"
mkdir -p "$PDB_DIR"

which hmmscan

echo "[Stage B] task=$SLURM_ARRAY_TASK_ID  fasta=$FASTA_DIR  pdb=$PDB_DIR"

"$PY_TCRBUILDER" "$SPTCE_ROOT/scripts/stage2_tcrbuilder2.py"
