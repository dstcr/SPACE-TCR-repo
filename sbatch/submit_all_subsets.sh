#!/bin/bash
# submit_all_subsets.sh — submit Stage A→B→C+D for each subset, with deps
set -euo pipefail
cd "$(dirname "$0")/.."

for i in 1; do
  CFG="$PWD/config.sh"
  echo "=== Submitting subset ${i} ==="

  A=$(sbatch --parsable \
        --export=ALL,SPTCE_CONFIG=$CFG \
        --job-name=sptce_A_s${i} \
        sbatch/run_A_slurm.sh)
  echo "  A: $A"

  B=$(sbatch --parsable --dependency=afterok:$A \
        --export=ALL,SPTCE_CONFIG=$CFG \
        --job-name=sptce_B_s${i} \
        sbatch/run_B_slurm.sh)
  echo "  B: $B (waits on A)"

  C=$(sbatch --parsable --dependency=afterok:$B \
        --export=ALL,SPTCE_CONFIG=$CFG \
        --job-name=sptce_C_s${i} \
        sbatch/run_C_slurm.sh)
  echo "  C: $C (waits on B)"

  D=$(sbatch --parsable --dependency=afterok:$B \
        --export=ALL,SPTCE_CONFIG=$CFG \
        --job-name=sptce_D_s${i} \
        sbatch/run_D_slurm.sh)
  echo "  D: $D (waits on B, parallel with C)"
done
