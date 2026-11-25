#!/bin/bash
#SBATCH --job-name=sptce_E_clustcr
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=3G
#SBATCH --output=logs/stage_E_%j.out
#SBATCH --error=logs/stage_E_%j.err
# Log paths are relative to the directory you sbatch from (the repository
# root). Check partition/time/cpus/mem against your cluster, and add
# --account=<project> if your cluster requires one.

set -euo pipefail

# Locate the repo root. Under sbatch, SLURM copies this script into a spool
# dir, so BASH_SOURCE/$0 point there, not at sbatch/. SLURM_SUBMIT_DIR is the
# directory you ran `sbatch` from (the repo root, per the documented workflow).
# When run by hand (no SLURM_SUBMIT_DIR), fall back to the script's own path.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    REPO_ROOT="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)"
fi
# If submitted from sbatch/ instead of the repo root, correct for it.
[[ -f "${REPO_ROOT}/config.sh" ]] || REPO_ROOT="$(dirname "${REPO_ROOT}")"

cd "${REPO_ROOT}"
echo "[stage E] running in: $(pwd)"
if [[ ! -f config.sh ]]; then
    echo "[stage E] ERROR: config.sh not found in ${REPO_ROOT}" >&2
    exit 1
fi
source config.sh

mkdir -p logs "${WORK_DIR}/clustcr"

# clusTCR/faiss live in clustcr_env. Keep ~/.local site-packages out of the
# path, so package versions come from the environment; this is exported on its
# own line, before Python starts.
export PYTHONNOUSERSITE=1

# Honour SLURM's CPU allocation rather than grabbing the whole node.
NCPUS="${CLUSTCR_NCPUS:-${SLURM_CPUS_PER_TASK:-all}}"

"${PY_CLUSTCR}" "${SPTCE_ROOT}/scripts/run_clustcr.py" \
    --input "${INPUT_CSV}" \
    --outprefix "${WORK_DIR}/clustcr/clustcr_${CLUSTCR_CHAIN:-beta}" \
    --method "${CLUSTCR_METHOD:-mcl}" \
    --n-cpus "${NCPUS}" \
    --chain "${CLUSTCR_CHAIN:-beta}" \
    --cdr3b-col "${CLUSTCR_CDR3B_COL:-cdr3b}" \
    --vb-col "${CLUSTCR_VB_COL:-v_b_gene}" \
    --cdr3a-col "${CLUSTCR_CDR3A_COL:-cdr3a}" \
    --va-col "${CLUSTCR_VA_COL:-v_a_gene}" \
    ${CLUSTCR_ID_COL:+--id-col "${CLUSTCR_ID_COL}"}

echo "Stage E (clustcr) complete -> ${WORK_DIR}/clustcr/"
