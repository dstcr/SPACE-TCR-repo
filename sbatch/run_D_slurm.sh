#!/bin/bash
#SBATCH --job-name=sptce_D
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=short
#SBATCH --output=logs/stage_D_%j.out
#SBATCH --error=logs/stage_D_%j.err

# Stage D: extract raw CDR3 backbone dihedral angles (phi/psi) for every PDB.
# Output: one row per (TCR x chain x CDR3 residue) with phi_rad, psi_rad,
#         phi_deg, psi_deg, atoms_present.
# Can run in parallel with Stage C — only needs Stage B (PDBs) to be complete.

set -euo pipefail
: "${SPTCE_CONFIG:=${SLURM_SUBMIT_DIR:-$PWD}/config.sh}"
source "$SPTCE_CONFIG"
echo "[config] $SPTCE_CONFIG"

module purge
module load Anaconda3

# Uses the space2tcr env (numpy available); no special dependencies needed.
PDB_DIR="$WORK_DIR/pdb"
OUT_DIR="$WORK_DIR/dihedrals"
OUT_TSV="$OUT_DIR/dihedrals.tsv.gz"

mkdir -p "$OUT_DIR"

N_PDBS=$(ls "$PDB_DIR"/*.pdb 2>/dev/null | wc -l)
echo "[Stage D] PDBs found: $N_PDBS"
echo "[Stage D] Output: $OUT_TSV"
echo "[Stage D] Workers: $SLURM_CPUS_PER_TASK"
echo "[Stage D] CDR3 residues: IMGT 105-117 (both chains)"

"$PY_SPACE2" "$SPTCE_ROOT/scripts/extract_dihedrals.py" \
  --pdb-dir        "$PDB_DIR" \
  --out-tsv        "$OUT_TSV" \
  --cdr3-start     105 \
  --cdr3-end       117 \
  --workers        "$SLURM_CPUS_PER_TASK" \
  --progress-every 500

echo ""
echo "[Stage D] Done."
echo "[Stage D] Output: $OUT_TSV"
echo "[Stage D] Rows: $(zcat $OUT_TSV | wc -l) (including header)"
echo ""
echo "Load in R:"
echo "  dihedrals <- readr::read_tsv('$OUT_TSV')"
echo "  # Columns: tcr_id, chain, resi, phi_rad, psi_rad, phi_deg, psi_deg, atoms_present"
