#!/usr/bin/env bash
# Build the constant-alpha arm of Figure S9.
#
#   bash scripts/s9_constant_alpha.sh [pdb_dir]
#
# Re-models the 1,523 receptors of bin CDRB(10,8,13) with one fixed alpha
# partner (tcr1's) and their own beta chain. The published Figure S9 used the
# SPACE2TCR distance tables of that run, which are committed in
# results/tables/extensions/s9_constant_alpha/; this script rebuilds the models
# and writes the beta-arm table and the dihedral table:
#
#   results/tables/extensions/s9_singlechain_constant_alpha.csv
#   results/tables/extensions/s9_constant_alpha_dihedrals.tsv.gz
#
# Needs: the modelled cohort structures (Zenodo record 23765886, 1.4 GB, 7z),
# ImmuneBuilder (pip install ImmuneBuilder; see install/TCRBUILDER2.md) and the
# SPACE2TCR environment. Modelling takes about 20-40 s per receptor on CPU, so
# expect several hours; `model` skips receptors already built, so the step can
# be rerun after an interruption.
set -euo pipefail
cd "$(dirname "$0")/.."

PDB_DIR="${1:-work/cohort_pdb}"
WORK=work/s9
EXT=results/tables/extensions
IDS=$EXT/s9_bin_CDRB_10_8_13_ids.txt
mkdir -p "$WORK"

if [[ ! -d "$PDB_DIR" ]]; then
  echo "No structures at $PDB_DIR. Download SPACE_paper_raw_pbs_2.7z from"
  echo "https://zenodo.org/records/23765886, extract it, and pass the folder"
  echo "holding the cohort's tcr*.pdb files as the first argument."
  exit 1
fi

# 1. Both chain sequences, recovered from the cohort models.
python3 scripts/dummy_alpha_control.py seqs --pdb-dir "$PDB_DIR" --out "$WORK/seqs.tsv"

# 2. Re-model the bin with a constant alpha partner (resumable).
python3 scripts/dummy_alpha_control.py model --seqs "$WORK/seqs.tsv" --ids "$IDS" \
  --mode dummy --dummy-id tcr1 --out-dir "$WORK/m_dummy"

# 3. Beta loops on the beta framework only, over the figure's thresholds.
PYTHONPATH=. python3 scripts/chain_framework_ablation.py --pdb-dir "$WORK/m_dummy" \
  --labels data/jouannet_labelled_input.tsv --arms beta_fwB \
  --subset-ids "$IDS" --thresholds 0.4,0.5,0.6,0.75,0.9,1.0 \
  --out $EXT/s9_singlechain_constant_alpha.csv

# 4. CDR3 backbone dihedrals of the constant-alpha models.
python3 scripts/extract_dihedrals.py --pdb-dir "$WORK/m_dummy" \
  --out-tsv $EXT/s9_constant_alpha_dihedrals.tsv.gz --workers 2

echo "Done. Now: python3 scripts/fig_s9_constant_alpha.py"
