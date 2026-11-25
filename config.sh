# ============================================================
# SPTCE config — edit these lines only
# ============================================================
# The Stage A input. Both benchmark inputs are tracked in this repository:
#   data/jouannet_labelled_input.tsv            4,779 paired receptors, the
#                                               curated CD8 cohort (Benchmark B)
#   data/benchmark_a/donor1_subset6_sptce.tsv   23,529 paired receptors, the
#                                               10X donor-1 superset whose
#                                               restrictions are the six noise
#                                               ladder rungs (Benchmark A)
# Stage A copies this single file to $WORK_DIR/stage_A_work/input_csv/input.csv,
# so a .tsv is fine -- stage1_stitchr.py sniffs the delimiter.
INPUT_CSV="data/jouannet_labelled_input.tsv"
# Your storage area on the cluster; everything below is relative to it.
CLUSTER_DATA="${CLUSTER_DATA:-$HOME}"
WORK_DIR="$CLUSTER_DATA/SPTCE/output/run1"
SPECIES="human"
CUTOFF=1.25

# Paths to the repository and the Python environments on your cluster
SPTCE_ROOT="$CLUSTER_DATA/SPTCE/"
PY_TCRBUILDER="$CLUSTER_DATA/myenv2TCR/bin/python"
PY_SPACE2="$CLUSTER_DATA/envs/space2tcr/bin/python"
STITCHR_ENV="$CLUSTER_DATA/stitchr"   # adjust if different
DATA="$CLUSTER_DATA"                                # adjust if different

# ============================================================
# Stage E — clusTCR (sequence-based clustering, alt. to Stage C)
# ============================================================
PY_CLUSTCR="$CLUSTER_DATA/clustcr_env/bin/python"
CLUSTCR_METHOD="mcl"          # "mcl" or "two-step"
CLUSTCR_NCPUS="all"           # "all" or an integer; SLURM alloc overrides if set
CLUSTCR_CHAIN="beta"          # "beta" or "paired"

# Column names in INPUT_CSV — EDIT to match your TSV header.
CLUSTCR_CDR3B_COL="cdr3_aa_beta"
CLUSTCR_VB_COL="v_call_beta"
CLUSTCR_CDR3A_COL="cdr3_aa_alpha"     # only used when CLUSTCR_CHAIN="paired"
CLUSTCR_VA_COL="v_call_alpha"         # only used when CLUSTCR_CHAIN="paired"
CLUSTCR_ID_COL=""             # optional ID column to carry through; "" = none
