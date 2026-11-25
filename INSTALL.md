# Installing SPACE-TCR

## Quick start

Prerequisites: **R ≥ 4.1** and **git**. If you don't have R yet, see
[Track A step 1](#1-r) below.

```bash
git clone --depth 1 https://github.com/dstcr/SPACE-TCR-repo.git
cd SPACE-TCR-repo
./install.sh                 # installs 7 R packages, then self-tests (a few minutes)
Rscript validate_outputs.R   # reproduces the paper's headline numbers (under a minute)
```

Success looks like `36 PASS   0 FAIL   0 SKIP` on the last line.

That is all of **Track A**. It needs only R and reproduces the paper's numbers
from the geometry in `data/`, without predicting any structures.

- **No git?** On GitHub, *Code → Download ZIP*, unzip, and run the same
  commands inside the folder. Everything needed is in the ZIP.
- **Windows:** `./install.sh` is a bash script. Use Git Bash or WSL, or replace
  it with `Rscript install/r_packages.R` in an ordinary terminal (or
  `source("install/r_packages.R")` from RStudio with the repository as working
  directory). `validate_outputs.R` works as-is.
- **Something failed?** `Rscript check_install.R` lists what is missing and the
  command to install it.

| You want to | Run | Time |
|---|---|---|
| Reproduce the paper | `./install.sh` | a few minutes |
| Cluster your own structures | `./install.sh --own-data` | ~8 min |
| Also the ML distillation | `./install.sh --everything` | ~13 min |
| Check an existing install | `./install.sh --check` | ~1 min |

`install.sh` can be re-run safely, exits non-zero on failure, and wraps the
installers under `install/`, which can also be used on their own:

| What | On its own | Needs |
|---|---|---|
| R — analysis and figures | `Rscript install/r_packages.R [--figures]` | R ≥ 4.1 |
| SPACE2-TCR — clustering (Stages C, D) | `bash install/space2tcr.sh` | Python 3.10–3.13 |
| Distillation model (`ml/`) | `bash install/distill.sh --venv .venv-ml` | Python 3.10–3.13 |
| TCRBuilder2 — prediction (Stage B) | [`install/TCRBUILDER2.md`](install/TCRBUILDER2.md) | conda |

`install.sh` does not install TCRBuilder2 or stitchr: each needs a conda
environment and a licence acceptance. Neither is needed to reproduce the paper;
Track B below covers them.

Each installer finishes with a check. The SPACE2-TCR installer computes RMSDs
and confirms they are zero for identical structures, unchanged under a random
rigid-body transform, and larger when a CDR loop is perturbed (22 checks). The
R installer loads the metric and pipeline files. The distillation installer
checks `ml/metrics.py` against `sptce_metrics.R` and then trains a model on
synthetic data, which must beat the Levenshtein baseline.

---

# The full picture

Two tracks, independent of each other:

| | **Track A — analysis** | **Track B — pipeline** |
|---|---|---|
| Where | a laptop | an HPC cluster with Slurm |
| What it does | scores clusters, reproduces the tables | sequences → structures → clusters |
| Needs | R and 7 packages | 4 conda environments, HMMER, ~50 GB, GPU optional |
| Setup time | a few minutes, plus installing R if needed | half a day |

Most people want Track A. It is also how you confirm a download is complete.

To check an install at any point:

```bash
Rscript check_install.R          # the R side; works on Windows
bash check_install.sh            # everything, including Python/SPACE2TCR (Linux/macOS)
bash check_install.sh --analysis # Track A only
bash check_install.sh --pipeline # Track B only
```

Both exit non-zero on failure and print the install command for anything missing.

> **On Windows**, use `check_install.R`. The `.sh` scripts need Git Bash or WSL
> and add only the Python/SPACE2TCR checks, which matter on the machine that
> runs Stage C. Track A is pure R and works on any platform. From RStudio you
> can also run `source("check_install.R")`.

---

## Track A — the analysis

### 1. R

R ≥ 4.1; tested with R 4.3.3.

```bash
# Debian / Ubuntu
sudo apt-get install -y r-base-core

# macOS
brew install r                     # or the installer from cran.r-project.org

# conda, to leave the system R alone
conda create -n sptce -c conda-forge r-base=4.3 && conda activate sptce
```

### 2. Seven R packages

```r
install.packages(c("readr","dplyr","tidyr","purrr","stringr","tibble","stringdist"))
```

On Debian/Ubuntu the distribution packages are faster and avoid compiling:

```bash
sudo apt-get install -y r-cran-dplyr r-cran-readr r-cran-tidyr r-cran-purrr \
                        r-cran-stringr r-cran-tibble r-cran-stringdist
```

These are all `sptce_metrics.R`, `validate_outputs.R` and the tests need.

### 3. Check

```bash
Rscript check_install.R            # any platform, including Windows
# or, on Linux/macOS:
bash check_install.sh --analysis
```

Expected:

```
  [ok]   all 7 required R packages present
  [ok]   sptce_metrics.R parses
  [ok]   metric unit tests: 41 passed, 0 failed
  [ok]   end-to-end tests: e2e: 0 failures
```

Or run the tests directly:

```bash
Rscript tests/test_metrics.R   # 41 unit tests against hand-computed values
Rscript tests/test_e2e.R       # 17 checks on a synthetic 60-receptor cohort
```

Neither needs any data. They cover the metric definitions: purity is `>= 0.9`,
retention counts distinct receptors, the sensitivity rule uses the reference's
`0.67`, and the self-checks reject an out-of-range sensitivity.

### 4. For the pipeline: the plotting packages

`sptce_pipeline.R` runs the benchmarks and draws the figures, which needs:

```r
install.packages(c("ggplot2","patchwork","scales","forcats","reshape2","igraph",
                   "ggraph","tidygraph","cowplot","ggseqlogo","ggforce",
                   "ggrepel","kableExtra"))
```

`./install.sh --own-data` installs these for you. `igraph`, `ggraph` and
`tidygraph` build against system libraries; on Debian/Ubuntu install
`libglpk-dev libxml2-dev libudunits2-dev libgdal-dev` first if they fail to
build. Validation and the tests do not need these packages.

### 5. Point it at your data

Edit `CFG` at the top of `sptce_pipeline.R`:

```r
CFG$jouannet_labelled   # labelled input TSV, one row per receptor (id, CDR3s, epitope)
CFG$jouannet_sptce_out  # folder holding distance_tables/*.tsv.gz and dihedrals.tsv.gz
CFG$jouannet_repo       # the Jouannet benchmark repository, for the published methods
CFG$jouannet_out_dir    # where tables and figures are written
CFG$focal_dihedral      # the dihedral gate, in radians
```

The paper uses a dihedral gate of 0.5 rad; 0.3 is a tighter alternative. The
value is part of the output filenames. The metric definitions themselves
(cluster size, purity cut-off, sensitivity rule) are `METRIC` in
`sptce_metrics.R`.

Then, from the repository root:

```r
source("sptce_pipeline.R")
run_benchmark_b()      # curated cohort: tables, accounting, consensus, figures
```

---

## Track B — the cluster pipeline

Only needed to go from sequences to structures to clusters. If you already have
`distance_tables/` and `dihedrals.tsv.gz`, you do not need it. To re-cluster
the paper's cohort without predicting structures, download the modelled PDBs
from [Zenodo record 23765886](https://zenodo.org/records/23765886) and start at
Stage C.

### 1. Clone and configure

```bash
git clone --depth 1 https://github.com/dstcr/SPACE-TCR-repo.git && cd SPACE-TCR-repo
nano config.sh      # set CLUSTER_DATA (defaults to $HOME), INPUT_CSV and the four environment paths
```

`config.sh` is the one file to edit; every job script reads it. To use a
different configuration for one submission:

```bash
sbatch --export=ALL,SPTCE_CONFIG=$PWD/config.other.sh sbatch/run_A_slurm.sh
```

### 2. Four conda environments

They are separate because their dependencies conflict: TCRBuilder2 needs torch,
clusTCR needs faiss, and SPACE2TCR needs a numba version matched to its numpy.

```bash
# Stage A — stitchr
conda create -n stitchr python=3.10 pandas -y
conda activate stitchr && pip install stitchr IMGTgeneDL

# Stage B — TCRBuilder2
conda create -n myenv2TCR python=3.10 -y
conda activate myenv2TCR && pip install ImmuneBuilder biopython
# ANARCI needs HMMER; on ARC, load the module
module load HMMER/3.3.2-gompi-2022a

# Stages C and D — SPACE2TCR
conda create -n space2tcr python=3.10 -y
conda activate space2tcr
pip install "numba>=0.51" "joblib>=1.0" scipy pandas "scikit-learn>=1.4"
# SPACE2TCR itself is included in this repository; there is nothing to pip install.

# Stage E — clusTCR
conda create -n clustcr_env -c bioconda -c conda-forge clustcr -y
```

Then set the four paths in `config.sh`:

```bash
STITCHR_ENV="/path/to/envs/stitchr"
PY_TCRBUILDER="/path/to/envs/myenv2TCR/bin/python"
PY_SPACE2="/path/to/envs/space2tcr/bin/python"
PY_CLUSTCR="/path/to/envs/clustcr_env/bin/python"
```

To record the exact versions in use, run `bash collect_from_cluster.sh` on the
cluster: it writes `pip freeze` and `conda env export` for all four
environments. `numba` and `numpy` matter most, because SPACE2TCR's compiled
cache is tied to both.

### 3. Check SPACE2TCR

```bash
PYTHONPATH=. $PY_SPACE2 tests/verify_space2tcr.py .
```

This computes RMSDs for identical structures, a rigid-body transform, a CDR
perturbation and a framework-only change, and checks the region definitions and
the binning. If it passes, Stage C will run. `run_C_slurm.sh` also checks the
import before it starts.

### 4. Input data

Your own input table is not tracked: `data/` is ignored by git apart from the
two benchmarks shipped with the paper (see [`data/README.md`](data/README.md)).
Put the Stage A input table there and name it in `config.sh` as `INPUT_CSV`
(a path relative to `SPTCE_ROOT`).

```bash
scp clonotypes_1.tsv user@cluster:/path/to/SPACE-TCR-repo/data/
```

Use filenames without spaces or parentheses.

### 5. Run

```bash
sbatch sbatch/run_A_slurm.sh     # CSV -> stitched FASTAs
sbatch sbatch/run_B_slurm.sh     # FASTAs -> PDBs (array of 100)
sbatch sbatch/run_C_slurm.sh     # PDBs -> clusters + distance tables
sbatch sbatch/run_D_slurm.sh     # PDBs -> CDR3 dihedrals (parallel with C)
sbatch sbatch/run_E_slurm.sh     # optional sequence-clustering arm

bash sbatch/submit_all_subsets.sh  # or chain A->B->C and A->B->D with dependencies
```

Submit from the repository root. Stage B takes longest: TCRBuilder2 needs
roughly a minute per receptor on CPU and about five seconds on GPU.

### 6. Bring the results back for Track A

The analysis reads two things per cohort:

```
<WORK_DIR>/clusters/distance_tables/*.tsv.gz   -> CFG$jouannet_sptce_out
<WORK_DIR>/dihedrals/dihedrals.tsv.gz          -> CFG$jouannet_sptce_out
```

Copy that folder to your machine and point `CFG$jouannet_sptce_out` at it.

---

## The Jouannet benchmark repository

The pipeline reads the published methods' cluster output and the 10X noise
ladder from it. It is about 382 MB, so keep it beside this repository rather
than inside it:

```
~/work/SPACE-TCR-repo/                    this repository
~/work/TCR_Unsupervised_Benchmark-main/   the benchmark repository
```

`bash third_party/fetch_jouannet.sh` downloads it there. Set `CFG$jouannet_repo`
in `sptce_pipeline.R` accordingly; [`third_party/JOUANNET.md`](third_party/JOUANNET.md)
lists the files that are read.

---

## Tested versions

October 2026, Linux x86-64, Python 3.13.16 and R 4.3.3. The requirements files
give minimum versions; these are the versions the results were checked with.

```
numpy 2.5.3   pandas 3.0.5   scipy 1.18.1   numba 0.68.0
joblib 1.6.0  scikit-learn 1.9.1   torch 2.14.1   rapidfuzz 3.14.6

readr 2.1.5   dplyr 1.1.4   tidyr 1.3.1   purrr 1.0.2   stringr 1.5.1
tibble 3.2.1  stringdist 0.9.12   ggplot2 3.4.4   igraph 1.6.0
ggraph 2.1.0  tidygraph 1.3.0     ggseqlogo 0.1   jsonlite 1.8.8
```

- `tidygraph` ≥ 1.3.1 requires `igraph` ≥ 2.0. With `igraph` 1.6.0 (the
  Debian package), install `tidygraph` 1.3.0 and `ggraph` 2.1.0; the R
  installer prints the commands.
- SPACE2TCR's compiled numba cache (`SPACE2TCR/__pycache__`) is tied to the
  numba and Python versions. Delete it after changing either.

---

## Troubleshooting

| Symptom | Cause |
|---|---|
| `cannot import SPACE2TCR` | run from the repository root with `PYTHONPATH=.`, or check `SPACE2TCR/util.py` exists. |
| numbers differ from the paper | check `CFG$focal_dihedral`: the paper uses 0.5 rad. Output filenames include the value. |
| `check_sensitivity failed` | a sensitivity exceeded the fraction of that epitope's receptors that were clustered, so the input or the denominator is wrong. |
| `Argument 'x' is not a vector: list` | `seqinr` or a Bioconductor package is masking a dplyr function. Start a clean R session. |
| numba cache errors after a version change | delete `SPACE2TCR/__pycache__`. |
| `no space left on device` on the cluster | Stage B writes one PDB per receptor: a few GB for 4,779 receptors, tens of GB for a large cohort. |
| `bash: command not found` / the `.sh` scripts won't run | on Windows, use `Rscript check_install.R`; the shell scripts need Git Bash or WSL. |
