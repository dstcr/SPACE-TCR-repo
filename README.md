# SPACE-TCR

Structure-guided clustering of paired αβ T-cell receptor repertoires. Receptors
are grouped by the predicted geometry of their CDR loops rather than by sequence
similarity, which recovers convergent families that sequence methods miss.

Devinder Kaur Srai, Oxford DPhil 2025.

---

## Reproduce the paper

**You need:** [R ≥ 4.1](https://cloud.r-project.org/) and git. No Python, GPU or
cluster is needed.

```bash
git clone --depth 1 https://github.com/dstcr/SPACE-TCR-repo.git
cd SPACE-TCR-repo
./install.sh                 # installs 7 R packages, then self-tests (a few minutes)
Rscript validate_outputs.R   # recomputes the headline numbers (under a minute)
```

The last line of output should read:

```
  36 PASS   0 FAIL   0 SKIP
```

The download is about 70 MB.
**No git?** On GitHub, *Code → Download ZIP*, unzip, and run the same two
commands inside the folder. **On Windows**, `./install.sh` needs Git Bash or WSL;
from a normal terminal use `Rscript install/r_packages.R` in its place.
If anything fails, `Rscript check_install.R` says what is missing and how to
install it.

Both benchmarks' inputs and their structural geometry (pairwise RMSD tables and
CDR3 dihedrals) are included, so the numbers reproduce without predicting any
structures. `validate_outputs.R` recomputes the headline figures from
[`data/`](data/README.md) with the same metric code the pipeline uses and
compares them with [`results/`](results/README.md): at the focal operating
point, 149 clusters, 681 receptors clustered, 94.63% of clusters at purity ≥ 0.9
and pooled purity 0.9868. Its exit status is the number of failures.

## Data

| What | Where |
|---|---|
| Labelled cohort, Stage C/D geometry for both benchmarks | [`data/`](data/README.md), in this repository |
| Raw data: the modelled structures (TCRBuilder2 PDBs) | [Zenodo record 23765886](https://zenodo.org/records/23765886) |
| Published methods' cluster output | the Jouannet benchmark repository, see [`third_party/JOUANNET.md`](third_party/JOUANNET.md) |

The structures are only needed to rerun Stages C and D, or to inspect
individual models.

## Run it on your own receptors

```bash
./install.sh --own-data
```

Then point `CFG` at the top of `sptce_pipeline.R` at your data —
`jouannet_labelled` (one row per receptor: `id`, CDR3s, `epitope`) and
`jouannet_sptce_out` (a folder holding `distance_tables/` and
`dihedrals.tsv.gz`) — and:

```r
source("sptce_pipeline.R"); run_benchmark_b()
```

To start from sequences rather than structures, use the five-stage cluster
pipeline described in Track B of [`INSTALL.md`](INSTALL.md).

| Stage | In | Out |
|---|---|---|
| A | `$INPUT_CSV` | stitched FASTAs (stitchr) |
| B | FASTAs | PDBs (100-task array, TCRBuilder2) |
| C | PDBs | RMSD clusters + pairwise distance tables |
| D | PDBs | CDR3 backbone dihedrals (parallel with C) |
| E | `$INPUT_CSV` | clusTCR sequence clusters (independent arm) |

```bash
sbatch sbatch/run_A_slurm.sh        # then B, C, D, E
bash sbatch/submit_all_subsets.sh   # or chain A→B→C and A→B→D with dependencies
```

`config.sh` is the one file to edit; every job script reads it.

## Figures

```bash
Rscript make_figures.R ~/paper/figures
```

Renders every figure your inputs support and writes a manifest of what was
produced and what was skipped.

---

## Where to go

| You want to | Read |
|---|---|
| Install anything | [`INSTALL.md`](INSTALL.md) |
| Understand the inputs | [`data/README.md`](data/README.md) |
| See the rendered figures | [`results/README.md`](results/README.md) |
| Predict structures (Stage B) | [`install/TCRBUILDER2.md`](install/TCRBUILDER2.md) |
| Distil the clustering into a sequence model | [`ml/README.md`](ml/README.md) |
| Run the noise ladder | [`docs/BENCHMARK_A.md`](docs/BENCHMARK_A.md), then `Rscript run_benchmark_a.R` |
| Get the Jouannet benchmark repo | [`third_party/JOUANNET.md`](third_party/JOUANNET.md) |

## Layout

```
install.sh              one-command install (wraps install/)
INSTALL.md              the full walkthrough, laptop and cluster
sptce_metrics.R         the metrics: clustering, purity, retention, sensitivity
sptce_pipeline.R        the pipeline: runs the benchmarks and draws the figures
config.sh               paths and parameters for the Slurm stages
validate_outputs.R      recompute the headline numbers and compare with results/
make_figures.R          render every paper figure
check_install.R  .sh    verify an install (the .R one works on Windows)

SPACE2TCR/              the clustering library (vendored, BSD-3)
scripts/  sbatch/       the five pipeline stages, and the figure scripts
data/                   both benchmarks' inputs and their Stage C/D geometry
results/                rendered figures and tables
ml/                     sequence-only distillation of the clustering
tests/                  R and Python tests; no data needed
third_party/            how to obtain the Jouannet benchmark repo
docs/                   benchmark notes
```

`sptce_metrics.R` follows the reference implementation of Jouannet et al., so
SPACE-TCR and the published methods are measured by the same equations. The
pipeline, `validate_outputs.R`, the scripts and the tests all use it.

The dihedral gate is `CFG$focal_dihedral`, in radians; the paper uses 0.5.
Output filenames include the value.

## Licence

BSD 3-Clause — see [`LICENSE`](LICENSE).

## Citing

The clustering library is a TCR adaptation of SPACE2 and carries its BSD 3-Clause
licence — see [`SPACE2TCR/_upstream/NOTICE.md`](SPACE2TCR/_upstream/NOTICE.md),
which records what the fork changed.

- Spoendlin et al. (2023), *SPACE2* — the clustering framework
- Abanades et al. (2023), *ImmuneBuilder / TCRBuilder2* — structure prediction
- Quast et al. (2025), *TCRBuilder2+* — the structure prediction model used here
- Heather et al. (2022), *Stitchr* — sequence reconstruction
- Jouannet et al. (2025) — the benchmark cohort and metric definitions
