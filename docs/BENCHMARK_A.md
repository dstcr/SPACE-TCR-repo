# Benchmark A — the 10X donor-1 noise ladder

Benchmark A tests whether structural clustering still recovers antigen-specific
families when most of the repertoire is unlabelled background. The labelled
donor-1 receptors are held fixed and unlabelled ones are added in increasing
numbers, so purity and retention can be read against a known noise fraction.

`run_benchmark_a()` in `sptce_pipeline.R` produces Figure 9 from it.

## One structural run, six rungs

Only the largest subset is modelled structurally. `score_paired_all()` loads
the pairwise distances once, from the `subset*_clusters/` folder that
`.detect_ref_subset()` finds, and derives every rung by restricting that pair
list to the rung's receptors before clustering:

```r
pairs <- .load_ref_pairs(pr$ref_subset, max(rmsd_thresholds), dih_map, ...)

for (i in ladder$table$subset) {
  uni <- intersect(ladder$subsets[[i]], pr$betas)   # this rung's membership
  ids <- pr$map$id[pr$map$cdr3b %in% uni]
  cl  <- .cluster_from_pairs(pairs, ids, rmsd_cut = thr, ...)
  ...
}
```

The six rungs are therefore six clusterings of one distance table over nested
sets of receptors, and the structures are predicted once. The published methods
are filtered to the same receptors (`uni`) at each rung, so both sides of every
comparison see the same set.

Because each rung is a restriction of the largest, the rungs are not
statistically independent of one another; Figure 9's points should be read as
a series, not as independent samples.

## The rungs

`jouannet_10x_ladder()` reads `$tenx_root/subsets/subset_{1..6}.txt`, which
defaults to the Jouannet repository's `Data/10X/Donor1/subsets/`:

| rung | unique CDR3β | increment | noise `(n − n₁)/n` |
|---|---|---|---|
| 1 | 2,876 | — | 0.0 % |
| 2 | 7,613 | 4,737 | 62.2 % |
| 3 | 12,350 | 4,737 | 76.7 % |
| 4 | 17,087 | 4,737 | 83.2 % |
| 5 | 21,824 | 4,737 | 86.8 % |
| 6 | 26,561 | 4,737 | 89.2 % |

Each rung contains the one below it. `jouannet_10x_ladder()` reports
"arithmetic ladder" when the increments are constant, and
`data/benchmark_a/verify_ladder.sh` checks a set of subset tables for nesting
and distinct rungs.

### Labelled counts

Three labelled counts appear, each correct for its own definition:

- **2,876** — rung 1 of the Jouannet ladder, the annotated CDR3β set.
- **2,374** — the count in the paper: the annotated receptors that also have a
  paired structural model (2,478 of the 2,876 β are paired at rung 1, carrying
  2,374 labels).
- **3,056** — a reconstruction from the dextramer matrix, used by the tables
  under `data/benchmark_a/`; it applies fewer contig filters than
  `10X_processing.R`.

`jouannet_10x_labels()` uses the benchmark's definition by default
(`align_to_subset1 = TRUE`) and prints the difference from the reconstruction.

## The tables under `data/benchmark_a/`

`donor1_subset6_sptce.tsv` and its labels are the Stage A input for the
structural run in `subset4_clusters/`. They are kept because that run was built
from them; `run_benchmark_a()` takes its rungs from the Jouannet files above,
not from these tables. In this construction rungs 4, 5 and 6 contain the same
23,529 receptors:

| rung | total | labelled | noise |
|---|---|---|---|
| 1 | 3,056 | 3,056 | 0 % |
| 2 | 8,089 | 3,056 | 62.2 % |
| 3 | 21,413 | 3,056 | 85.7 % |
| 4–6 | 23,529 | 3,056 | 87.0 % |

`subset4_clusters/` (distance tables and cluster assignments) and
`subset4_dihedrals/` are the reference structural run every rung is derived
from. `.detect_ref_subset()` picks the run with the most distance files, warns
if it has no matching dihedral table, and reports how many of its receptor ids
the dihedral table covers.

## Replicate runs of the published methods

`Data/10X_Output_methods/` holds three replicate runs of each method at each
rung (`subset_1a`, `subset_1b`, `subset_1c` … `subset_6c`). For the
deterministic methods (Hamming, Levenshtein, iSMART, GIANA, TCRMatch) the
replicates are identical; for clusTCR and DeepTCR they differ:

```
DeepTCR subset_2a   7,372 sequences   2,478 clusters
DeepTCR subset_2b   7,200 sequences   2,452 clusters
DeepTCR subset_2c   6,953 sequences   2,582 clusters
```

`run_benchmark_a()` uses replicate `a` by default. To see run-to-run spread for
the stochastic methods:

```r
run_benchmark_a(versions = c("a", "b", "c"))
```

## Running it

```r
source("sptce_pipeline.R")
CFG$tenx_root  <- ".../TCR_Unsupervised_Benchmark/Data/10X/Donor1"
CFG$tenx_out   <- ".../TCR_Unsupervised_Benchmark/Data/10X_Output_methods"
CFG$sptce_root <- "<folder holding subset4_clusters/ and subset4_dihedrals/>"
CFG$binmat_file <- ".../vdj_v1_hs_aggregated_donor1_binarized_matrix.csv"
CFG$contig_file <- ".../vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv"
run_benchmark_a()
```

### Inputs from 10X Genomics

The epitope labels are rebuilt from two 10X Genomics donor-1 tables, which are
not included here:

| `CFG$` | file |
|---|---|
| `binmat_file` | `vdj_v1_hs_aggregated_donor1_binarized_matrix.csv` |
| `contig_file` | `vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv` |

The Jouannet repository's README links them at
<https://zenodo.org/records/6952657>. `jouannet_10x_labels()` first looks for
`$tenx_root/labelled_cdr3b.csv`; otherwise it rebuilds the labels as
`10X_processing.R` does (mean dextramer score per CDR3β across its cells, binder
if > 10).

## Results

The run produces the labelled counts per rung:

```
subset 1 ( 0.0% noise):  2478 of  2876 betas paired (86.2%), 2374 labelled
subset 2 (62.2% noise):  6492 of  7613 betas paired (85.3%), 2374 labelled
subset 3 (76.7% noise): 10500 of 12350 betas paired (85.0%), 2374 labelled
subset 4 (83.2% noise): 14469 of 17087 betas paired (84.7%), 2374 labelled
subset 5 (86.8% noise): 18497 of 21824 betas paired (84.8%), 2374 labelled
subset 6 (89.2% noise): 22560 of 26561 betas paired (84.9%), 2374 labelled
```

Figures are in [`results/figures/benchmark_a/`](../results/figures/benchmark_a/).
At 0.5 Å SPACE-TCR's purity goes from 0.95 at 0% noise to 0.63 at 89%, while
every sequence method falls from ~0.9 to ~0.15 by the second rung. Purity has
to be read alongside yield:

| RMSD | clusters across the ladder | signal recovery |
|---|---|---|
| 0.50 Å | 18–37 | 1.8–2.1 % |
| 0.75 Å | 30–164 | 3.3–5.1 % |
| 1.00 Å | 51–844 | 7.1–17.8 % |

At 0.5 Å purity is high but only about 2% of the signal is recovered, in 18–37
clusters. At 1.0 Å recovery is higher and the purity advantage is largely gone.
0.75 Å keeps both, which supports the focal threshold used in the paper.

## Provenance

10X Genomics donor 1, via the 10X arm of the Jouannet benchmark
(`Data/10X/Donor1` and `Data/10X_Output_methods`; see
[`third_party/JOUANNET.md`](../third_party/JOUANNET.md)). The subset files are
built by `Functions/subset_functions.R` there.
