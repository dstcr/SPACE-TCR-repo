# The Jouannet benchmark repository

`CFG$jouannet_repo` points at this repository. It supplies the nine published
methods' cluster output, the pooled human database the cohort is drawn from, and
the reference metric definitions that `../sptce_metrics.R` follows.

It is not included here, for three reasons:

1. **Licence.** The repository has no licence file, so it cannot be
   redistributed.
2. **Third-party databases.** 292 MB of its 382 MB is seven files that
   redistribute IEDB, VDJdb and McPAS-TCR, each with its own terms:

   ```
   74 MB  Data/Database/IEDB_tcell_receptor_table_MHCI.csv
   73 MB  Data/Database/vdjdb.txt
   72 MB  Data/Database/IEDB_tcell_assay_table_MHCI.csv
   30 MB  Data/Database/database_pooled_human_2023_03_15.txt
   22 MB  .RData
   13 MB  Data/Database/IEDB_tcell_assay_table_MHCII.csv
   8 MB   Data/Database/McPAS-TCR.csv
   ```
3. **Size.** It would make this repository more than ten times larger.

## Getting it

```bash
bash third_party/fetch_jouannet.sh            # clones beside this repository
bash third_party/fetch_jouannet.sh ~/benchmarks
```

Then set `CFG$jouannet_repo` in `../sptce_pipeline.R` to its location. The fetch
script checks that the eight paths the pipeline reads are present and names any
that are missing.

## What this repository reads from it

| path | used for |
|---|---|
| `Data/Output/*` | the nine methods' cluster assignments, parsed by `parse_published_outputs()` |
| `Functions/purity_functions.R` | reference metric definitions, reimplemented in `../sptce_metrics.R` |
| `Functions/similarity_functions.R` | as above |
| `Functions/subset_functions.R` | the Benchmark A noise ladder |
| `Analysis/benchmark_tools_pooled_database_CD8.Rmd` | the metric conventions, reimplemented in `../sptce_metrics.R` |
| `Data/Database/database_pooled_human_2023_03_15.txt` | rebuilding the cohort from the public database |

`../sptce_metrics.R` reimplements the metrics so that SPACE-TCR and the
published methods are measured with the same equations on the same units:

- a cluster counts when `size >= 2`
- high purity is `purity >= 0.9`
- retention counts distinct receptors
- sensitivity counts a cluster when `size > 3` and `n_e/size >= 0.67`, and
  divides the distinct matching receptors by the number of input receptors
  carrying the epitope

## Citing

Jouannet et al. (2025), *Benchmarking unsupervised methods for inferring TCR
specificity*. Cite IEDB, VDJdb and McPAS-TCR as well; the cohort labels come
from them.
