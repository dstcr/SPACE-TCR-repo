# Tests

None of the tests need data; each builds its own fixture.

```bash
Rscript tests/test_metrics.R           # 41 unit tests on sptce_metrics.R
Rscript tests/test_e2e.R               # 17 checks, from files on disk to scores
Rscript tests/test_dihedral_distance.R # 13 checks on the dihedral distance, R and Python
Rscript tests/test_pipeline.R          # 44 checks on sptce_pipeline.R
PYTHONPATH=. python3 tests/verify_space2tcr.py .   # 22 checks on the clustering library
```

To run them all, with the toolchain checks:

```bash
bash check_install.sh
```

The R tests need `dplyr readr tidyr purrr stringr tibble stringdist`.

## `test_metrics.R`

Checks each metric definition against hand-computed values:

- global purity is `sum(n_dom)/sum(size)` over clusters of size ≥ 2
- high purity is `>= 0.9` (a cluster at exactly 0.9 counts)
- retention counts distinct receptors, not summed cluster sizes
- a cluster qualifies for sensitivity when `size > 3` and `n_e/size >= 0.67`,
  so a 4-of-6 cluster (0.6667) does not
- the sensitivity numerator counts distinct receptors, so a receptor in two
  qualifying clusters counts once
- a unit with two epitope labels counts towards each, as in the reference
- a published method is scored on its own cluster sizes, and clusters with no
  annotated member are left out, as in `purity_functions.R`
- the circular dihedral distance wraps (0.1 and 2π−0.1 are 0.2 apart)

It also checks that `check_metrics` and `check_sensitivity` reject out-of-range
values.

## `test_e2e.R`

Writes a 60-receptor cohort (3 epitopes, 2 CDR3-length bins) to disk in the
formats Stages C and D produce, with small RMSDs and similar dihedrals within
an epitope. It clusters and scores the cohort with the same `sptce_metrics.R`
calls the pipeline uses, and checks that:

- all five thresholds are swept, with a gated mode
- retention and purity stay within range, and the clustered count never exceeds
  the denominator
- retention rises with the threshold, and purity is about 1 at a tight threshold
- the gate never clusters more than RMSD alone, and a tighter gate never
  clusters more than a looser one
- sensitivity never exceeds the clustered fraction
- gene homogeneity and published-method scoring run on the same clusters

## `test_dihedral_distance.R`

Checks `circ_dihedral()` against hand-computed values and against
`scripts/dihedral_util.py`, the Python version Stages C and D use.

## `test_pipeline.R`

Loads `sptce_pipeline.R` and checks that:

- every stage function and the figure helpers they call are defined
- `run_stage5d_cohort_accounting` counts each operating point separately
- the HuARdb functions are defined but not called from `main()`
- `CFG$focal_dihedral` is set, is one of `CFG$dihedral_modes`, and is read by
  the code
- the Benchmark A methods figure draws no gate by default, and the thresholds
  figure takes its gates (0.3 and 0.5 rad) from a named constant
- no gate value is written into the code as a literal
- the noise ladder has six rungs, and the thresholds and focal threshold match
  the paper
