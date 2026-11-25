# Benchmark A inputs — the 10X donor-1 noise ladder

See [`../../docs/BENCHMARK_A.md`](../../docs/BENCHMARK_A.md) for how the
benchmark works.

These tables are the input to the structural run in `subset4_clusters/`, which
every rung of Figure 9 is derived from. `run_benchmark_a()` takes the rungs
themselves from the Jouannet repository's `Data/10X/Donor1/subsets/`
(2,876 to 26,561 CDR3β in steps of 4,737). In the construction stored here,
rungs 4, 5 and 6 contain the same 23,529 receptors.

## Contents

| file | rows | what |
|---|---|---|
| `donor1_subset6_sptce.tsv` | 23,529 | SPACE-TCR input for the last rung: `id`, V/J and CDR3 for both chains |
| `donor1_subset6_labels.tsv` | 23,529 | epitope labels (3,056 labelled) |
| `ladder_membership.tsv` | 23,529 | `id` → first rung the receptor appears in |
| `subset_summary.tsv` | 6 | the ladder: totals, labelled, noise fraction |
| `rebuild_rungs.R` | | writes out the six per-rung tables from the two above |
| `verify_ladder.sh` | | checks a set of rung tables for nesting and distinct rungs |
| `subset4_clusters/cluster_assignments.tsv` | 19,275 | Stage C output for the reference run |
| `subset4_clusters/member_to_bin.tsv`, `bins.jsonl` | | receptor → structural bin, bin definitions |
| `subset4_clusters/distance_tables/` | 1,704 files | pairwise RMSD tables, 15 MB |
| `subset4_dihedrals/dihedrals.tsv.gz` | | Stage D output, 8 MB |

```bash
Rscript data/benchmark_a/rebuild_rungs.R build   # write out the six rungs
bash    data/benchmark_a/verify_ladder.sh build  # and check them
```

Only the reference run has structural output. `score_paired_all()` loads its
pair list once and derives each rung by restricting to that rung's receptors, so
the structures are predicted once.

## Labelled counts

The 3,056 labelled receptors here, the 2,876 of Jouannet's rung 1 and the 2,374
in the paper differ in definition; [`../../docs/BENCHMARK_A.md`](../../docs/BENCHMARK_A.md)
explains each.

## Provenance

10X Genomics donor 1, through the 10X arm of the Jouannet benchmark
(`Data/10X/Donor1`; see [`../../third_party/JOUANNET.md`](../../third_party/JOUANNET.md)).
The subsets are built by `Functions/subset_functions.R` there.
