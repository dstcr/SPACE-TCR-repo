# data/

The inputs for both benchmarks and the Stage C/D geometry they were clustered
into. With these, the paper's numbers reproduce without predicting any
structures:

```bash
Rscript validate_outputs.R        # under a minute, 36 checks, exit status = failures
```

`validate_outputs.R` reads the two input tables and the two geometry folders
below, recomputes the headline numbers and compares them with
[`../results/`](../results/README.md): at the focal operating point, 149
clusters, 681 receptors clustered, 94.63% of clusters at purity ≥ 0.9 and
pooled purity 0.9868 (98.7% in the paper), plus the six sensitivity
denominators.

The predicted structures themselves (the TCRBuilder2 PDBs) are on Zenodo:
<https://zenodo.org/records/23765886>. They are needed only to rerun Stages C
and D or to inspect individual models.

---

## Benchmark B — the curated CD8 cohort

### `jouannet_labelled_input.tsv` — 535 KB, the Stage A input

**4,779 paired αβ receptors**, one per row, with epitope labels. This is the
table every figure and table in [`../results/`](../results/README.md) was
produced from. `CFG$jouannet_labelled` in `sptce_pipeline.R` defaults to it.

```
id  v_call_alpha  j_call_alpha  cdr3_aa_alpha
    v_call_beta   j_call_beta   cdr3_aa_beta
    epitope  antigen  organism
```

`id` is `tcr1 … tcr4779` and is the join key against the cluster assignments,
the distance tables and the dihedral table. 4,292 receptors carry a unique
CDR3β (the figure in §3.9 of the paper); there are 367 distinct epitopes.

`config.sh` points `INPUT_CSV` at this file for Stage A.

**Provenance.** Derived from the pooled human database of Jouannet et al.
(2025), itself pooled from IEDB, VDJdb and McPAS-TCR. The epitope labels come
from those resources and carry their terms of use; cite them alongside this
repository. The cohort can be rebuilt from the public database (4,758 unique
CDR3β, 454 for GILGFVFTL).

### `cohort_structural/` — 5.9 MB, the Stage C and D output

```
cohort_structural/
  cluster_assignments.tsv              4,031 rows, the 1.25 A assignment
  clusters/member_to_bin.tsv           receptor -> structural bin
  dihedrals.tsv.gz                     CDR3 backbone phi/psi, both chains
  distance_tables/bin_*.pairwise.tsv.gz    828 files, the pairwise geometry
```

Every RMSD threshold up to 1.25 Å, including the focal 0.75 Å, can be derived
from the pairwise tables. Columns are `i j id_i id_j status rmsd dihedral_alpha
dihedral_beta clustering_dist error`; `id_i` and `id_j` are the PDB paths on the
machine that produced them, and the receptor id is the file stem. 9,642 pairs
fall within 1.25 Å, over 1,555 receptors; the other 3,224 receptors have no
neighbour within 1.25 Å and are singletons at every threshold.

Bins with no pairs within the cut-off are header-only files, so read the tables
with explicit column types.

---

## Benchmark A — the 10X donor-1 noise ladder

See [`../docs/BENCHMARK_A.md`](../docs/BENCHMARK_A.md) and
[`benchmark_a/README.md`](benchmark_a/README.md) for how the ladder is built.

```
benchmark_a/
  donor1_subset6_sptce.tsv         Stage A input for the last rung: 23,529 receptors
  donor1_subset6_labels.tsv        epitope labels, 3,056 labelled
  ladder_membership.tsv            id -> first rung the receptor appears in
  rebuild_rungs.R                  writes out all six rungs from those two
  subset_summary.tsv               the ladder as specified
  subset4_clusters/
    cluster_assignments.tsv        Stage C output for the reference rung
    member_to_bin.tsv  bins.jsonl
    distance_tables/               1,704 pairwise tables, 15 MB
  subset4_dihedrals/dihedrals.tsv.gz    Stage D output, 8.2 MB
  verify_ladder.sh
```

The ladder is nested, so each rung is the largest one restricted by receptor
id, and the epitope labels are the same throughout. Only the last rung is
stored; rung N is every receptor whose `first_rung` is at most N:

| rung | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| receptors | 3,056 | 8,089 | 21,413 | 23,529 | 23,529 | 23,529 |

`validate_outputs.R` checks these six sizes and the nesting on every run, and
`.rung_table()` in `sptce_pipeline.R` derives the rungs in memory. A folder
holding six per-rung files is also read as-is.

```bash
Rscript data/benchmark_a/rebuild_rungs.R build   # write out all six rungs
bash    data/benchmark_a/verify_ladder.sh build  # and check them
```

Rungs 4, 5 and 6 are the same 23,529 receptors, so one structural run serves
all six: `score_paired_all()` loads the reference pair list once and restricts
it per rung. That is why only `subset4_clusters/` and `subset4_dihedrals/`
exist.

The Jouannet ladder (`Data/10X/Donor1/subsets/subset_{1..6}.txt`, 2,876 →
26,561 in steps of 4,737) lists CDR3β only, so it cannot be a SPACE-TCR input
directly. The paired tables here were built by joining that list back to the
10X donor-1 contig annotations.

---

## Not in this repository

| needed for | what | where |
|---|---|---|
| inspecting or re-clustering the structures | the modelled PDBs | [Zenodo record 23765886](https://zenodo.org/records/23765886) |
| the nine published methods' own cluster output | the Jouannet benchmark repository | [`../third_party/JOUANNET.md`](../third_party/JOUANNET.md) |
| the ladder's noise fraction from raw counts | `vdj_v1_hs_aggregated_donor1_{binarized_matrix,all_contig_annotations}.csv` | 10X Genomics |
| Figure 8 (training-set leakage) | `igblast_tcr_parsed.csv` | derived from STCRDab |

The last three are read from `<SPTCE_EXTERNAL>/<filename>`, where
`SPTCE_EXTERNAL` is an environment variable defaulting to `~/sptce_external`:

```bash
export SPTCE_EXTERNAL=~/sptce_external
bash third_party/fetch_jouannet.sh      # clones the benchmark repository into it
```
