# Rendered results

60 figures and 33 tables, produced by

```bash
Rscript make_figures.R
```

at `CFG$focal_threshold = 0.75` Å and `CFG$focal_dihedral = 0.5` rad, from the
4,779-receptor curated CD8 cohort.

## Benchmark A and training-set leakage

| | |
|---|---|
| `benchmark_a/fig_paired_methods.png` | Figure 9 — purity and cluster size across the six-rung noise ladder, SPACE-TCR against all nine published methods at three RMSD thresholds |
| `benchmark_a/fig_paired_thresholds.png` | the same cohort across all five RMSD thresholds and the dihedral gate |
| `stage7/fig31_stcrdab_overlap.png` | Figure 8 — TCRBuilder2 training-set leakage |

Figure 8: **45 of 4,292 cohort CDR3β match a training PDB exactly, 1.05%.**

Figure 9, with the yield alongside purity:

| RMSD | purity, 0% → 89% noise | clusters | signal recovery |
|---|---|---|---|
| 0.50 Å | 0.95 → 0.63 | 18–37 | 1.8–2.1 % |
| 0.75 Å | 0.96 → 0.31 | 30–164 | 3.3–5.1 % |
| 1.00 Å | 0.88 → 0.15 | 51–844 | 7.1–17.8 % |

Every sequence method falls from ~0.9 to ~0.15 by the second rung (62% noise).
At 0.5 Å SPACE-TCR's purity is well above all of them, but it forms only 18
clusters; at 1.0 Å it recovers more and the purity advantage has largely gone.
0.75 Å keeps both, which supports the focal threshold.

Reproduce with:

```bash
Rscript run_benchmark_a.R <jouannet_repo> <dir with subsetN_clusters/>
```

See [`../docs/BENCHMARK_A.md`](../docs/BENCHMARK_A.md) for the two 10X inputs it needs.

## Sensitivity

Sensitivity follows Eq. 4: distinct qualifying units carrying an epitope,
divided by the **input** units carrying it. Figures are in
`figures/corrected/` and tables in `tables/corrected/`.

Mean over the six focal epitopes, top three:

| method | mean sensitivity | rank |
|---|---|---|
| DeepTCR | 0.340 | 1 |
| SPACE-TCR RMSD+dih (1 Å) | 0.307 | 2 |
| SPACE-TCR RMSD (0.75 Å) | 0.265 | 3 |

Across all 4,037 method × epitope rows, sensitivity is at most the fraction of
that epitope's receptors that were clustered; `check_sensitivity()` in
`sptce_metrics.R` checks this.

Small epitopes give unstable values. At LLFGYPVYV SPACE-TCR recovers 0.90 in
both modes while seven of the nine sequence methods score 0.00, but the epitope
has 10 input receptors. Three of the six focal epitopes (LLFGYPVYV, NYNYLYRLF,
FLCMKALLL) have fewer than 90 input receptors; GILGFVFTL is the largest (705
paired / 455 β / 938 chains).

### Units

Five of the nine published methods were run on a chain-pooled single-chain list,
so their clustering unit is **one chain**; four were run paired; SPACE-TCR is
paired. Sensitivity and retention are therefore divided by different pools —
938 chains, 455 β keys or 705 paired receptors for GILGFVFTL.
`figures/corrected/fig_unit_mismatch.png` shows this. Sensitivity is computed
per unit here; retention for the five chain-pooled methods is reported in the
unit their runs used.

### Regenerating

From the repository root:

```bash
Rscript -e 'source("sptce_pipeline.R"); saveRDS(run_stage5(), "stage5.rds")'
Rscript scripts/corrected_sensitivity.R data/jouannet_labelled_input.tsv stage5.rds results/tables/corrected
Rscript scripts/corrected_sensitivity_figures.R
```

The first step takes about 6 minutes; the other two are quick.

## What reproduces

`sptce_all_threshold_metrics.csv` reproduces all eight of Table 1's cluster
counts from the pairwise RMSDs:

```
0.25 RMSD      21 /   65      0.75 RMSD+dih  139 /  658
0.50 RMSD      95 /  393      1.00 RMSD+dih  172 /  903
0.75 RMSD     149 /  681      1.25 RMSD+dih  221 / 1186
1.00 RMSD     202 / 1000
1.25 RMSD     321 / 1555
```

The gated rows are at a dihedral gate of 0.5 rad.

## Layout

```
figures/               54 PNG, 300 dpi
figures/corrected/      5 PNG            sensitivity (Eq. 4) and units
tables/                15 CSV + 15 under extensions/
tables/corrected/       sensitivity tables + sens_corrected.rds
```

`figure_manifest_r0.75_d0.5.csv` lists every stage, its status and its runtime;
`figure_files_r0.75_d0.5.csv` lists every file written with its size. Both
filenames include the threshold and gate, so runs at different settings are
kept apart.

## Inputs

The inputs are named in `CFG` at the top of `../sptce_pipeline.R`:

| input | what | in this repo? |
|---|---|---|
| `jouannet_labelled` | the 4,779-row labelled input TSV | yes, `data/jouannet_labelled_input.tsv` |
| `jouannet_sptce_out` | `distance_tables/` and `dihedrals.tsv.gz` | yes, `data/cohort_structural/` |
| `jouannet_repo` | the Jouannet benchmark repository (~382 MB) | no — `bash third_party/fetch_jouannet.sh` |
| `stcrdab_parsed_csv` | `igblast_tcr_parsed.csv`, for Figure 8 only | no |
