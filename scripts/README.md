# scripts/

## The paper figures

These three draw the figures in the paper, in one shared style
(`fig_common.R`). They read only the tables in this repository, so they need no
structure prediction or pipeline run, and each checks its headline values
against the numbers in the paper before writing a file. `make_figures.R` runs
all three, last, into `<out_dir>/rebuild`.

| Script | Produces |
|---|---|
| `fig_common.R` | shared theme, type scale, palette and checking helper, sourced by the others |
| `fig01_rebuild.R` | Figure 1 |
| `fig_rebuild_all.R` | Figures 2 to 6 |
| `fig_rebuild_supp.R` | Supplementary Figures S1 to S9 |

## Analysis scripts

`corrected_sensitivity.R`, `corrected_sensitivity_figures.R`,
`fig_method_comparison.R` and `common_denominator_retention.R` produce tables
and checks used in the paper; they are not the paper's figure scripts. Some
write files with the same names as pipeline output
(`corrected_sensitivity_figures.R` writes `fig03_sensitivity_heatmap.png`), so
run them into a separate folder rather than over `results/figures/`.

## Pipeline stages

`stage1_stitchr.py`, `stage2_tcrbuilder2.py`, `make_bins.py`,
`cluster_one_bin.py`, `extract_dihedrals.py` and `run_clustcr*.py` are the
Stage A–E programs that the job scripts in `../sbatch/` call.
