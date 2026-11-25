# results/figures/

`rebuild/` holds the figures as they appear in the paper. The other files in
this folder are drawn by the pipeline's own `fig_*` functions.

To redraw the paper figures:

```bash
Rscript scripts/fig01_rebuild.R      results/figures/rebuild
Rscript scripts/fig_rebuild_all.R    results/figures/rebuild
Rscript scripts/fig_rebuild_supp.R   results/figures/rebuild
python3 scripts/fig_s8_tcrvdb.py     results/figures/rebuild
python3 scripts/fig_s9_constant_alpha.py results/figures/rebuild
```

The two Python scripts need matplotlib and pandas. Figure S8 reads the TCRvdb
summary tables (the TCRvdb labels themselves are request-only). Figure S9
reads the cohort's distance and dihedral tables for its paired arm, and for its
single-chain arm the constant-alpha run in
`results/tables/extensions/s9_constant_alpha/`: the bin's 1,523 receptors
re-modelled with one fixed alpha partner. `scripts/s9_constant_alpha.sh`
rebuilds that run from the Zenodo structures.

`make_figures.R` also runs all of these, last, into `<out_dir>/rebuild`. They read
only the tables in this repository, so they need no structure prediction or
pipeline run, and each checks its headline values against the numbers in the
paper before writing a file.
