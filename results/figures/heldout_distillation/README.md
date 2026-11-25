# Held-out distillation figures

These figures score methods on the **1,389 held-out receptors** of the
sequence-only distillation experiment, and each includes a `distilled model`
series. The distilled model appears in the paper only as Supplementary Note S12,
a negative result; these figures support that note.

The held-out subset has a different denominator from the whole cohort, so these
numbers should not be set beside whole-cohort numbers. On the held-out
receptors, SPACE-TCR's sequence-divergent share is 74.6% over 40 high-purity
clusters; on the whole cohort at the same threshold it is 22.7% over 141, the
figure in Section 3.4 and Table 3. Restricting to held-out epitopes removes a
cluster's members from other epitopes, which raises purity and changes the
divergent share.

The paper's Figure 2C is
[`../house_style/fig02c_divergent_share.png`](../house_style/fig02c_divergent_share.png),
built from the whole-cohort table by `scripts/fig_rebuild_supp.R`.

| file | what it shows |
|---|---|
| `fig_divergent_share_heldout.png` | divergent share by method, held-out receptors |
| `fig_method_comparison_heldout.png` | purity/retention by method, held-out receptors |

Source table: `results/tables/heldout/comparison_focal075.csv`.
Code: `ml/compare_all_methods.py`, then `scripts/fig_method_comparison.R`.
