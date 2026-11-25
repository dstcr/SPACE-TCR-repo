# Distilling the structural clustering into a sequence-only model

## The idea

SPACE-TCR needs a predicted structure for every receptor, and structure
prediction limits how many receptors it can process (§4 of the paper).

`distill.py` trains a siamese encoder on **sequence alone**, supervised by the
structural clusters SPACE-TCR produced. If it generalises, structural
co-membership could be predicted for receptors that were never modelled.

```bash
python3 ml/distill.py --synthetic                      # self-test, no data
python3 ml/distill.py --assignments cluster_assignments.tsv \
                      --labels data/jouannet_labelled_input.tsv \
                      --out ml_out
```

It needs only the cluster assignments and the labelled input; the distance
tables are not required.

## Baseline

Every number is reported next to **Levenshtein distance on CDR3β at a comparable
operating point**. A model that beats chance but not edit distance has learned
edit distance. The script's closing summary states which way the comparison
went.

## What it measures

1. **Pair-level ranking** on held-out epitopes: AUROC and average precision for
   "did SPACE-TCR group these two?", model versus Levenshtein.
2. **Cluster-level metrics**, after clustering the learned embedding by connected
   components of a distance-threshold graph (the procedure SPACE-TCR uses),
   swept across thresholds as Table 1 sweeps RMSD, and scored by
   `ml/metrics.py`.
3. **Recovery of sequence-diverse clusters**: positives whose CDR3β Levenshtein
   distance is ≥ 4, the 22.7% of clusters the paper describes as out of reach of
   sequence methods. This is the comparison the distillation question depends
   on.

It also reports ARI against the structural clustering for both methods.

## Splits

Epitopes are held out whole (`epitope_split`), so the test is on unseen
specificities; a receptor-level split would let near-clonal receptors of the
same epitope fall on both sides. Test receptors whose CDR3β appears in training
are also dropped, and the number dropped is reported.

## Agreement with the R metrics

`ml/metrics.py` is a Python copy of `sptce_metrics.R` (training runs in Python),
so a learned model can be scored on the same terms as SPACE-TCR. The two are
checked against each other:

```bash
python3 ml/tests/test_metrics_match_r.py
```

It generates random assignment tables, with duplicate receptors, singletons,
unlabelled rows and repeated rows, scores each in both languages and compares
the aggregates and the per-cluster sizes and purities (284 agreements,
0 disagreements).

## The synthetic fixture

`--synthetic` builds a cohort modelled on §3.6 of the paper: clusters coherent
in CDR3α and the V/J genes while diverse in CDR3β. In most clusters CDR3β is
drawn at random, so CDR3β edit distance is at chance by construction; a
minority also share a β motif, covering the cases edit distance can reach.
Simpler designs (one CDR3 length per cluster, fixed shared residues) let edit
distance identify the clusters and so could not test the model.

On this fixture the model reaches **0.959 AUROC against 0.730** for Levenshtein,
and ARI **0.520 against 0.020**.

Two points from the synthetic run carry over to the real data:

- **A good pairwise ranker is not a good clustering.** At 0.959 pair AUROC the
  model still does not match the structural purity at any threshold, because
  connected components chain through borderline pairs. Average linkage or a
  density method would be less prone to this.
- **The cohort is small for metric learning.** 4,779 receptors over 368
  epitopes, with the top three holding 969. Holding out 30% of epitopes leaves
  little to train on, so results vary between seeds.

## Results on the curated cohort

Run at three supervision levels. `cluster_assignments.tsv` is the 1.25 Å
clustering (321 clusters over 1,555 receptors, 78.4% purity); the 0.75 Å and
1.00 Å assignments were rebuilt from `run_stage5()`'s `all_sptce` table.

```
supervision  train  test  held-out epi |  model     LD  |  model(div)  LD(div) |  mARI  lARI
0.75 A         528   152       8       |  0.910  0.906  |     0.784    0.757   | 0.000  0.339
1.00 A         850   148      21       |  0.876  0.915  |     0.866    0.810   | 0.195  0.655
1.25 A         874   674      53       |  0.684  0.647  |     0.599    0.532   | 0.936  0.938
```

`model` and `LD` are pair-level AUROC on held-out epitopes; `(div)` restricts to
positives with CDR3β Levenshtein distance ≥ 4.

### Five seeds at 1.25 Å

The three supervision levels share a cohort and are nested, so they are not
replicates. Five seeds at 1.25 Å are (`--seed 0..4`, 53 held-out epitopes each):

```
seed  bestEp  test | model     LD    diff | model(div) LD(div)   diff |  mARI  lARI
   0       2   674 | 0.684  0.647  +0.038 |     0.599   0.532  +0.067 | 0.936 0.938
   1       1   310 | 0.711  0.750  -0.040 |     0.654   0.679  -0.025 | 0.877 0.781
   2       1   405 | 0.839  0.829  +0.010 |     0.739   0.714  +0.025 | 0.896 0.891
   3       1   204 | 0.788  0.747  +0.040 |     0.739   0.669  +0.070 | 0.592 0.600
   4       1   474 | 0.887  0.861  +0.026 |     0.801   0.735  +0.066 | 0.968 0.967
```

Paired across seeds:

| | mean | sd | t(4) | p | 95% CI |
|---|---|---|---|---|---|
| **divergent-pair advantage** | **+0.041** | 0.041 | 2.21 | **0.092** | **[−0.010, +0.092]** |
| overall advantage | +0.015 | 0.033 | 1.01 | 0.371 | [−0.026, +0.056] |

The advantage over edit distance is not statistically significant: both
confidence intervals include zero, and the model is behind on seed 1. The
direction is consistent (4 of 5 seeds on both measures), and the effect is
larger on sequence-divergent pairs than overall (+0.041 against +0.015). Test
sets range from 204 to 674 receptors depending on which epitopes are held out,
and overall AUROC from 0.684 to 0.887, so single-seed numbers vary widely.

### Limitations

1. **The clustering does not reach the structural purity.** `reached_truth_purity`
   is `False` at all three supervision levels. In these runs the encoder had
   collapsed (maximum pairwise distance 0.065 on training receptors, against a
   threshold sweep starting at 0.05), so every threshold merged everything.
   `distill.py` applies a `MIN_SPREAD` floor in model selection to prevent this;
   the numbers above predate it.
2. **The tighter runs are small.** 8 and 21 held-out epitopes, 152 and 148 test
   receptors, so the 0.910 at 0.75 Å is not comparable with the 0.684 at 1.25 Å.
3. **The selected epoch is early.** Across the five seeds the restored epoch was
   1, 1, 1, 1 and 2 of 15–40, while training loss keeps falling. Without early
   stopping the 1.25 Å run scores 0.656 against Levenshtein's 0.647; with it,
   0.684.

### Summary

Structural co-membership among sequence-divergent receptors appears partly
predictable from sequence (+0.041 ± 0.041 AUROC over edit distance across five
seeds), but the effect is not significant and the encoder does not yet produce a
usable clustering. This is reported as a negative result in Supplementary Note
S12. A larger cohort is the main requirement.

### Reproducing

```bash
# rebuild assignments at a chosen threshold from a saved stage-5 object
Rscript -e 'source("sptce_pipeline.R"); s5 <- run_stage5()
            x <- subset(s5$all_sptce, mode=="SPACE-TCR RMSD" & threshold==0.75)
            write.table(data.frame(tcr_id=x$tcr_id, cluster_id=paste0("c",x$cluster_id)),
                        "assign_0.75.tsv", sep="\t", row.names=FALSE, quote=FALSE)'

python3 ml/distill.py --assignments assign_0.75.tsv \
                      --labels data/jouannet_labelled_input.tsv \
                      --epochs 40 --out ml_out_0.75
```

## Files

| | |
|---|---|
| `distill.py` | data loading, model, training, evaluation, command line |
| `metrics.py` | the Jouannet metrics in Python |
| `summarise_runs.py` | reduces a set of runs to the paired test above |
| `tests/test_metrics_match_r.py` | checks `metrics.py` against `sptce_metrics.R` |
| `results/*.json` | the runs reported above |

```bash
python3 ml/summarise_runs.py ml/results/{125es,seed1,seed2,seed3,seed4}.json
```

reproduces the seed table and the t-test. Pass only true replicates; the three
supervision levels are not.

Requires `torch`, `numpy`, `pandas` and `rapidfuzz`. A CPU is enough; the
synthetic run trains in under a minute.
