# =============================================================================
# fig_pairing_summary.R — paired against single-chain structural clustering,
#                         on the metrics used in the benchmark literature
#
#   Rscript scripts/fig_pairing_summary.R [out_dir]
#
# Four panels, all from tables in this repository; the structures are on Zenodo
# (record 23765886).
#
# Jouannet et al. and the methods they benchmark report pooled purity alongside
# the share of clusters whose own purity is >= 0.9 and the share of clustered
# receptors in those clusters. Pooled purity is a weighted average that a few
# large clean clusters can dominate, so panel A shows both, computed with the
# same equations as the pipeline (threshold >= 0.9) for the structural arms and
# the nine published methods alike.
#
# A  Both metrics against unique CDR3-beta clustered over the 4,292-sequence
#    universe, with the nine published methods on the same axes.
#
# B  Pooled purity under repertoire noise, 0% / 62% / 85.7%.
#
# C  Polyspecificity as error counts rather than as a proportion near its
#    ceiling, which shows it agrees with purity.
#
# D  Alpha loops do not vary less than beta loops, so unequal spread does not
#    drive the single-chain comparison.
#
# Every panel uses matched binning (see fig_chain_framework.R).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/pairing_summary"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))
suppressMessages(library(patchwork))

T_ <- function(f) file.path(repo, "results/tables/extensions", f)
pub  <- read_csv(file.path(repo, "results/tables/published_method_metrics.csv"), show_col_types = FALSE)
comm <- read_csv(file.path(repo, "results/tables/retention_common_denominator.csv"), show_col_types = FALSE)
abl  <- read_csv(T_("chain_ablation_matched_binning.csv"), show_col_types = FALSE)
r2   <- read_csv(T_("ladder_rung2_matched.csv"), show_col_types = FALSE)
r3   <- read_csv(T_("ladder_rung3_matched.csv"), show_col_types = FALSE)
dist <- read_csv(T_("chain_rmsd_distribution.csv"), show_col_types = FALSE)
stopifnot("perc_clusters_gt09" %in% names(abl))

ARM <- c(paired_fw = "Paired αβ",
         alpha_fw  = "CDR3α only",
         beta_fw   = "CDR3β only")
PAL <- setNames(c(SPTCE_RED, "#2166ac", "#7fa8cf"), unname(ARM))
tidy <- function(d) d %>% filter(arm %in% names(ARM)) %>%
  mutate(lab = factor(recode(arm, !!!ARM), levels = unname(ARM)))
abl <- tidy(abl); r2 <- tidy(r2); r3 <- tidy(r3)

# ---- assertions ------------------------------------------------------------
expect_close("paired purity at 0.75 A (Table 1)",
             abl$purity[abl$arm == "paired_fw" & abl$threshold == 0.75], 0.9868, 0.001)
expect_close("ClusTCR purity (Table 2)", pub$purity[pub$method == "clusTCR"], 0.993, 0.001)
expect_close("ClusTCR %% clusters >= 0.9 (Table 2)",
             pub$perc_clusters_gt09[pub$method == "clusTCR"], 0.973, 0.001)
expect_close("beta polyspecificity at 0.75 A",
             abl$polyspecific[abl$arm == "beta_fw" & abl$threshold == 0.75], 0.2987, 0.002)

# The two high-purity fractions are computed for the structural arms by
# scripts/chain_framework_ablation.py, in Python, and for the published methods
# by compute_metrics() in sptce_pipeline.R, in R. The paired arm is the one
# point where both implementations score the same clustering, so it is the
# parity check between them. Agreement is exact to six decimals, which is why
# the tolerance here is 1e-6 rather than the 1e-3 used above.
pipe <- read_csv(file.path(repo, "results/tables/sptce_all_threshold_metrics.csv"),
                 show_col_types = FALSE) %>%
  filter(mode == "SPACE-TCR RMSD", threshold == 0.75)
stopifnot(nrow(pipe) == 1)
foc <- abl %>% filter(arm == "paired_fw", threshold == 0.75)
expect_close("R/Python parity, %% clusters >= 0.9",
             foc$perc_clusters_gt09, pipe$perc_clusters_gt09, 1e-6)
expect_close("R/Python parity, %% TCRs >= 0.9",
             foc$perc_tcrs_gt09, pipe$perc_tcrs_gt09, 1e-6)

cat("\n  at the focal 0.75 Å, matched binning:\n")
cat(sprintf("    %-14s %8s %8s %10s %10s\n",
            "arm", "purity", "retent.", "%clus>=.9", "%tcrs>=.9"))
for (a in names(ARM)) {
  r <- abl %>% filter(arm == a, threshold == 0.75)
  cat(sprintf("    %-14s %8.4f %8.4f %10.4f %10.4f\n",
              ARM[[a]], r$purity, r$retention_common,
              r$perc_clusters_gt09, r$perc_tcrs_gt09))
}

# =============================== A ==========================================
MET <- c(purity             = "Pooled cluster purity",
         perc_clusters_gt09 = "% of clusters with purity ≥ 0.9")

pubx <- pub %>% mutate(method = str_replace_all(method, "_", " ")) %>%
  left_join(comm %>% mutate(method = str_replace_all(method, "_", " ")) %>%
              select(method, retention_common), by = "method") %>%
  filter(!is.na(retention_common)) %>%
  select(method, retention_common, purity, perc_clusters_gt09) %>%
  pivot_longer(c(purity, perc_clusters_gt09),
               names_to = "metric", values_to = "y") %>%
  mutate(metric = factor(recode(metric, !!!MET), levels = unname(MET)))
stopifnot(nrow(pubx) == 18)

arms_long <- abl %>%
  select(lab, threshold, retention_common, purity, perc_clusters_gt09) %>%
  pivot_longer(c(purity, perc_clusters_gt09),
               names_to = "metric", values_to = "y") %>%
  mutate(metric = factor(recode(metric, !!!MET), levels = unname(MET)))

# DeepTCR clusters 92% of the universe. Drawn to scale it pushes every other
# method into the left quarter of the panel, so it is carried as an off-scale
# note at the axis edge and the axis stops at 55%.
XMAX <- 0.55
deep <- pubx %>% filter(method == "DeepTCR")
stopifnot(nrow(deep) == 2, all(deep$retention_common > XMAX))
pubx_in <- pubx %>% filter(method != "DeepTCR")

pA <- ggplot() +
  geom_line(data = arms_long, aes(retention_common, y, colour = lab), linewidth = 1.1) +
  geom_point(data = arms_long, aes(retention_common, y, colour = lab), size = 1.9) +
  geom_point(data = pubx_in, aes(retention_common, y), shape = 21,
             fill = "white", colour = "grey25", stroke = 1.0, size = 3.0) +
  geom_text_repel(data = pubx_in, aes(retention_common, y, label = method),
                  size = 3.4 * FS, colour = "grey20", seed = 7,
                  box.padding = 0.45, point.padding = 0.3,
                  min.segment.length = 0.2, segment.colour = "grey70",
                  segment.size = 0.3, max.overlaps = Inf) +
  geom_point(data = arms_long %>% filter(threshold == 0.75),
             aes(retention_common, y, colour = lab), size = 4.2, shape = 18) +
  # Parked in the empty top-right corner of each facet rather than beside its
  # own y value, which lands in the crowded middle of the left facet.
  geom_text(data = deep, aes(x = XMAX - 0.005, y = 1.02,
                             label = sprintf("DeepTCR off scale →\n92%%, %.0f%%",
                                             100 * y)),
            hjust = 1, vjust = 1, lineheight = 0.95,
            size = 3.3 * FS, colour = "grey45", inherit.aes = FALSE) +
  facet_wrap(~metric, nrow = 1) +
  scale_colour_manual(values = PAL, name = NULL,
                      guide = guide_legend(nrow = 1, override.aes =
                                             list(linewidth = 1.6, size = 0))) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, XMAX),
                     breaks = seq(0, 0.5, 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.03),
                     breaks = seq(0, 1, 0.2)) +
  labs(x = "Unique CDR3β clustered, of the cohort's 4,292", y = NULL,
       title = "A  One chain's geometry already beats the sequence methods",
       subtitle = paste0("Both metrics Jouannet et al. report. Pooled purity ",
                         "is receptor-weighted, so a few large clean\n",
                         "clusters can carry it; the 0.9 fraction asks how ",
                         "many clusters a user would trust.\nArms sweep ",
                         "0.4–1.25 Å, diamonds the focal 0.75 Å. The ",
                         "right facet separates α from β.")) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.text = element_text(size = 12.5 * FS, face = "bold"),
        legend.position = "bottom", legend.text = element_text(size = 12 * FS),
        plot.title = element_text(size = 15 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 8)))

# =============================== B ==========================================
noise <- bind_rows(
  abl %>% mutate(cohort = "Clean cohort"),
  r2  %>% mutate(cohort = "62% noise"),
  r3  %>% mutate(cohort = "85.7% noise")) %>%
  mutate(cohort = factor(cohort, c("Clean cohort", "62% noise", "85.7% noise")))

pB <- ggplot(noise, aes(retention, purity, colour = lab)) +
  geom_line(linewidth = 1.1) + geom_point(size = 1.9) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.78),
                     breaks = seq(0, 0.75, 0.25)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.5, 1.02)) +
  labs(x = "Retention", y = "Pooled cluster purity",
       title = "B  The same under repertoire noise",
       subtitle = paste0("Pairing's margin over the better single chain is the ",
                         "same order at 85.7% noise as when clean — noise ",
                         "does not\nwiden it. Absolute purity does fall, which ",
                         "is what Benchmark A reports.")) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.text = element_text(size = 12.5 * FS, face = "bold"),
        plot.title = element_text(size = 15 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 8)))

# =============================== C ==========================================
# Counted at the focal threshold.
tab <- tribble(
  ~arm,            ~clusters, ~receptors, ~npoly, ~bad,
  "Paired αβ",  149,  681,  8,  9,
  "CDR3β only",      154,  822, 46, 58,
  "CDR3α only",      212,  973, 59, 94) %>%
  mutate(poly = npoly / clusters, impurity = bad / receptors,
         poly_r = poly / poly[1], imp_r = impurity / impurity[1])
expect_close("beta polyspecificity ratio", tab$poly_r[2], 5.6, 0.15)
expect_close("beta impurity ratio",        tab$imp_r[2],  5.3, 0.15)

cmp <- tab %>% filter(arm != "Paired αβ") %>%
  select(arm, `Polyspecific clusters` = poly_r, `Mislabelled receptors` = imp_r) %>%
  pivot_longer(-arm, names_to = "metric", values_to = "ratio") %>%
  mutate(arm = factor(arm, c("CDR3α only", "CDR3β only")),
         metric = factor(metric, c("Polyspecific clusters",
                                   "Mislabelled receptors")))

pC <- ggplot(cmp, aes(ratio, arm, fill = metric)) +
  geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.7) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("×%.1f", ratio)),
            position = position_dodge(width = 0.7), hjust = -0.18,
            size = 3.9 * FS, colour = "grey20") +
  annotate("text", x = 1.06, y = 2.52, label = "paired = 1", hjust = 0,
           size = 3.5 * FS, colour = "grey40") +
  scale_fill_manual(values = c("Polyspecific clusters" = SPTCE_RED,
                               "Mislabelled receptors" = "#6e6e6e"), name = NULL,
                    guide = guide_legend(nrow = 2)) +
  scale_x_continuous(limits = c(0, 9), breaks = 0:8) +
  coord_cartesian(clip = "off") +
  labs(x = "Times worse than paired", y = NULL,
       title = "C  Errors counted, not rated",
       subtitle = paste0("Purity 0.987→0.929 looks flat only\n",
                         "because it is a rate near its ceiling.\n",
                         "As counts, β agrees with polyspecificity\n",
                         "(×5.6, ×5.3). α does not: its impure\n",
                         "clusters carry more wrong receptors each.")) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 12.5 * FS, face = "bold", colour = "grey15"),
        legend.position = "bottom", legend.justification = "left",
        plot.margin = margin(8, 26, 8, 8),
        plot.title = element_text(size = 15 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 8)))

# =============================== D ==========================================
D_LV <- c("CDR3α only", "CDR3β only", "Paired (87 residues)")
db <- dist %>% mutate(lab = factor(recode(selection, alpha = D_LV[1],
                                          beta = D_LV[2], paired = D_LV[3]),
                                   levels = D_LV))
medb <- db %>% group_by(lab) %>% summarise(m = median(rmsd), .groups = "drop")
expect_close("median within-bin RMSD, alpha", medb$m[medb$lab == D_LV[1]], 1.848, 0.01)

pD <- ggplot(db, aes(rmsd, lab, fill = lab)) +
  geom_violin(colour = NA, alpha = 0.55, scale = "width", trim = TRUE) +
  geom_point(data = medb, aes(m, lab), size = 3, colour = "grey15", inherit.aes = FALSE) +
  geom_text(data = medb, aes(m, lab, label = sprintf("%.3f Å", m)),
            vjust = -1.4, size = 3.6 * FS, colour = "grey20", inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(c("#2166ac", "#7fa8cf", OTHER_GREY), D_LV),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 8), breaks = 0:8) +
  labs(x = "Within-bin pairwise Cα RMSD (Å)", y = NULL,
       title = "D  Not a variance artefact",
       subtitle = paste0("α loops vary as much as β loops,\n",
                         "so the single-chain comparison is not\n",
                         "driven by unequal spread.\n\n")) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 12 * FS, face = "bold", colour = "grey15"),
        plot.margin = margin(8, 46, 8, 8),
        plot.title = element_text(size = 15 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 8)))

combo <- (pA / pB / (pC | pD)) + plot_layout(heights = c(1.25, 0.95, 0.9))
save_plot("fig_pairing_summary.png", combo, 14, 17.5)
cat("\nFigure written to", OUT, "\n")
