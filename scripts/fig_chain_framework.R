# =============================================================================
# fig_chain_framework.R — what the second chain adds, with and without matched
#                         binning
#
#   Rscript scripts/fig_chain_framework.R [out_dir]
#
# Draws the Ca RMSD ablation computed by scripts/chain_framework_ablation.py on
# the 4,031 predicted structures. Reads only tables in this repository; the
# structures are on Zenodo (record 23765886).
#
# Receptors are compared only within a bin, so coarser binning puts more
# candidates in front of the distance function and more chances to link
# receptors of different epitopes. Binning on three loops instead of six gives
# 26 bins for beta and 68 for alpha, against 477 for the paired set.
#
#   Panel A  each arm binned on its own loops: compares whole methods, since a
#            single-chain method cannot know the other chain's loop lengths.
#   Panel B  every arm binned on the 477 six-loop bins, so only the distance
#            measure varies: compares representations.
#
# The two panels share axes so the difference between them is visible. Panel D
# shows that the within-bin RMSD distributions of the alpha and beta loops are
# near-identical.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/chain_framework"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))
suppressMessages(library(patchwork))

F_PER  <- file.path(repo, "results/tables/extensions/chain_framework_ablation.csv")
F_MAT  <- file.path(repo, "results/tables/extensions/chain_ablation_matched_binning.csv")
F_DIST <- file.path(repo, "results/tables/extensions/chain_rmsd_distribution.csv")
stopifnot(file.exists(F_PER), file.exists(F_MAT), file.exists(F_DIST))
per  <- read_csv(F_PER,  show_col_types = FALSE)
mat  <- read_csv(F_MAT,  show_col_types = FALSE)
dist <- read_csv(F_DIST, show_col_types = FALSE)

LV  <- c("Paired αβ", "CDR3α only", "CDR3β only",
         "Paired, no framework")
PAL <- setNames(c(SPTCE_RED, "#2166ac", "#7fa8cf", "#fb9a99"), LV)

tidy <- function(d) d %>%
  mutate(lab = factor(recode(arm,
      "paired_fw" = LV[1], "alpha_fw" = LV[2], "beta_fw" = LV[3],
      "paired_noframe" = LV[4]), levels = LV)) %>%
  filter(!is.na(lab))
per <- tidy(per); mat <- tidy(mat)

# The no-framework arm already bins on the six loops, because its SELECTION is
# the six loops -- the control changes nothing for it. So its curve is the same
# in both panels, and carrying it into B costs nothing and lets A and B share
# one legend instead of showing two slightly different ones.
mat <- bind_rows(mat, per %>% filter(arm == "paired_noframe") %>%
                        mutate(binned_on = "CDR_all"))

# ---- assertions ------------------------------------------------------------
gp <- function(d, a, t) d$purity[d$arm == a & d$threshold == t]
expect_close("paired purity at 0.75 A (published 0.9868)", gp(per, "paired_fw", 0.75), 0.9868, 0.001)
expect_close("per-arm   beta purity at 0.75 A", gp(per, "beta_fw", 0.75),  0.6735, 0.002)
expect_close("matched   beta purity at 0.75 A", gp(mat, "beta_fw", 0.75),  0.9294, 0.002)
expect_close("matched  alpha purity at 0.75 A", gp(mat, "alpha_fw", 0.75), 0.9034, 0.002)

interp <- function(d, a, t) {
  x <- d %>% filter(arm == a) %>% arrange(retention)
  approx(x$retention, x$purity, xout = t, ties = "ordered")$y
}
cat("\n  purity at matched retention:\n")
cat(sprintf("    %-9s %22s %22s\n", "retention", "per-arm binning", "matched binning"))
for (t in c(0.20, 0.25, 0.30, 0.38))
  cat(sprintf("    %9.2f   p %.3f a %.3f b %.3f   p %.3f a %.3f b %.3f\n", t,
      interp(per,"paired_fw",t), interp(per,"alpha_fw",t), interp(per,"beta_fw",t),
      interp(mat,"paired_fw",t), interp(mat,"alpha_fw",t), interp(mat,"beta_fw",t)))

# ---- panels A and B --------------------------------------------------------
curve <- function(d, ttl, sub) {
  ggplot(d, aes(retention, purity, colour = lab)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
    scale_colour_manual(values = PAL, name = NULL,
                        guide = guide_legend(nrow = 2, byrow = TRUE,
                          override.aes = list(linewidth = 1.6, size = 0))) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0.05, 0.88), breaks = seq(0.1, 0.8, 0.2)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0.25, 1.04), breaks = seq(0.3, 1, 0.1)) +
    labs(x = "Retention", y = "Cluster purity", title = ttl, subtitle = sub) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          legend.position = "bottom", legend.text = element_text(size = 11.5 * FS),
          plot.title = element_text(size = 13.5 * FS, face = "bold"),
          plot.subtitle = element_text(size = 11 * FS, colour = "grey35",
                                       margin = margin(t = 2, b = 8)))
}

pA <- curve(per, "A  Binned per arm",
            paste0("What a real single-chain method does.\n",
                   "26 bins for β, 68 for α, 477 paired, so\n",
                   "single chains are handicapped twice."))
pB <- curve(mat, "B  Binning held fixed",
            paste0("All arms on the same 477 bins, so only\n",
                   "the distance varies. The α-over-β gap goes\n",
                   "and the pairing gap shrinks to a third."))

# ---- panel C: how much of each gap was binning ----------------------------
gaps <- bind_rows(
  tibble(comparison = "Paired − CDR3β", binning = "Per-arm",
         gap = interp(per,"paired_fw",0.25) - interp(per,"beta_fw",0.25)),
  tibble(comparison = "Paired − CDR3β", binning = "Matched",
         gap = interp(mat,"paired_fw",0.25) - interp(mat,"beta_fw",0.25)),
  tibble(comparison = "Paired − CDR3α", binning = "Per-arm",
         gap = interp(per,"paired_fw",0.25) - interp(per,"alpha_fw",0.25)),
  tibble(comparison = "Paired − CDR3α", binning = "Matched",
         gap = interp(mat,"paired_fw",0.25) - interp(mat,"alpha_fw",0.25)),
  tibble(comparison = "CDR3α − CDR3β", binning = "Per-arm",
         gap = interp(per,"alpha_fw",0.25) - interp(per,"beta_fw",0.25)),
  tibble(comparison = "CDR3α − CDR3β", binning = "Matched",
         gap = interp(mat,"alpha_fw",0.25) - interp(mat,"beta_fw",0.25))) %>%
  mutate(binning = factor(binning, c("Per-arm", "Matched")),
         comparison = factor(comparison, rev(c("Paired − CDR3β",
                                               "Paired − CDR3α",
                                               "CDR3α − CDR3β"))))

pC <- ggplot(gaps, aes(gap, comparison, fill = binning)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.6) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%+.3f", gap)),
            position = position_dodge(width = 0.72),
            hjust = ifelse(gaps$gap >= 0, -0.15, 1.15),
            size = 3.8 * FS, colour = "grey20") +
  scale_fill_manual(values = c("Per-arm" = "#bdbdbd", "Matched" = SPTCE_RED),
                    name = NULL) +
  scale_x_continuous(limits = c(-0.02, 0.175), breaks = seq(0, 0.15, 0.05)) +
  labs(x = "Purity advantage at 25% retention", y = NULL,
       title = "C  Most of the chain effect was the binning, not the geometry",
       subtitle = paste0("67% of the paired-over-β gap and 59% of ",
                         "paired-over-α vanish once binning is held fixed.\n",
                         "The α-over-β advantage all but vanishes, and ",
                         "reverses sign below 22% retention.")) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 12.5 * FS, face = "bold", colour = "grey15"),
        legend.position = "bottom",
        plot.title = element_text(size = 13.5 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 8)))

# ---- panel D: within-bin RMSD spread of the alpha and beta loops ----------
D_LV <- c("CDR3α only", "CDR3β only", "Paired (87 residues)")
db <- dist %>% mutate(lab = factor(recode(selection, alpha = D_LV[1],
                                          beta = D_LV[2], paired = D_LV[3]),
                                   levels = D_LV))
medb <- db %>% group_by(lab) %>% summarise(m = median(rmsd), .groups = "drop")
expect_close("median within-bin RMSD, alpha", medb$m[medb$lab == D_LV[1]], 1.848, 0.01)
expect_close("median within-bin RMSD, beta",  medb$m[medb$lab == D_LV[2]], 1.861, 0.01)

pD <- ggplot(db, aes(rmsd, lab, fill = lab)) +
  geom_violin(colour = NA, alpha = 0.55, scale = "width", trim = TRUE) +
  geom_point(data = medb, aes(m, lab), size = 3.2, colour = "grey15",
             inherit.aes = FALSE) +
  geom_text(data = medb, aes(m, lab, label = sprintf("%.3f Å", m)),
            vjust = -1.5, size = 3.7 * FS, colour = "grey20", inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(c("#2166ac", "#7fa8cf", OTHER_GREY), D_LV),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 6), breaks = 0:6) +
  labs(x = "Within-bin pairwise Cα RMSD (Å)", y = NULL,
       title = "D  Within-bin spread of α and β loops",
       subtitle = paste0("α loops vary as much as β loops, so the single-chain ",
                         "difference\nis not explained by unequal spread.")) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 12 * FS, face = "bold", colour = "grey15"),
        plot.title = element_text(size = 13.5 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 8)))

# One legend for A and B, not two identical ones.
top <- (pA | pB) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combo <- top / pC / pD + plot_layout(heights = c(1.2, 0.85, 0.8))
save_plot("fig_chain_framework.png", combo, 13, 15.5)
cat("\nFigure written to", OUT, "\n")
