# =============================================================================
# fig_pairing_evidence.R — paired structural input, five panels
#
#   Rscript scripts/fig_pairing_evidence.R [out_dir]
#
# Five panels from four tables in this repository. Reads no structures (those
# are on Zenodo, record 23765886).
#
# A  Purity against retention for paired, single-chain and no-framework
#    clustering, all binned identically so only the distance varies. Pairing
#    gives a modest improvement.
#
# B  Why A uses identical binning. Receptors are compared only within a bin, so
#    coarser binning puts more candidates in front of the distance function and
#    more chances to link different epitopes. Binning on three loops gives 26
#    bins for beta and 68 for alpha, against 477 paired. With each arm binned
#    on its own loops the paired-over-beta gap is +0.136; with identical
#    binning it is +0.045, and the apparent alpha-over-beta advantage
#    disappears.
#
# C  Under 62% and 85.7% repertoire noise the pairing advantage is of the same
#    order as on the clean cohort, and slightly negative at 30% retention.
#
# D  Polyspecificity (the fraction of clusters carrying more than one epitope)
#    is lower for paired input in the high-purity regime: 0.060 against
#    0.113 (CDR3α) and 0.118 (CDR3β) at 15% retention. The difference narrows
#    as retention rises and is gone by 30%, so the benefit of pairing is
#    concentrated at the operating point, and shows more clearly in
#    polyspecificity than in purity.
#
# E  The within-bin RMSD distributions of the alpha and beta loops are
#    near-identical, so the single-chain comparison is not affected by one
#    chain simply varying less than the other.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/pairing"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))
suppressMessages(library(patchwork))

T_ <- function(f) file.path(repo, "results/tables/extensions", f)
need <- c("chain_framework_ablation.csv", "chain_ablation_matched_binning.csv",
          "ladder_rung2_matched.csv", "ladder_rung3_matched.csv",
          "chain_rmsd_distribution.csv")
stopifnot(all(file.exists(T_(need))))
per  <- read_csv(T_(need[1]), show_col_types = FALSE)
mat  <- read_csv(T_(need[2]), show_col_types = FALSE)
r2   <- read_csv(T_(need[3]), show_col_types = FALSE)
r3   <- read_csv(T_(need[4]), show_col_types = FALSE)
dist <- read_csv(T_(need[5]), show_col_types = FALSE)

LV  <- c("Paired αβ", "Paired, no framework",
         "CDR3α only", "CDR3β only")
PAL <- setNames(c(SPTCE_RED, "#fb9a99", "#2166ac", "#7fa8cf"), LV)
tidy <- function(d) d %>%
  mutate(lab = factor(recode(arm, "paired_fw" = LV[1], "paired_noframe" = LV[2],
                             "alpha_fw" = LV[3], "beta_fw" = LV[4]),
                      levels = LV)) %>% filter(!is.na(lab))
mat <- tidy(mat); per <- tidy(per); r2 <- tidy(r2); r3 <- tidy(r3)

ip <- function(d, a, t, col = "purity") {
  x <- d %>% filter(arm == a) %>% arrange(retention)
  approx(x$retention, x[[col]], xout = t, ties = "ordered")$y
}

# ---- assertions ------------------------------------------------------------
expect_close("paired purity 0.75 A (published 0.9868)",
             per$purity[per$arm == "paired_fw" & per$threshold == 0.75], 0.9868, 0.001)
expect_close("matched beta purity 0.75 A", mat$purity[mat$arm == "beta_fw" & mat$threshold == 0.75], 0.9294, 0.002)
expect_close("paired polyspecificity at 15% retention", ip(mat, "paired_fw", 0.15, "polyspecific"), 0.060, 0.005)
expect_close("beta polyspecificity at 15% retention",   ip(mat, "beta_fw",   0.15, "polyspecific"), 0.118, 0.005)
expect_close("rung 3 paired purity at 15% retention",   ip(r3,  "paired_fw", 0.15), 0.838, 0.005)

base_thm <- function(ttl, sub) list(
  labs(title = ttl, subtitle = sub),
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.title = element_text(size = 13 * FS, face = "bold"),
        plot.subtitle = element_text(size = 10.5 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 7))))

# ---- A: purity-retention ---------------------------------------------------
pA <- ggplot(mat, aes(retention, purity, colour = lab)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.4) +
  scale_colour_manual(values = PAL, name = NULL,
                      guide = guide_legend(nrow = 2, byrow = TRUE,
                        override.aes = list(linewidth = 1.6, size = 0))) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0.05, 0.68)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.5, 1.03)) +
  labs(x = "Retention", y = "Cluster purity") +
  base_thm("A  Pairing helps, modestly",
           "All arms on the same 477 bins, so only the\ndistance varies. 4,031 receptors.") +
  theme(legend.position = "bottom", legend.text = element_text(size = 10.5 * FS))

# ---- B: the binning correction --------------------------------------------
gaps <- bind_rows(lapply(
  list(c("Paired − CDR3β", "paired_fw", "beta_fw"),
       c("Paired − CDR3α", "paired_fw", "alpha_fw"),
       c("CDR3α − CDR3β", "alpha_fw", "beta_fw")),
  function(v) tibble(comparison = v[1],
    Binning = c("Per-arm", "Matched"),
    gap = c(ip(per, v[2], 0.25) - ip(per, v[3], 0.25),
            ip(mat, v[2], 0.25) - ip(mat, v[3], 0.25))))) %>%
  mutate(Binning = factor(Binning, c("Per-arm", "Matched")),
         comparison = factor(comparison, rev(c("Paired − CDR3β",
                                               "Paired − CDR3α",
                                               "CDR3α − CDR3β"))))
pB <- ggplot(gaps, aes(gap, comparison, fill = Binning)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.6) +
  geom_col(position = position_dodge(width = 0.72), width = 0.6) +
  geom_text(aes(label = sprintf("%+.3f", gap)),
            position = position_dodge(width = 0.72), hjust = -0.12,
            size = 3.3 * FS, colour = "grey20") +
  scale_fill_manual(values = c("Per-arm" = "#bdbdbd", "Matched" = SPTCE_RED), name = NULL) +
  scale_x_continuous(limits = c(0, 0.175), breaks = seq(0, 0.15, 0.05)) +
  labs(x = "Purity gain at 25% retention", y = NULL) +
  base_thm("B  Two thirds of it was the binning",
           "Coarser binning puts more candidates in each\nbin. Controlling it shrinks every gap.") +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        axis.text.y = element_text(size = 11 * FS, face = "bold", colour = "grey15"))

# ---- C: noise ---------------------------------------------------------------
noise <- bind_rows(lapply(
  list(list("Clean", mat), list("62% noise", r2), list("85.7% noise", r3)),
  function(z) tibble(cohort = z[[1]], retention = c(0.15, 0.20, 0.25, 0.30),
    adv = vapply(c(0.15, 0.20, 0.25, 0.30), function(t)
      ip(z[[2]], "paired_fw", t) -
      max(ip(z[[2]], "alpha_fw", t), ip(z[[2]], "beta_fw", t), na.rm = TRUE),
      numeric(1))))) %>%
  mutate(cohort = factor(cohort, c("Clean", "62% noise", "85.7% noise")))
pC <- ggplot(noise, aes(factor(retention), adv, fill = cohort)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.6) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68) +
  scale_fill_manual(values = c("Clean" = SPTCE_RED, "62% noise" = "#f09d9c",
                               "85.7% noise" = "#9e9e9e"), name = NULL) +
  scale_y_continuous(limits = c(-0.02, 0.055)) +
  labs(x = "Retention", y = "Purity gain from pairing") +
  base_thm("C  Noise does not rescue it",
           "Pairing's advantage over the better single chain\nis the same order at 85.7% noise as when clean.") +
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom")

# ---- D: polyspecificity -----------------------------------------------------
pD <- ggplot(mat, aes(retention, polyspecific, colour = lab)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.4) +
  annotate("rect", xmin = 0.13, xmax = 0.21, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.18) +
  annotate("text", x = 0.17, y = 0.78, label = "operating\nregime",
           size = 3.3 * FS, colour = "grey35", lineheight = 0.95) +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0.05, 0.68)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.88)) +
  labs(x = "Retention", y = "Clusters carrying >1 epitope") +
  base_thm("D  Where pairing earns its keep",
           sprintf("Polyspecificity at 15%% retention: paired %.1f%%,\nα %.1f%%, β %.1f%%. The gap closes by 30%%.",
                   100 * ip(mat, "paired_fw", 0.15, "polyspecific"),
                   100 * ip(mat, "alpha_fw",  0.15, "polyspecific"),
                   100 * ip(mat, "beta_fw",   0.15, "polyspecific")))

# ---- E: the variance control ------------------------------------------------
D_LV <- c("CDR3α only", "CDR3β only", "Paired (87 residues)")
db <- dist %>% mutate(lab = factor(recode(selection, alpha = D_LV[1],
                                          beta = D_LV[2], paired = D_LV[3]),
                                   levels = D_LV))
medb <- db %>% group_by(lab) %>% summarise(m = median(rmsd), .groups = "drop")
pE <- ggplot(db, aes(rmsd, lab, fill = lab)) +
  geom_violin(colour = NA, alpha = 0.55, scale = "width", trim = TRUE) +
  geom_point(data = medb, aes(m, lab), size = 3, colour = "grey15", inherit.aes = FALSE) +
  geom_text(data = medb, aes(m, lab, label = sprintf("%.3f Å", m)),
            vjust = -1.4, size = 3.4 * FS, colour = "grey20", inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(c("#2166ac", "#7fa8cf", OTHER_GREY), D_LV),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 6), breaks = 0:6) +
  labs(x = "Within-bin pairwise Cα RMSD (Å)", y = NULL) +
  base_thm("E  Within-bin spread of α and β loops",
           paste0("α loops vary as much as β loops, so the single-chain ",
                  "comparison is not\ndriven by unequal spread.")) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 11 * FS, face = "bold", colour = "grey15"))

combo <- ((pA | pD) / (pB | pC) / pE) +
  plot_layout(heights = c(1.08, 0.92, 0.72)) +
  plot_annotation(
    title = "Does paired structural input earn its place?",
    subtitle = paste0("Cα RMSD over the CDR loops, single linkage, on the ",
                      "4,031 modelled receptors and two noise-ladder rungs. ",
                      "Every comparison is at\nmatched retention and matched ",
                      "binning, because the arms otherwise differ in two ways ",
                      "at once."),
    theme = theme(plot.title = element_text(size = 19 * FS, face = "bold"),
                  plot.subtitle = element_text(size = 12 * FS, colour = "grey35",
                                               margin = margin(t = 4, b = 12))))

save_plot("fig_pairing_evidence.png", combo, 14.5, 17)
cat("\nFigure written to", OUT, "\n")
