# =============================================================================
# fig_chain_ablation.R — does native alpha-beta pairing carry joint information?
#
#   Rscript scripts/fig_chain_ablation.R [out_dir]
#
# SPACE-TCR keeps native alpha-beta pairing. This figure measures what pairing
# adds over clustering one chain.
#
# Each receptor's CDR3 loop is represented by its backbone dihedrals (IMGT
# 105-117, 13 positions x phi/psi = 26 angles per chain, from
# data/cohort_structural/dihedrals.tsv.gz). Distance is the pipeline's: mean
# circular distance over the angle vector, with the paired distance the mean of
# the two per-chain distances, as dihedral_distance() and composite_distance()
# in scripts/dihedral_util.py define them. Four representations are clustered
# with the same greedy procedure over the same threshold sweep:
#
#   alpha-only, beta-only, paired, and paired with alpha shuffled across
#   receptors (5 replicates).
#
# Averaging two distances reduces variance, so paired clustering could beat a
# single chain through denoising alone. The shuffled arm uses the same two
# distance matrices and the same averaging and changes only which alpha is
# paired with which beta, so the gap between paired and shuffled is the
# information carried by the pairing itself. The shuffled arm also falls below
# either single chain.
#
# Scope: this is the dihedral representation, not the Ca RMSD used for the
# paper's headline numbers. Absolute purity here is well below the 98.7% of the
# focal RMSD setting and is not comparable with Table 1; the figure shows the
# relative contribution of the two chains within one representation.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/chain_ablation"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))

F_IN <- file.path(repo, "results/tables/extensions/chain_ablation_dihedral.csv")
stopifnot(file.exists(F_IN))
ab <- read_csv(F_IN, show_col_types = FALSE)

LV <- c("Paired αβ", "CDR3β only", "CDR3α only",
        "α shuffled within epitope", "α shuffled, any receptor")
ab <- ab %>%
  mutate(rep_lab = factor(recode(representation,
      "paired"                                 = LV[1],
      "beta-only"                              = LV[2],
      "alpha-only"                             = LV[3],
      "paired (alpha shuffled within epitope)" = LV[4],
      "paired (alpha shuffled)"                = LV[5]), levels = LV))
stopifnot(!any(is.na(ab$rep_lab)))

# The two shuffled arms are 5 replicates each: mean line plus a range band.
SHUF_LV <- LV[4:5]
shuf <- ab %>% filter(rep_lab %in% SHUF_LV) %>%
  group_by(rep_lab, threshold) %>%
  summarise(retention = mean(retention), purity = mean(purity),
            lo = min(purity), hi = max(purity), .groups = "drop")
real <- ab %>% filter(!rep_lab %in% SHUF_LV)

PAL <- setNames(c(SPTCE_RED, "#2166ac", "#7fa8cf", "#6e6e6e", "#bdbdbd"), LV)

# ---- headline numbers, asserted so the figure cannot drift from the text ----
interp <- function(d, target) {
  d <- arrange(d, retention)
  approx(d$retention, d$purity, xout = target, ties = "ordered")$y
}
cat("\n  purity at matched retention:\n")
cat(sprintf("    %-9s %8s %8s %8s %12s %10s\n", "retention", "alpha", "beta",
            "paired", "shuf w/in ep", "shuf any"))
for (tg in c(0.05, 0.10, 0.20, 0.30, 0.50)) {
  a  <- interp(filter(real, rep_lab == LV[3]), tg)
  b  <- interp(filter(real, rep_lab == LV[2]), tg)
  p  <- interp(filter(real, rep_lab == LV[1]), tg)
  s4 <- interp(filter(shuf, rep_lab == LV[4]), tg)
  s5 <- interp(filter(shuf, rep_lab == LV[5]), tg)
  cat(sprintf("    %9.2f %8.3f %8.3f %8.3f %12.3f %10.3f\n", tg, a, b, p, s4, s5))
}
expect_close("paired purity at 10% retention",
             interp(filter(real, rep_lab == LV[1]), 0.10), 0.895, 0.01)
expect_close("CDR3beta-only purity at 10% retention",
             interp(filter(real, rep_lab == LV[2]), 0.10), 0.841, 0.01)
expect_close("alpha shuffled within epitope, purity at 10% retention",
             interp(filter(shuf, rep_lab == LV[4]), 0.10), 0.755, 0.02)
expect_close("alpha shuffled across any receptor, purity at 10% retention",
             interp(filter(shuf, rep_lab == LV[5]), 0.10), 0.581, 0.02)

# ---- panel ------------------------------------------------------------------
first_pt <- function(d) d %>% group_by(rep_lab) %>%
  filter(retention == min(retention[retention > 0.08])) %>%
  ungroup() %>% select(retention, purity, rep_lab)
lab_pts <- bind_rows(first_pt(real), first_pt(shuf))

p <- ggplot() +
  geom_ribbon(data = shuf, aes(retention, ymin = lo, ymax = hi,
                               group = rep_lab, fill = rep_lab), alpha = 0.25) +
  geom_line(data = shuf, aes(retention, purity, colour = rep_lab),
            linewidth = 1.0, linetype = "dashed") +
  scale_fill_manual(values = PAL, guide = "none") +
  geom_line(data = real, aes(retention, purity, colour = rep_lab),
            linewidth = 1.2) +
  geom_point(data = real, aes(retention, purity, colour = rep_lab), size = 2.1) +
  # A legend, not direct labels. Direct labels were tried: with five series the
  # curves converge above 40% retention and every label landed on a line.
  scale_colour_manual(values = PAL, name = NULL, breaks = LV,
                      guide = guide_legend(nrow = 2, byrow = TRUE,
                        override.aes = list(linewidth = 1.6, size = 0))) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.1), limits = c(0, 1.02)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0.1, 1, 0.1), limits = c(0.1, 1.02)) +
  labs(x = "Retention (% of receptors in non-singleton clusters)",
       y = "Cluster purity",
       title = "It matters which α goes with which β",
       subtitle = paste0("CDR3 backbone dihedrals, 4,031 labelled receptors.\n",
                         "Re-pairing β with another receptor's α costs purity ",
                         "even when that receptor\nsees the same epitope ",
                         "(bands = range over 5 shuffles)")) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "bottom",
        legend.text = element_text(size = 12.5 * FS),
        plot.title = element_text(size = 15 * FS, face = "bold"),
        plot.subtitle = element_text(size = 12 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 10)))

save_plot("fig_chain_ablation.png", p, 12.5, 8)
cat("\nFigure written to", OUT, "\n")
