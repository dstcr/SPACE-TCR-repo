# =============================================================================
# fig_cluster_controls.R — is structural purity just V genes, or just length?
#
#   Rscript scripts/fig_cluster_controls.R [out_dir]
#
# Draws the permutation controls computed by scripts/cluster_controls.py. The
# clusters are held fixed and only the epitope labels move, within strata; each
# null asks what the clusters would look like if the geometry knew nothing but
# that stratifying variable. The distance from a null to the observed purity is
# the part of purity that variable does not explain.
#
# Each null has sd ~0.006 against gaps of 0.16-0.49, narrower than a line on a
# shared purity axis, so each null is drawn as a point with a segment to the
# observed value, and its sd printed beside it.
#
# Stratifying on V gene and CDR3 length together leaves 2,609 of 4,031
# receptors (65%) alone in a stratum, where they keep their true label, so that
# null sits close to the observed value by construction. Its small gap (+0.007)
# reflects a test with little power, not an absence of geometric signal; it is
# drawn in grey and labelled accordingly.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/cluster_controls"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))

F_IN <- file.path(repo, "results/tables/extensions/cluster_controls.csv")
stopifnot(file.exists(F_IN))
cc <- read_csv(F_IN, show_col_types = FALSE)

OBS <- unique(cc$purity_observed)
stopifnot(length(OBS) == 1)
N_TOT <- 4031

# A null is uninformative when most receptors cannot be permuted at all.
cc <- cc %>%
  mutate(stuck_frac = n_receptors_in_strata_of_one / N_TOT,
         degenerate = stuck_frac > 0.5,
         short = recode(null,
           "Chance (labels shuffled freely)" = "Chance",
           "Matched on V gene"               = "V gene",
           "Matched on CDR3 length"          = "CDR3 length",
           "Matched on V gene + CDR3 length" = "V gene + CDR3 length")) %>%
  arrange(null_mean) %>%
  mutate(short = factor(short, levels = short))

cat("\n  observed purity:", sprintf("%.4f", OBS), "\n")
print(as.data.frame(cc %>% select(short, null_mean, null_sd, gain_over_null,
                                  z_vs_null, stuck_frac, degenerate)))

expect_close("chance null",        cc$null_mean[cc$short == "Chance"],      0.2987, 0.01)
expect_close("CDR3 length null",   cc$null_mean[cc$short == "CDR3 length"], 0.4745, 0.01)
expect_close("V gene null",        cc$null_mean[cc$short == "V gene"],      0.6221, 0.01)
expect_close("observed purity",    OBS,                                     0.7839, 0.005)

lab <- cc %>% mutate(
  mid  = (null_mean + OBS) / 2,
  gtxt = ifelse(degenerate, "no power", sprintf("+%.3f", gain_over_null)),
  rtxt = sprintf("%.3f ± %.3f", null_mean, null_sd),
  note = ifelse(degenerate,
                sprintf("%.0f%% of receptors unpermutable", 100 * stuck_frac),
                sprintf("%.0f%% unpermutable", 100 * stuck_frac)))

p <- ggplot(cc, aes(y = short)) +
  geom_vline(xintercept = OBS, colour = SPTCE_RED, linewidth = 0.9) +
  # Anchored to the bottom row, not a numeric y: a numeric y here forces the
  # y scale continuous and the discrete `short` mapping then fails. The bottom
  # row is the Chance null at 0.30, so the space right of the line is empty.
  annotate("text", x = OBS + 0.014, y = levels(cc$short)[1], hjust = 0,
           vjust = 0.4, label = sprintf("observed\n%.3f", OBS),
           colour = SPTCE_RED, fontface = "bold", size = 4.1 * FS,
           lineheight = 0.95) +
  geom_segment(aes(x = null_mean, xend = OBS, yend = short,
                   colour = degenerate), linewidth = 1.5,
               arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  geom_point(aes(x = null_mean, colour = degenerate), size = 4.6) +
  geom_text(data = lab, aes(x = null_mean - 0.012, label = rtxt,
                            colour = degenerate),
            hjust = 1, size = 3.9 * FS, show.legend = FALSE) +
  geom_text(data = lab, aes(x = mid, label = gtxt, colour = degenerate),
            vjust = -1.0, fontface = "bold", size = 4.1 * FS,
            show.legend = FALSE) +
  geom_text(data = lab, aes(x = null_mean, label = note), colour = "grey45",
            hjust = 0.5, vjust = 2.4, size = 3.4 * FS, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "#2166ac", `TRUE` = "#b0b0b0"),
                      guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.08, 0.95), breaks = seq(0.1, 0.9, 0.1)) +
  labs(x = "Cluster purity", y = NULL,
       title = "Purity is not explained by V gene or CDR3 length",
       subtitle = paste0("Clusters held fixed; only epitope labels permuted, ",
                         "within strata.\n200 permutations each; the arrow is ",
                         "the purity that variable does NOT\naccount for. ",
                         "1,555 receptors in 321 clusters.")) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 14 * FS, face = "bold",
                                   colour = "grey15"),
        plot.title = element_text(size = 15 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 12)),
        plot.margin = margin(10, 16, 10, 8))

save_plot("fig_cluster_controls.png", p, 12.5, 7)
cat("\nFigure written to", OUT, "\n")
