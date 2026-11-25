# =============================================================================
# fig_singlechain_vs_published.R — structural single-chain clustering beside
#                                  the nine published methods
#
#   Rscript scripts/fig_singlechain_vs_published.R [out_dir]
#
# Places the chain ablation (scripts/chain_framework_ablation.py) on the axes
# used for the method comparison: how does clustering one chain's loop geometry
# compare with the published sequence methods, and what does the second chain
# add?
#
# Units: unique CDR3-beta sequences clustered, over the 4,292-sequence universe
# (results/tables/retention_common_denominator.csv). Purity is each method's
# own, as in Table 2. The structural arms' numerators are derived from the
# cluster membership; the published points are drawn from
# retention_common_denominator.csv.
#
# Purity falls with retention for every method. The structural arms trace a
# curve because their threshold can be swept; the published methods are single
# points. Up and to the right is better.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/singlechain"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))

PUB  <- file.path(repo, "results/tables/published_method_metrics.csv")
COMM <- file.path(repo, "results/tables/retention_common_denominator.csv")
ABL  <- file.path(repo, "results/tables/extensions/chain_ablation_matched_binning.csv")
stopifnot(file.exists(PUB), file.exists(COMM), file.exists(ABL))
pub  <- read_csv(PUB,  show_col_types = FALSE)
comm <- read_csv(COMM, show_col_types = FALSE)
abl  <- read_csv(ABL,  show_col_types = FALSE)

# Published methods: purity from the metrics file, retention on the shared
# unique-CDR3b denominator rather than each method's own input size.
pubx <- pub %>%
  mutate(method = str_replace_all(method, "_", " ")) %>%
  left_join(comm %>% mutate(method = str_replace_all(method, "_", " ")) %>%
              select(method, retention_common), by = "method") %>%
  filter(!is.na(retention_common))
stopifnot(nrow(pubx) == 9)

ARM <- c(paired_fw = "SPACE-TCR, paired αβ",
         alpha_fw  = "SPACE-TCR, CDR3α loops only",
         beta_fw   = "SPACE-TCR, CDR3β loops only")
arms <- abl %>% filter(arm %in% names(ARM), !is.na(retention_common)) %>%
  mutate(lab = factor(recode(arm, !!!ARM), levels = unname(ARM)))
PAL <- setNames(c(SPTCE_RED, "#2166ac", "#7fa8cf"), unname(ARM))

# ---- assertions ------------------------------------------------------------
expect_close("ClusTCR purity (Table 2: 0.993)",
             pubx$purity[pubx$method == "clusTCR"], 0.993, 0.001)
expect_close("paired purity at 0.75 A (Table 1: 0.9868)",
             arms$purity[arms$arm == "paired_fw" & arms$threshold == 0.75], 0.9868, 0.001)
expect_close("paired retention_common at 0.75 A (derived 457/4292)",
             arms$retention_common[arms$arm == "paired_fw" & arms$threshold == 0.75],
             457 / 4292, 0.001)

cat("\n  structural arms at the focal 0.75 A, on the shared denominator:\n")
for (a in names(ARM)) {
  r <- arms %>% filter(arm == a, threshold == 0.75)
  cat(sprintf("    %-32s retention %.4f  purity %.4f  polyspec %.3f\n",
              ARM[[a]], r$retention_common, r$purity, r$polyspecific))
}

# ---- panel -----------------------------------------------------------------
focal <- arms %>% filter(threshold == 0.75)

p <- ggplot() +
  geom_line(data = arms, aes(retention_common, purity, colour = lab),
            linewidth = 1.2) +
  geom_point(data = arms, aes(retention_common, purity, colour = lab),
             size = 2.2) +
  # the published methods: one point each, named
  geom_point(data = pubx, aes(retention_common, purity),
             shape = 21, fill = "white", colour = "grey25",
             stroke = 1.0, size = 3.4) +
  geom_text_repel(data = pubx, aes(retention_common, purity, label = method),
                  size = 3.9 * FS, colour = "grey20", seed = 7,
                  box.padding = 0.55, point.padding = 0.35,
                  min.segment.length = 0.2, segment.colour = "grey70",
                  segment.size = 0.3, max.overlaps = Inf) +
  geom_point(data = focal, aes(retention_common, purity, colour = lab),
             size = 4.6, shape = 18) +
  geom_text_repel(data = focal, aes(retention_common, purity,
                                    label = sprintf("%.2f Å", threshold),
                                    colour = lab),
                  size = 3.6 * FS, fontface = "bold", show.legend = FALSE,
                  nudge_y = -0.045, seed = 9, min.segment.length = 0.2,
                  segment.colour = "grey70", segment.size = 0.3) +
  scale_colour_manual(values = PAL, name = NULL,
                      guide = guide_legend(nrow = 3, override.aes =
                                             list(linewidth = 1.6, size = 0))) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 0.95), breaks = seq(0, 0.9, 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.38, 1.02), breaks = seq(0.4, 1, 0.1)) +
  labs(x = "Unique CDR3β clustered, of the cohort's 4,292",
       y = "Cluster purity",
       title = "One chain's geometry against the nine published methods",
       subtitle = paste0("Structural arms sweep 0.4–1.25 Å; diamonds ",
                         "mark the focal 0.75 Å. Published methods are not ",
                         "tunable here,\nso each is one point. Up and to the ",
                         "right is better. Shared unique-CDR3β denominator ",
                         "throughout.")) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = "grey80"),
        legend.text = element_text(size = 12 * FS),
        plot.title = element_text(size = 15 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35",
                                     margin = margin(t = 2, b = 10)))

save_plot("fig_singlechain_vs_published.png", p, 13, 9)
cat("\nFigure written to", OUT, "\n")
