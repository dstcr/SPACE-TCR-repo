# =============================================================================
# fig01_rebuild.R — Figure 1
#
#   Rscript scripts/fig01_rebuild.R [out_dir]
#
# Both panels are drawn from tables in this repository, so no pipeline run is
# needed:
#   results/tables/published_method_metrics.csv
#   results/tables/sptce_all_threshold_metrics.csv
#
# 1A: point labels are placed with ggrepel so they do not overlap each other or
# the lines they annotate.
#
# 1B: a faceted bar chart. Position shows magnitude; colour shows only
# SPACE-TCR against the published methods; each value is printed outside the
# end of its bar in neutral ink, so it stays legible at any bar length.
#
# Colour: the house pair #e31a1c (RMSD) and #b30000 (RMSD+dih) differ by only
# dE 10.6 for a reader with normal colour vision. In 1B both SPACE-TCR rows are
# named on the y axis, so a single red is used; 1A keeps both, since the two
# modes there also differ in shape and line type.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/fig01_rebuild"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))

PUB_F <- file.path(repo, "results/tables/published_method_metrics.csv")
SW_F  <- file.path(repo, "results/tables/sptce_all_threshold_metrics.csv")
stopifnot(file.exists(PUB_F), file.exists(SW_F))
pub   <- read_csv(PUB_F, show_col_types = FALSE)
sweep <- read_csv(SW_F,  show_col_types = FALSE)

pretty_method <- function(x) str_replace_all(x, "_", " ")

save_plot <- function(name, plot, w, h, dpi = 300) {
  ggsave(file.path(OUT, name), plot, width = w, height = h, dpi = dpi)
  cat("  saved:", name, "\n")
}

# =============================================================================
# Panel A — purity against retention
# =============================================================================
sw <- sweep %>%
  mutate(mode_lab = ifelse(mode == "SPACE-TCR RMSD", "SPACE-TCR RMSD",
                           "SPACE-TCR RMSD+dih"),
         thr_lab  = sprintf("%.2f \u00c5", threshold))

clus_purity <- pub$purity[pub$method == "clusTCR"]

# Below 1.0 A the two modes sit almost on top of one another, so labelling both
# draws every threshold twice into the same few millimetres — which is what made
# the top-left corner unreadable. Label the RMSD arc throughout, and the gated
# arc only where it visibly separates from it.
sw_lab <- sw %>%
  filter(mode == "SPACE-TCR RMSD" | threshold >= 1.00)

labs_all <- bind_rows(
  pub %>% transmute(retention, purity, label = pretty_method(method),
                    lab_group = "published", face = "plain"),
  sw_lab %>% transmute(retention, purity, label = thr_lab,
                       lab_group = mode_lab, face = "bold")
)

pA <- ggplot() +
  geom_hline(yintercept = clus_purity, linetype = "dashed",
             colour = "grey55", linewidth = 0.45) +
  annotate("text", x = 0.995, y = clus_purity + 0.022,
           label = "ClusTCR purity", hjust = 1, size = 4.0 * FS,
           colour = "grey30", fontface = "italic") +
  geom_point(data = pub, aes(retention, purity),
             colour = OTHER_GREY, size = 3.2) +
  geom_line(data = sw, aes(retention, purity, colour = mode_lab,
                           linetype = mode_lab), linewidth = 1.1) +
  geom_point(data = sw, aes(retention, purity, colour = mode_lab,
                            shape = mode_lab), size = 3.6) +
  # ALL labels go through ONE repel layer. Two separate geom_text_repel calls do
  # not repel each other — they each solve their own layout — which is how
  # "GLIPH2" and the "1.25 A" annotation ended up printed against one another.
  geom_text_repel(data = labs_all,
                  aes(retention, purity, label = label,
                      colour = lab_group, fontface = face),
                  size = 3.8 * FS, show.legend = FALSE,
                  box.padding = 0.6, point.padding = 0.4,
                  min.segment.length = 0.25, segment.colour = "grey65",
                  segment.size = 0.3, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = c("SPACE-TCR RMSD" = SPTCE_RED,
                                 "SPACE-TCR RMSD+dih" = SPTCE_DARK,
                                 "published" = "grey25"),
                      breaks = c("SPACE-TCR RMSD", "SPACE-TCR RMSD+dih"),
                      name = NULL) +
  scale_discrete_identity(aesthetics = "fontface") +
  scale_linetype_manual(values = c("SPACE-TCR RMSD" = "solid",
                                   "SPACE-TCR RMSD+dih" = "dashed"), name = NULL) +
  scale_shape_manual(values = c("SPACE-TCR RMSD" = 16,
                                "SPACE-TCR RMSD+dih" = 15), name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.35, 1.09),
                     breaks = seq(0.4, 1.0, 0.1)) +
  labs(x = "Retention (% of input TCRs clustered)", y = "Global cluster purity") +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "bottom")

save_plot("fig01a_purity_retention.png", pA, 12, 8)

# =============================================================================
# Panel B — the four benchmark metrics
# =============================================================================
# The two SPACE-TCR operating points in the paper: RMSD at the focal 0.75 A,
# and RMSD+dih at 1.0 A (.sptce_bubble_dih_thr in the pipeline).
sptce_rows <- bind_rows(
  sweep %>% filter(mode == "SPACE-TCR RMSD",     threshold == 0.75) %>%
    mutate(method = "SPACE-TCR RMSD (0.75 \u00c5)"),
  sweep %>% filter(mode == "SPACE-TCR RMSD+dih", threshold == 1.00) %>%
    mutate(method = "SPACE-TCR RMSD+dih (1.0 \u00c5)")
) %>% select(method, retention, purity, perc_clusters_gt09, perc_tcrs_gt09)

stopifnot(nrow(sptce_rows) == 2)

METRIC_LABELS <- c(
  retention          = "Retention",
  purity             = "Purity",
  perc_clusters_gt09 = "% clusters\npurity \u2265 0.9",
  perc_tcrs_gt09     = "% TCRs in\nhigh-purity clusters"
)

allm <- bind_rows(
  sptce_rows %>% mutate(family = "SPACE-TCR"),
  pub %>% mutate(method = pretty_method(method), family = "Published method")
)

# One order for every facet, so a method sits on the same row throughout and the
# eye can track it across the four metrics. Purity is the ordering key because it
# is the metric the paper leads on.
ord <- allm %>% arrange(purity) %>% pull(method)

long <- allm %>%
  pivot_longer(c(retention, purity, perc_clusters_gt09, perc_tcrs_gt09),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, names(METRIC_LABELS), METRIC_LABELS),
         method = factor(method, levels = ord),
         family = factor(family, c("SPACE-TCR", "Published method")))

pB <- ggplot(long, aes(value, method, fill = family)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.18,
            size = 3.9 * FS, colour = "grey20") +
  # 2x2, not 1x4. On the page this figure is 6.2 in wide whatever the canvas,
  # so four facets side by side left each one ~1.5 in — not enough for the
  # strip label and the value label, both of which clipped. Two rows give each
  # facet twice the width and use the height a portrait figure already has.
  facet_wrap(~metric, nrow = 2) +
  scale_fill_manual(values = c("SPACE-TCR" = SPTCE_RED,
                               "Published method" = OTHER_GREY), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1"),
                     expand = expansion(mult = c(0, 0.30))) +
  labs(x = NULL, y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 13 * FS),
        axis.text.x = element_text(size = 12 * FS),
        strip.text = element_text(size = 13 * FS, face = "bold"),
        panel.spacing.x = unit(1.1, "lines"),
        legend.position = "bottom")

save_plot("fig01b_metric_bars.png", pB, 12, 9)

# ---- parity check against the tables the panels were built from -------------
chk_df <- long %>% filter(method == "SPACE-TCR RMSD (0.75 \u00c5)")
chk <- setNames(chk_df$value, as.character(chk_df$metric))
cat("\n  parity check, SPACE-TCR at the focal 0.75 A:\n")
for (n in names(chk)) cat(sprintf("    %-28s %.4f\n", gsub("\n", " ", n), chk[[n]]))
stopifnot(abs(chk[["Retention"]] - 0.14249843063402384) < 1e-9,
          abs(chk[["Purity"]]    - 0.986784140969163)   < 1e-9)
cat("  values match results/tables/sptce_all_threshold_metrics.csv exactly\n")

# =============================================================================
# Combined Figure 1 — A over B, panel tags baked in
# =============================================================================
# The published Figure 1 was an assembly: two separate images plus floating text
# boxes carrying the "A" and "B" labels. That is fragile — the tags are
# positioned absolutely, so they drift if either image is resized. The
# replacement is a single image with the tags drawn by patchwork.
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combo <- (pA / pB) +
    plot_layout(heights = c(0.82, 1.18)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 20 * FS, face = "bold"),
          plot.tag.position = c(0, 1))
  save_plot("fig01_combined.png", combo, 12, 15.0)
} else {
  cat("  [skip] fig01_combined: patchwork not installed\n")
}

cat("\nFigures written to", OUT, "\n")
