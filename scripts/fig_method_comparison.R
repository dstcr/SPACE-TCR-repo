# =============================================================================
# fig_method_comparison.R — purity vs retention for every method, on one axis
#
#   Rscript scripts/fig_method_comparison.R <comparison.csv> [out_dir]
#
# Input is the CSV written by ml/compare_all_methods.py, in which SPACE-TCR,
# the nine published methods, Levenshtein and the distilled model are all
# scored on the same held-out receptors against the same denominator.
#
# The figure makes the trade-off the table only implies: the methods that reach
# 100% purity do it by clustering almost nothing, and none of them forms a
# single sequence-divergent cluster. Point area is that divergent share, which
# is the paper's actual claim, so the claim is visible rather than tabulated.
# =============================================================================

suppressMessages({library(readr); library(dplyr); library(ggplot2); library(ggrepel)})

args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "ml_compare/comparison.csv"
OUT <- if (length(args) >= 2) args[2] else "results/figures/corrected"
if (!file.exists(IN))
  stop("not found: ", IN, "\n  run ml/compare_all_methods.py first", call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read_csv(IN, show_col_types = FALSE) %>%
  mutate(
    divergent = ifelse(is.na(pct_hi_divergent), 0, pct_hi_divergent),
    # The divergent share is a proportion of high-purity clusters, so for a
    # method with very few of them it rests on one or two clusters (the
    # distilled model shows 80% from a single cluster). Methods with fewer
    # than 5 high-purity clusters are marked as unreliable.
    n_hi = if ("n_clusters_hi" %in% names(.)) n_clusters_hi
           else round(n_clusters * pct_clusters_hi),
    unreliable = n_hi < 5,
    family = case_when(
      grepl("^SPACE-TCR", method)        ~ "SPACE-TCR (structural)",
      grepl("^distilled", method)        ~ "distilled model",
      method %in% c("DeepTCR", "TCRdist3") ~ "paired, learned or biophysical",
      TRUE                               ~ "CDR3 similarity"),
    family = factor(family, c("SPACE-TCR (structural)", "distilled model",
                              "paired, learned or biophysical", "CDR3 similarity")))

pal <- c("SPACE-TCR (structural)"         = "#C8553D",
         "distilled model"                = "#2E7D5B",
         "paired, learned or biophysical" = "#5B8FA8",
         "CDR3 similarity"                = "#9E7BB5")

# The denominator is the same for every row by construction, so recover it from
# any row rather than guessing from the largest clustered count.
denom <- round(stats::median(d$clustered / d$retention, na.rm = TRUE))
sub <- sprintf(paste0("every method scored on the same %s held-out receptors, ",
                      "same denominator, same equations.\n",
                      "point area = share of high-purity clusters too CDR3b-divergent ",
                      "for an edit distance to have formed (0%% for all CDR3-similarity methods)"),
               format(denom, big.mark = ","))

cap <- ""
if (any(d$unreliable)) {
  cap <- paste0(
    "Open circle: fewer than 5 high-purity clusters, so that method's divergent share rests on too small a base to read (",
    paste(sprintf("%s, n=%d", d$method[d$unreliable], d$n_hi[d$unreliable]),
          collapse = "; "), ").")
}

p <- ggplot(d, aes(100 * retention, 100 * purity)) +
  # Open circle where the divergent share rests on fewer than 5 high-purity
  # clusters, so an unreadable share cannot masquerade as a large filled point.
  geom_point(aes(size = 100 * divergent, fill = family, shape = unreliable),
             colour = "grey25", stroke = 0.6, alpha = 0.9) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 1), guide = "none") +
  ggrepel::geom_text_repel(aes(label = ifelse(unreliable,
                                              sprintf("%s (n=%d)", method, n_hi),
                                              method)), size = 2.9, min.segment.length = 0.1,
                           box.padding = 0.7, point.padding = 0.3, max.overlaps = 30,
                           force = 8, seed = 1,
                           segment.colour = "grey60", segment.size = 0.3) +
  scale_size_area(max_size = 15, name = "sequence-divergent\nclusters (%)",
                  breaks = c(0, 25, 50, 75)) +
  scale_fill_manual(values = pal, name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 4.5, shape = 21),
                             order = 1),
         size = guide_legend(order = 2)) +
  labs(title = "Purity is cheap if you cluster nothing",
       subtitle = sub,
       caption = cap,
       x = "retention (% of held-out receptors placed in a cluster of size 2 or more)",
       y = "purity (%)") +
  theme_minimal(base_size = 10) +
  theme(plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        plot.caption = element_text(size = 7, colour = "grey40", hjust = 0),
        panel.grid.minor = element_blank(),
        legend.position = "right")

f <- file.path(OUT, "fig_method_comparison_heldout.png")
ggsave(f, p, width = 8.6, height = 5.4, dpi = 300)
cat("wrote", f, "\n")

# ---- companion: divergent share as a plain ranked bar ------------------------
p2 <- ggplot(d %>% filter(!is.na(pct_hi_divergent)),
             aes(reorder(method, divergent), 100 * divergent, fill = family)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * divergent)), hjust = -0.15, size = 2.9) +
  scale_fill_manual(values = pal, name = NULL) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Sequence-divergent clusters, on the held-out receptors",
       subtitle = "share of each method's high-purity clusters whose median within-cluster CDR3b\nLevenshtein distance is 4 or more",
       x = NULL, y = "% of high-purity clusters") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        legend.position = "bottom")
f2 <- file.path(OUT, "fig_divergent_share_heldout.png")
ggsave(f2, p2, width = 7.4, height = 5.0, dpi = 300)
cat("wrote", f2, "\n")
