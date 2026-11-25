# =============================================================================
# fig_rebuild_all.R — main-text Figures 2 to 6
#
#   Rscript scripts/fig_rebuild_all.R [out_dir]
#
# Every panel is drawn from a table in this repository, so this needs no
# pipeline run, and each figure checks its headline values against the numbers
# in the paper before it is written.
#
# Style, shared with Figure 1:
#   - type at 1.3x, because the figure is scaled to ~0.52 on the page
#   - one repel layer per panel, so all labels repel one another
#   - no colour ramp on a quantity a bar length already shows
#   - SPACE-TCR in the house red, everything else muted
#   - values outside the mark in neutral ink, not inside a coloured fill
#   - panel tags drawn into the figure
#
# Figure 3B (mono-, bi-, poly- and highly polyspecific cluster fractions) is
# built from results/tables/extensions/polyspecific_classes.csv, whose
# monospecific column matches polyspecific_summary.csv for all ten methods.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/rebuild"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))

T <- function(...) file.path(repo, "results/tables", ...)

# =============================================================================
# FIGURE 2 — cluster-size behaviour
# =============================================================================
cat("\nFigure 2 - cluster-size behaviour\n")
ld <- read_csv(T("ld_all_clusters.csv"), show_col_types = FALSE)

FOCAL  <- "SPACE-TCR RMSD (0.75A)"
GATED  <- "SPACE-TCR RMSD+dih (1A)"
# One SPACE-TCR operating point per panel plus the published methods. The other
# thresholds live in the sweep figures; putting them here would invite a reader
# to compare a sweep against nine fixed methods on one axis.
keep <- c(FOCAL, unique(ld$method[ld$mode == "published"]))

# ---- 2A: composition of each method's clusters by size ----------------------
BINS <- c("2", "3", "4-5", "6-10", "11-20", ">20")
sizes <- ld %>% filter(method %in% keep) %>%
  mutate(bin = cut(n_members, c(1, 2, 3, 5, 10, 20, Inf), labels = BINS)) %>%
  count(method, bin) %>% group_by(method) %>%
  mutate(frac = n / sum(n), n_clusters = sum(n)) %>% ungroup() %>%
  mutate(bin = factor(bin, levels = BINS))

# Pair share decides the row order and is the number Section 3.3 quotes.
pair_share <- sizes %>% filter(bin == "2") %>% select(method, pair = frac)
expect_close("2A SPACE-TCR focal: share of clusters that are pairs",
             pair_share$pair[pair_share$method == FOCAL], 0.550, 0.002)
expect_close("2A GLIPH2: share of clusters that are pairs",
             pair_share$pair[pair_share$method == "GLIPH2"], 0.382, 0.002)

sizes <- sizes %>%
  mutate(method = factor(pretty_method(method),
                         levels = pretty_method(pair_share$method[order(-pair_share$pair)])))

lab2a <- sizes %>% filter(bin == "2")

p2a <- ggplot(sizes, aes(frac, method, fill = bin)) +
  # reverse = TRUE puts the SMALLEST bin at the origin. Without it ggplot
  # stacks ">20" first, so the pair share -- the number Section 3.3 quotes --
  # had to be measured between two interior boundaries instead of off the axis.
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  # Only the pair share is labelled. A number on every segment is unreadable at
  # print size and the legend plus axis carry the rest.
  geom_text(data = lab2a, aes(x = frac, label = sprintf("%.0f%%", 100 * frac)),
            hjust = -0.15, size = 3.6 * FS, colour = "grey20") +
  scale_fill_manual(values = setNames(ORDINAL_6, BINS), name = "Cluster size",
                    guide = guide_legend(nrow = 1)) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Fraction of clusters", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))

# ---- 2B: global purity as the minimum cluster size rises --------------------
MINS <- c(2, 3, 4, 5, 6, 8, 10, 15, 20)
purity_at <- function(df, k) {
  d <- df %>% filter(n_members >= k)
  if (!nrow(d)) return(NA_real_)
  sum(d$epitope_purity * d$n_members) / sum(d$n_members)
}
curve <- ld %>% filter(method %in% c(keep, GATED)) %>%
  group_by(method) %>%
  group_modify(~tibble(min_size = MINS,
                       purity = vapply(MINS, function(k) purity_at(.x, k), 0))) %>%
  ungroup() %>% filter(!is.na(purity))

expect_close("2B SPACE-TCR focal at minimum size 2 (Section 3.3: 98.7%)",
             curve$purity[curve$method == FOCAL & curve$min_size == 2], 0.987, 0.002)

curve <- curve %>% mutate(lab = pretty_method(method), sp = is_sptce(method))
ends <- curve %>% group_by(lab) %>% slice_max(min_size, n = 1) %>% ungroup()

p2b <- ggplot() +
  geom_line(data = filter(curve, !sp),
            aes(min_size, purity, group = lab),
            colour = OTHER_GREY, linewidth = 0.7, alpha = 0.85) +
  geom_line(data = filter(curve, sp),
            aes(min_size, purity, group = lab, colour = lab), linewidth = 1.4) +
  geom_point(data = filter(curve, sp),
             aes(min_size, purity, colour = lab), size = 2.6) +
  # Eleven categorical hues is an anti-pattern; the published methods are one
  # recessive grey band and are identified by a direct label at their endpoint.
  geom_text_repel(data = ends,
                  aes(min_size, purity, label = lab,
                      colour = ifelse(is_sptce(lab), lab, NA)),
                  size = 3.5 * FS, direction = "y", hjust = 0,
                  nudge_x = 0.9, segment.colour = "grey70", segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.42,
                  point.padding = 0.1, force = 3, max.iter = 20000,
                  max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c(SPTCE_RED, SPTCE_DARK),
                                        pretty_method(c(FOCAL, GATED))),
                      na.value = "grey25", breaks = pretty_method(c(FOCAL, GATED)),
                      name = NULL) +
  scale_x_continuous(breaks = MINS, limits = c(2, 27)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Minimum cluster size", y = "Global cluster purity") +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "bottom")

save_plot("fig02_cluster_size.png", stack_panels(p2a, p2b, heights = c(1, 1.1)),
          12, 13.5)

# =============================================================================
# FIGURE 3 — sensitivity and polyspecificity (panel B restored)
# =============================================================================
cat("\nFigure 3 - sensitivity and polyspecificity\n")
sens <- read_csv(T("corrected/sensitivity_CORRECTED_focal.csv"), show_col_types = FALSE)
rank <- read_csv(T("corrected/sensitivity_CORRECTED_ranking.csv"), show_col_types = FALSE)
ord3 <- rank$method[order(-rank$mean_corrected)]

# A and B are both 11-method row charts over the SAME methods, so they are set
# side by side on ONE shared ordering rather than stacked. Three panels down a
# page left the heatmap and the bars each using half the width and pushed the
# figure past a full page; this is roughly square and fits with its caption.
#
# The shared order is by mean corrected sensitivity, which is A's own sort. It
# also earns its keep: read across, the methods that recover the most receptors
# are the ones with the most polyspecific clusters, which is the trade-off
# Section 3.5 describes in words ("the methods that cluster most aggressively
# recover the most receptors").
ORD3 <- pretty_method(ord3)

# Section 3.5 quotes the GILGFVFTL column method by method. Pin the values the
# prose quotes so the text and this panel cannot disagree.
.g <- sens %>% filter(epitope == "GILGFVFTL")
.gv <- setNames(.g$sensitivity_corrected, .g$method)
for (nm in names(SEC35 <- c("DeepTCR" = 0.413, "SPACE-TCR RMSD+dih (1A)" = 0.407,
                            "GLIPH2" = 0.367, "Hamming_Distance" = 0.358,
                            "SPACE-TCR RMSD (0.75A)" = 0.339, "GIANA" = 0.310,
                            "iSMART" = 0.302, "Levenshtein_Distance" = 0.278,
                            "clusTCR" = 0.141, "TCRMatch" = 0.096,
                            "TCRdist3" = 0.044)))
  expect_close(paste("3A GILGFVFTL, Section 3.5 quotes:", nm), .gv[[nm]], SEC35[[nm]], 0.001)

s3 <- sens %>% mutate(method = factor(pretty_method(method), levels = rev(ORD3)),
                      lab = ifelse(is.na(sensitivity_corrected), "",
                                   sprintf("%.2f", sensitivity_corrected)))

p3a <- ggplot(s3, aes(epitope, method, fill = sensitivity_corrected)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = lab, colour = sensitivity_corrected > 0.5),
            size = 3.5 * FS, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_gradient(low = "#fff5eb", high = "#99000d", limits = c(0, 1),
                      breaks = c(0, 0.5, 1), na.value = "grey92",
                      name = "Sensitivity",
                      guide = guide_colourbar(barwidth = 9, barheight = 0.8,
                                              title.position = "top")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = "Benchmark epitope", y = NULL, title = "Epitope-level sensitivity") +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text.x = element_text(size = 12 * FS),
        axis.title.x = element_text(size = 13 * FS),
        plot.margin = margin(6, 8, 0, 30),
        plot.title = element_text(size = 14 * FS, face = "bold",
                                  margin = margin(t = 20, b = 6, l = 0)))

# ---- 3B: polyspecificity, on A's ordering --------------------------------
pc <- read_csv(T("extensions/polyspecific_classes.csv"), show_col_types = FALSE)
CLS <- c("Monospecific (1 epitope)", "Bispecific (2 epitopes)",
         "Polyspecific (3-5 epitopes)", "Highly polyspecific (>5 epitopes)")
stopifnot(all(pc$class %in% CLS))
mono <- pc %>% filter(class == CLS[1]) %>% select(method, mono = frac)
expect_close("3B SPACE-TCR focal: monospecific share (caption: >94%)",
             mono$mono[mono$method == FOCAL], 0.9463, 0.002)
expect_close("3B clusTCR: polyspecific share (Section 3.6: 6.3%)",
             1 - mono$mono[mono$method == "clusTCR"], 0.063, 0.002)
expect_close("3B TCRdist3: polyspecific share (Section 3.6: 92.7%)",
             1 - mono$mono[mono$method == "TCRdist3"], 0.927, 0.002)

pcx <- pc %>% mutate(class = factor(class, levels = CLS),
                     method = factor(pretty_method(method), levels = rev(ORD3)))
lab3b <- pcx %>% filter(class == CLS[1])

p3b <- ggplot(pcx, aes(frac, method, fill = class)) +
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  geom_text(data = lab3b, aes(x = frac, label = sprintf("%.0f%%", 100 * frac)),
            hjust = 1.12, size = 3.3 * FS, colour = "grey20") +
  scale_fill_manual(values = setNames(ORDINAL_4, CLS), name = NULL,
                    labels = c("Mono (1 epitope)", "Bi (2)", "Poly (3-5)",
                               "Highly poly (>5)"),
                    guide = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Fraction of clusters", y = NULL, title = "Cluster polyspecificity") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(size = 13 * FS),
        plot.margin = margin(6, 8, 0, 30),
        plot.title = element_text(size = 14 * FS, face = "bold",
                                  margin = margin(t = 20, b = 6, l = 0)))

# ---- 3C: GILGFVFTL coverage on a shared CDR3b sequence graph --------------
# Promoted from the supplement (2026-10-11). Section 3.5 discusses it, which is
# the section Figure 3 serves, so it belongs here rather than at the back.
suppressMessages({library(igraph); library(stringdist)})
om <- read_csv(T("overlay_membership_GILGFVFTL.csv"), show_col_types = FALSE)
stopifnot(nrow(om) == 454)
expect_close("3C receptors recovered by SPACE-TCR (Section 3.5: 141)",
             sum(om$capture %in% c("both", "SPACE-TCR only")), 141, 0.5)
expect_close("3C receptors recovered by ClusTCR (Section 3.5: 95)",
             sum(om$capture %in% c("both", "clusTCR only")), 95, 0.5)
expect_close("3C recovered by both (Section 3.5: 70)",
             sum(om$capture == "both"), 70, 0.5)

dmx2 <- stringdistmatrix(om$cdr3b, om$cdr3b, method = "lv")
g <- graph_from_adjacency_matrix(dmx2 <= 2 & dmx2 > 0, mode = "undirected",
                                 diag = FALSE)
set.seed(1)                       # the force layout is random; pin it
lay <- layout_with_fr(g)
nodes <- as_tibble(lay, .name_repair = ~c("x", "y")) %>%
  bind_cols(om %>% select(capture, sptce_cluster, clustcr_cluster, ld_degree)) %>%
  mutate(capture = factor(capture,
          c("both", "SPACE-TCR only", "clusTCR only", "neither")),
         # Connected component of the LD <= 2 graph = one "sequence island".
         # A cluster drawn across two components is a cluster that no
         # CDR3b-similarity criterion at this threshold could have formed.
         comp = components(g)$membership)
ed <- as_edgelist(g)
edges <- tibble(x = lay[ed[, 1], 1], y = lay[ed[, 1], 2],
                xend = lay[ed[, 2], 1], yend = lay[ed[, 2], 2])
CAP_COL <- c("both" = "#7f2704", "SPACE-TCR only" = SPTCE_RED,
             "clusTCR only" = "#4575b4", "neither" = "#cccccc")
cnt <- nodes %>% count(capture) %>% mutate(lab = sprintf("%s (%d)", capture, n))

# Three small multiples on ONE shared layout, not a single four-colour overlay.
#
# The overlay asked the reader to do the set algebra: to see what SPACE-TCR
# recovers they had to mentally union two of the four colours ("both" plus
# "SPACE-TCR only"), and the same again for ClusTCR, while the 288 grey
# "neither" nodes dominated the ink. Splitting it gives each method its own
# read at one hue, and keeps the four-way panel for the comparison itself.
#
# The node coordinates and the edge list are computed ONCE and reused, so the
# three panels are the same graph with different ink. That is what makes them
# comparable -- a per-panel force layout would move every node between panels
# and the figure would be worthless.
NOT_COL  <- "#e8e8e8"
CLUS_COL <- "#4575b4"
# The same grey must mean "not recovered" in all three panels, so the overlay
# panel's "neither" takes NOT_COL rather than the darker grey it used when it
# stood alone.
CAP_COL[["neither"]] <- NOT_COL

graph_frame <- function(title, subtitle, show_legend = FALSE) {
  list(
    geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
                 colour = "grey85", linewidth = 0.22),
    coord_equal(),
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle),
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = if (show_legend) "bottom" else "none",
          plot.margin = margin(0, 6, 4, 6),
          plot.title = element_text(size = 13 * FS, face = "bold",
                                    margin = margin(t = 4, b = 1)),
          plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35",
                                       margin = margin(t = 0, b = 4)))
  )
}

# ---- cluster membership, drawn as receptor-to-receptor links ---------------
# With 35 clusters there are too many to colour apart, so membership is shown
# by lines. Each cluster is drawn as a minimum spanning tree over its members in
# layout space: n-1 lines for n members, every endpoint a real receptor, and a
# cluster that crosses a gap draws the one long edge needed to bridge it. (A
# full clique would be 703 lines for the 38-member cluster alone; lines to a
# centroid would end at points where no receptor sits.)
#
# Every cluster is drawn, including those that stay within one island, so
# ClusTCR's panel shows its clusters as local rather than absent.
mst_links <- function(cluster_col) {
  nd <- nodes
  nd$cl <- nodes[[cluster_col]]
  nd <- filter(nd, !is.na(cl))
  stat <- nd %>% group_by(cl) %>%
    summarise(n = n(), n_island = n_distinct(comp), .groups = "drop")
  seg <- lapply(stat$cl[stat$n > 1], function(k) {
    m <- filter(nd, cl == k)
    gk <- graph_from_adjacency_matrix(as.matrix(dist(cbind(m$x, m$y))),
                                      mode = "undirected", weighted = TRUE,
                                      diag = FALSE)
    e <- as_edgelist(mst(gk), names = FALSE)
    tibble(x = m$x[e[, 1]], y = m$y[e[, 1]],
           xend = m$x[e[, 2]], yend = m$y[e[, 2]], cl = k)
  })
  seg <- if (length(seg)) bind_rows(seg) else
    tibble(x = numeric(), y = numeric(), xend = numeric(),
           yend = numeric(), cl = character())
  span_ids <- stat$cl[stat$n_island > 1]
  list(seg = mutate(seg, spans = cl %in% span_ids),
       nodes = mutate(nd, spans = cl %in% span_ids),
       stat = stat, span_ids = span_ids)
}

# One hue per island-spanning cluster, so scattered circles of the same cluster
# can be matched. Red (SPACE-TCR in this row) and ClusTCR's blue are excluded,
# so cluster colour is not confused with method colour.
SPAN_PAL <- c("#6a51a3", "#009E73", "#CC79A7", "#E69F00", "#252525", "#8c510a")

# One method, one hue. Non-recovered nodes are drawn FIRST so the signal sits
# on top of the context rather than under it.
one_method <- function(is_hit, col, title, cluster_col) {
  # Assigned outside the data mask: inside mutate() a bare `hit` would refer
  # to a column of that name if one existed.
  nd <- nodes
  nd$hit <- is_hit
  n_hit <- sum(nd$hit)
  sp <- mst_links(cluster_col)
  stopifnot(nrow(sp$nodes) == n_hit)
  n_cl <- nrow(sp$stat); n_span <- length(sp$span_ids)
  n_orph <- sum(sp$nodes$ld_degree == 0)
  cat(sprintf("    3C %-16s %d recovered, %d clusters, %d spanning >1 island, %d degree-0\n",
              cluster_col, n_hit, n_cl, n_span, n_orph))

  # NO CLUSTER-MEMBERSHIP LINKS. Three encodings were tried -- within-cluster
  # cliques, spokes to a cluster centroid, and a minimum spanning tree over
  # each cluster's members -- and all three were rejected, the last one only
  # after checking what the "islands" actually are.
  #
  # The LD <= 2 graph is not a set of comparable sequence families. It is ONE
  # component of 212, a second of 24, and 199 scraps of which 186 are single
  # sequences. Every SPACE-TCR cluster that "spans >= 2 islands" spans the
  # giant component plus one or two isolated sequences. So the lines were not
  # drawing bridges between sequence families -- there is only one family to
  # bridge from -- they were drawing the recruitment of orphan sequences, which
  # the node positions already show, and dressing a modest effect as a dramatic
  # one. The spanning counts are still computed and printed as an audit, and
  # live in results/tables/extensions/gilgfvftl_cluster_island_span.csv.
  ggplot() +
    graph_frame(title,
                sprintf("%d recovered in %d clusters\n%d sequence-orphans",
                        n_hit, n_cl, n_orph)) +
    geom_point(data = filter(nd, !hit), aes(x, y), shape = 21, fill = NOT_COL,
               colour = "grey72", stroke = 0.2, size = 2.1) +
    geom_point(data = filter(nd, hit), aes(x, y), shape = 21, fill = col,
               colour = "grey25", stroke = 0.3, size = 2.6)
}

p3c_sptce <- one_method(nodes$capture %in% c("both", "SPACE-TCR only"),
                        SPTCE_RED, "SPACE-TCR, 0.75 Å", "sptce_cluster")
p3c_clus  <- one_method(nodes$capture %in% c("both", "clusTCR only"),
                        CLUS_COL, "ClusTCR", "clustcr_cluster")

# "neither" goes through the fill scale here rather than being drawn as plain
# ink, so the key carries all four classes and the grey used in ALL THREE
# panels is named once for the row.
p3c_both <- ggplot() +
  graph_frame("Overlap", "70 shared, 96 unique",
              show_legend = TRUE) +
  geom_point(data = filter(nodes, capture == "neither"),
             aes(x, y, fill = capture), shape = 21, colour = "grey72",
             stroke = 0.2, size = 2.1) +
  geom_point(data = filter(nodes, capture != "neither"),
             aes(x, y, fill = capture), shape = 21,
             colour = "grey25", stroke = 0.3, size = 2.6) +
  scale_fill_manual(values = CAP_COL, name = NULL, drop = FALSE,
                    labels = setNames(cnt$lab, cnt$capture),
                    # NOT a right-hand legend: patchwork shares column widths
                    # down the figure, so putting C's key beside it stole
                    # width from A and B above and clipped their titles.
                    guide = guide_legend(nrow = 2,
                                         override.aes = list(size = 4.6)))

# Sub-panels carry bold titles rather than their own tags, so the row is
# "Figure 3C" as a whole, matching the caption.
#
# guides = "collect" moves the key out of the third panel and centres it under
# all three, so it is laid out against the full width.
p3c <- (p3c_sptce | p3c_clus | p3c_both) +
  patchwork::plot_layout(widths = c(1, 1, 1), guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(t = 2, b = 0),
        legend.text = element_text(size = 12 * FS))

save_plot("fig03_sensitivity_polyspecificity.png",
          (p3a | p3b) / p3c + patchwork::plot_layout(heights = c(0.95, 1.0)) +
            patchwork::plot_annotation(tag_levels = list(c("A", "B", "C"))) &
            # The default legend.box.spacing (2 * half_line) is ~0.75 in at
            # this base size, which leaves a gap between A/B's keys and row C;
            # it is reduced here.
            theme(plot.tag = element_text(size = 19 * FS, face = "bold"),
                  plot.tag.position = c(0, 1),
                  legend.box.spacing = unit(0.3, "lines")), 14, 10.0)
# Canvas width sets the on-page type size: the figure is placed at 6.2 in
# whatever it is drawn at, so 14 in lands type at 0.44 scale against ~0.52 for
# the other figures. 12.5 in was tried to bring it into line and is WORSE --
# panel A's method labels run into the heatmap and both legends clip. Three
# panels over two columns need the width more than they need the extra 15% of
# type size, so 14 in stands.

# =============================================================================
# FIGURE 4 — within-cluster CDR3 diversity
# =============================================================================
cat("\nFigure 4 - within-cluster CDR3 diversity\n")
beta <- ld %>% filter(method %in% c(keep, GATED), !is.na(median_ld))

expect_close("4A SPACE-TCR focal: median of per-cluster median CDR3b LD",
             median(beta$median_ld[beta$method == FOCAL]), 2, 0.001)

mk_violin <- function(df, xlab, order_by_median = TRUE) {
  med <- df %>% group_by(method) %>%
    summarise(m = median(median_ld), .groups = "drop") %>% arrange(m)
  d <- df %>% mutate(method = factor(pretty_method(method),
                                     levels = pretty_method(med$method)),
                     family = fam(method))
  ggplot(d, aes(median_ld, method, fill = family)) +
    geom_violin(scale = "width", linewidth = 0.3, colour = "grey35", alpha = 0.85,
                show.legend = FALSE) +
    geom_boxplot(width = 0.16, fill = "white", alpha = 0.9,
                 outlier.size = 0.5, linewidth = 0.35) +
    # No legend on the violin panels. Their alpha makes the key a different
    # spec from the bar key in panel C, so patchwork could not merge the two
    # guides and drew the legend twice, the second copy running off the page.
    # Panel C carries the only legend and it serves the whole figure.
    scale_fill_manual(values = FAM_FILL, name = NULL, guide = "none") +
    scale_x_continuous(breaks = seq(0, 16, 2)) +
    labs(x = xlab, y = NULL) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          legend.position = "bottom")
}
p4a <- mk_violin(beta, "Median within-cluster CDR3β Levenshtein distance")

AL <- T("extensions/within_cluster_ld_alpha.csv")
if (file.exists(AL)) {
  alpha <- read_csv(AL, show_col_types = FALSE)
  nm <- names(alpha)
  mcol <- intersect(c("median_ld", "median_ld_alpha", "median"), nm)[1]
  alpha <- alpha %>% rename(median_ld = !!mcol) %>%
    filter(method %in% c(keep, GATED), !is.na(median_ld))
  p4b <- mk_violin(alpha, "Median within-cluster CDR3\u03b1 Levenshtein distance")

  # ---- 4C: the sequence-divergent share -----------------------------------
  # Promoted from the supplement (2026-10-11). It measures the same thing the
  # violins do -- how far apart cluster members sit in CDR3b -- reduced to the
  # single number Section 3.4 quotes, so it belongs beside them.
  hi <- ld %>% filter(n_members >= 2, epitope_purity >= 0.9, !is.na(median_ld),
                      !str_starts(method, "SPACE-TCR") | method == FOCAL)
  dv <- hi %>% group_by(method) %>%
    summarise(n_hi = n(), pct = mean(median_ld > 3), .groups = "drop")
  expect_close("4C SPACE-TCR divergent share (Section 3.4: 22.7%)",
               dv$pct[dv$method == FOCAL], 0.227, 0.002)
  dv <- dv %>% mutate(family = fam(method), method = pretty_method(method)) %>%
    arrange(pct) %>% mutate(method = factor(method, levels = method))

  p4c <- ggplot(dv, aes(100 * pct, method, fill = family)) +
    geom_col(width = 0.7, key_glyph = "rect") +
    geom_text(aes(label = sprintf("%.1f%%  (n=%d)", 100 * pct, n_hi)),
              hjust = -0.1, size = 3.4 * FS, colour = "grey20") +
    scale_fill_manual(values = FAM_FILL, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.34))) +
    labs(x = "% of high-purity clusters with median within-cluster CDR3\u03b2 LD > 3",
         y = NULL,
         subtitle = "n = high-purity clusters the share is taken over") +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          axis.title.x = element_text(size = 12.5 * FS),
          plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35"))

  save_plot("fig04_within_cluster_ld.png",
            stack_panels(p4a, p4b, p4c, heights = c(1.5, 1, 1.35),
                         collect = TRUE), 12, 18)
} else {
  cat("  [warn] alpha table absent; Figure 4 is panel A only\n")
  save_plot("fig04_within_cluster_ld.png", p4a, 12, 9)
}

# =============================================================================
# FIGURE 5 — epitope purity and TRBV/TRBJ homogeneity
# =============================================================================
cat("\nFigure 5 - V/J gene homogeneity\n")
gc <- read_csv(T("extensions/gene_consensus_by_method.csv"), show_col_types = FALSE)
LBL <- c("Epitope (purity)", "V gene (consensus)", "J gene (consensus)")
gcx <- gc %>% filter(label %in% LBL) %>%
  mutate(label = factor(label, levels = LBL))

vg <- gcx %>% filter(label == "V gene (consensus)")
expect_close("5 SPACE-TCR focal: V-gene homogeneity (Section 3.7: 90%)",
             vg$perc_clusters_hi[vg$method == FOCAL], 0.90, 0.01)
jg <- gcx %>% filter(label == "J gene (consensus)")
expect_close("5 SPACE-TCR focal: J-gene homogeneity (Section 3.7: 66%)",
             jg$perc_clusters_hi[jg$method == FOCAL], 0.66, 0.01)

ordv <- vg %>% arrange(perc_clusters_hi) %>% pull(method)
gcx <- gcx %>% mutate(method = factor(pretty_method(method),
                                      levels = pretty_method(ordv)),
                      family = fam(method))

p5 <- ggplot(gcx, aes(perc_clusters_hi, method, fill = family)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * perc_clusters_hi)),
            hjust = -0.15, size = 3.5 * FS, colour = "grey20") +
  facet_wrap(~label, nrow = 1) +
  scale_fill_manual(values = FAM_FILL, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.32))) +
  labs(x = "% of non-singleton clusters in which one label covers ≥ 90% of members",
       y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.title.x = element_text(size = 13 * FS, face = "plain"),
        panel.spacing.x = unit(1.1, "lines"),
        legend.position = "bottom")

save_plot("fig05_gene_consensus.png", p5, 13, 8)

# =============================================================================
# FIGURE 6 — robustness to repertoire noise
# =============================================================================
cat("\nFigure 6 - repertoire-noise robustness\n")
lad <- read_csv(T("benchmark_a_ladder.csv"), show_col_types = FALSE,
                na = c("", "NA")) %>%
  mutate(pct_noise = as.numeric(pct_noise))

# RMSD-only SPACE-TCR at the focal threshold against the nine sequence
# methods; the paper's methods comparison is RMSD-only, so no gate is drawn.
focal_arm <- lad %>%
  filter(method == "SPACE-TCR", param == "0.75", is.na(dih)) %>%
  mutate(lab = "SPACE-TCR (0.75 Å)")
# The ladder table abbreviates the two edit-distance baselines; every other
# figure spells them out, so the names are made to match.
LADDER_NAMES <- c(HD = "Hamming Distance", LD = "Levenshtein Distance")
pub_arm <- lad %>% filter(!method %in% c("SPACE-TCR", "SPACE-TCR+dih")) %>%
  mutate(lab = ifelse(method %in% names(LADDER_NAMES),
                      LADDER_NAMES[method], pretty_method(method)))

expect_close("6 SPACE-TCR focal purity at the undiluted rung (Section 3.10: 96.7%)",
             focal_arm$purity[focal_arm$pct_noise == 0], 0.967, 0.003)
expect_close("6 SPACE-TCR focal purity at the densest rung (Section 3.10: 30.4%)",
             focal_arm$purity[which.max(focal_arm$pct_noise)], 0.304, 0.003)

# The ladder's rungs are not evenly spaced in noise (0, 62, 77, 83, 87, 89%),
# so on a linear axis the last four tick labels would overlap. They are six
# discrete subsets, so they are spaced evenly and each tick shows the
# unannotated share with the subset size beneath it, as in the figure caption.
rungs <- lad %>% distinct(subset, pct_noise, n_universe) %>% arrange(pct_noise) %>%
  mutate(step = row_number(),
         tick = sprintf("%.0f%%\n%s", pct_noise, format(n_universe, big.mark = ",")))
step_of <- setNames(rungs$step, as.character(rungs$pct_noise))

focal_arm <- focal_arm %>% mutate(step = step_of[as.character(pct_noise)])
pub_arm   <- pub_arm   %>% mutate(step = step_of[as.character(pct_noise)])

mk_ladder <- function(yvar, ylab, pct = TRUE, show_x_title = TRUE) {
  # ONE repel layer for every label, SPACE-TCR included. An annotate() label
  # cannot repel against a geom_text_repel layer, so labels would overlap.
  ends <- bind_rows(
    pub_arm   %>% group_by(lab) %>% slice_max(step, n = 1) %>% ungroup() %>%
      mutate(grp = "published"),
    focal_arm %>% slice_max(step, n = 1) %>% mutate(grp = "sptce")
  )
  p <- ggplot() +
    geom_line(data = pub_arm, aes(step, .data[[yvar]], group = lab),
              colour = OTHER_GREY, linewidth = 0.7, alpha = 0.8) +
    geom_line(data = focal_arm, aes(step, .data[[yvar]]),
              colour = SPTCE_RED, linewidth = 1.6) +
    geom_point(data = focal_arm, aes(step, .data[[yvar]]),
               colour = SPTCE_RED, size = 2.8) +
    geom_text_repel(data = ends,
                    aes(step, .data[[yvar]], label = lab,
                        colour = grp, fontface = ifelse(grp == "sptce", "bold", "plain")),
                    size = 3.4 * FS, direction = "y", hjust = 0, nudge_x = 0.28,
                    segment.colour = "grey70", segment.size = 0.3,
                    min.segment.length = 0, box.padding = 0.3, force = 2,
                    max.iter = 20000, max.overlaps = Inf, seed = 1,
                    show.legend = FALSE) +
    scale_colour_manual(values = c(published = "grey25", sptce = SPTCE_RED),
                        guide = "none") +
    scale_discrete_identity(aesthetics = "fontface") +
    scale_x_continuous(breaks = rungs$step, labels = rungs$tick,
                       limits = c(0.9, 9.2)) +
    labs(x = if (show_x_title)
           "Unannotated share of the repertoire (second line: sequences in the subset)"
         else NULL,
         y = ylab) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          axis.text.x = element_text(size = 12 * FS),
          axis.title.x = element_text(size = 13 * FS))
  if (pct) p <- p + scale_y_continuous(labels = percent_format(accuracy = 1),
                                       limits = c(0, 1))
  p
}
p6a <- mk_ladder("purity", "Cluster purity", show_x_title = FALSE)
p6b <- mk_ladder("mean_size", "Mean cluster size", pct = FALSE)

save_plot("fig06_noise_ladder.png", stack_panels(p6a, p6b, heights = c(1, 1)),
          12, 13.5)

cat("\nAll figures written to", OUT, "\n")
