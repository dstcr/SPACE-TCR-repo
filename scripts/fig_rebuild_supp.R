# =============================================================================
# fig_rebuild_supp.R — supplementary Figures S1 to S9, rebuilt for legibility
#
#   Rscript scripts/fig_rebuild_supp.R [out_dir]
#
# Same principles and the same shared style as the main-text rebuilds
# (scripts/fig_common.R). Every panel is built
# from a committed table and asserts its headline values against the numbers
# the supplementary captions report.
#
# Three source tables were derived once from the stage-5 objects and committed
# so this reproduces from a clean clone:
#   extensions/dihedral_mode_metrics.csv        (S3)
#   extensions/nlvpmvatv_paired_members.csv     (S7)
#   extensions/cluster_overview_membership.csv  (S2)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/rebuild"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))
T <- function(...) file.path(repo, "results/tables", ...)

FOCAL <- "SPACE-TCR RMSD (0.75A)"
CFG6  <- c("GILGFVFTL","NLVPMVATV","GLCTLVAML","LLFGYPVYV","FLCMKALLL","NYNYLYRLF")
GATED <- "SPACE-TCR RMSD+dih (1A)"

# =============================================================================
# S1 — distinct epitopes captured at three quality thresholds
# =============================================================================
cat("\nS1 - epitope breadth\n")
ec <- read_csv(T("extensions/epitope_capture_stats.csv"), show_col_types = FALSE)
# Wrapped onto two lines: the one-line labels were clipped by the facet strip,
# losing the first word of the left panel and the last of the right.
LV <- c(n_distinct_epitopes_seen  = "Any non-singleton\ncluster",
        n_epitopes_hi_purity      = "Purity ≥ 0.9",
        n_epitopes_hi_purity_big  = "Purity ≥ 0.9\nand ≥ 5 members")
expect_close("S1 DeepTCR epitopes in any cluster (Section S2: 266)",
             ec$n_distinct_epitopes_seen[ec$method == "DeepTCR"], 266, 0.5)

s1 <- ec %>% select(method, names(LV)) %>%
  pivot_longer(-method, names_to = "level", values_to = "n") %>%
  mutate(level = factor(LV[level], levels = unname(LV)),
         family = fam(method), method = pretty_method(method))
ord1 <- s1 %>% filter(level == unname(LV)[1]) %>% arrange(n) %>% pull(method)
s1 <- s1 %>% mutate(method = factor(method, levels = ord1))

pS1 <- ggplot(s1, aes(n, method, fill = family)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5 * FS, colour = "grey20") +
  facet_wrap(~level, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = FAM_FILL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Distinct dominant epitopes captured", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing.x = unit(1.1, "lines"))
save_plot("figS1_epitope_breadth.png", pS1, 14, 7)

# =============================================================================
# S2 — cluster-level overview
# =============================================================================
# A circle packing: one dot per clustered TCR, one circle per cluster, coloured
# by epitope, three method rows against two columns (all epitopes / GILGFVFTL
# highlighted).
#
# There is no epitope legend: the panel shows that each cluster is internally
# one colour, which needs no key, and a legend of about a hundred epitopes
# would not be readable. Colours are recycled across epitopes (the subtitle
# says so), and GILGFVFTL is fixed to the house red so it can be found in both
# columns. Titles and type are sized for page width.
cat("\nS2 - cluster-level overview (original form, rebuilt)\n")
suppressMessages(library(ggforce))
ov <- read_csv(T("extensions/cluster_overview_membership.csv"), show_col_types = FALSE)

cl <- ov %>% filter(!is.na(epitope), epitope != "") %>%
  group_by(method, cluster_id) %>%
  summarise(size = n(), .groups = "drop") %>% filter(size >= 2)
chk <- cl %>% group_by(method) %>% summarise(n = n(), tcr = sum(size), .groups = "drop")
expect_close("S2 focal clusters (caption: 149)",  chk$n[chk$method == FOCAL], 149, 0.5)
expect_close("S2 focal receptors (caption: 681)", chk$tcr[chk$method == FOCAL], 681, 0.5)
expect_close("S2 clusTCR clusters (caption: 111)", chk$n[chk$method == "clusTCR"], 111, 0.5)

# A deterministic circle packer (packcircles is not installable for this R
# version, and a force layout differs between runs). The largest circle goes at
# the origin; each next one takes the free position nearest the origin, found
# by scanning rings outward.
pack <- function(r) {
  n <- length(r); ord <- order(-r)
  x <- numeric(n); y <- numeric(n); placed <- integer(0)
  step <- max(min(r) * 0.30, max(r) / 60)
  for (k in ord) {
    if (!length(placed)) { x[k] <- 0; y[k] <- 0; placed <- k; next }
    found <- FALSE; rad <- 0
    while (!found) {
      rad <- rad + step
      na <- max(12L, as.integer(ceiling(2 * pi * rad / step)))
      ang <- seq(0, 2 * pi, length.out = na + 1L)[-(na + 1L)]
      cx <- rad * cos(ang); cy <- rad * sin(ang)
      for (j in seq_along(cx)) {
        if (all(sqrt((x[placed] - cx[j])^2 + (y[placed] - cy[j])^2) -
                (r[placed] + r[k]) > 0.01)) {
          x[k] <- cx[j]; y[k] <- cy[j]; found <- TRUE; break
        }
      }
      if (rad > 200 * max(r)) { x[k] <- cx[1]; y[k] <- cy[1]; found <- TRUE }
    }
    placed <- c(placed, k)
  }
  data.frame(x = x, y = y, r = r)
}
# members inside a cluster circle, on concentric rings
ring_xy <- function(n) {
  if (n == 1) return(cbind(0, 0))
  out <- matrix(0, 0, 2); left <- n; ring <- 0
  while (left > 0) {
    k <- if (ring == 0) min(left, 1) else min(left, round(2 * pi * ring))
    th <- seq(0, 2 * pi, length.out = k + 1)[seq_len(k)]
    out <- rbind(out, cbind(ring * cos(th), ring * sin(th)))
    left <- left - k; ring <- ring + 1
  }
  out[seq_len(n), , drop = FALSE]
}

MLAB <- c("SPACE-TCR RMSD (0.75A)" = "SPACE-TCR RMSD 0.75 Å",
          "SPACE-TCR RMSD (1A)"    = "SPACE-TCR RMSD 1.0 Å",
          "clusTCR"                = "ClusTCR")
EPAL <- c("#e31a1c", "#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a", "#b15928",
          "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#8dd3c7")

build_row <- function(meth) {
  cc <- cl %>% filter(method == meth) %>% arrange(desc(size))
  lay <- pack(0.46 + 0.30 * sqrt(cc$size))
  cc <- bind_cols(cc, lay)
  mem <- ov %>% filter(method == meth, !is.na(epitope), epitope != "") %>%
    inner_join(cc %>% select(cluster_id, x, y, r), by = "cluster_id") %>%
    group_by(cluster_id) %>% mutate(i = row_number(),
                                    off = list(ring_xy(n()))) %>% ungroup() %>%
    rowwise() %>%
    mutate(px = x + off[i, 1] * (r / (max(2, sqrt(n())) * 1.9)),
           py = y + off[i, 2] * (r / (max(2, sqrt(n())) * 1.9))) %>%
    ungroup() %>% select(-off)
  eps <- sort(unique(mem$epitope))
  mem$col <- ifelse(mem$epitope == "GILGFVFTL", EPAL[1],
                    EPAL[2:12][(match(mem$epitope, setdiff(eps, "GILGFVFTL")) - 1) %% 11 + 1])
  cnt <- sprintf("%d clusters, %s TCRs", nrow(cc), format(sum(cc$size), big.mark = ","))
  ttl <- unname(MLAB[meth])
  base <- function(d, sub) ggplot() +
    geom_circle(data = cc, aes(x0 = x, y0 = y, r = r),
                colour = "grey60", fill = NA, size = 0.22) +
    geom_point(data = d, aes(px, py), colour = d$col, size = 1.25) +
    coord_equal() + labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 13 * FS, face = "bold"),
          plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35"))
  # Wrapped over several lines so it stays within the panel.
  left  <- base(mem, paste0(cnt, "\nDot = one TCR, circle = one cluster.",
                                 "\nColour = epitope (recycled)."))
  hl <- mem; hl$col <- ifelse(hl$epitope == "GILGFVFTL", EPAL[1], "grey78")
  right <- base(hl, paste0(cnt, "\nGILGFVFTL in red, all other epitopes grey.", "\n "))
  left | right
}
suppressMessages(library(patchwork))
pS2 <- build_row(FOCAL) / build_row("SPACE-TCR RMSD (1A)") / build_row("clusTCR")
save_plot("figS2_cluster_overview.png", pS2, 15, 17)

# =============================================================================
# S3 — dihedral-gate stringency across the RMSD sweep
# =============================================================================
cat("\nS3 - dihedral gate stringency\n")
dm   <- read_csv(T("extensions/dihedral_mode_metrics.csv"), show_col_types = FALSE)
rmsd <- read_csv(T("sptce_all_threshold_metrics.csv"), show_col_types = FALSE) %>%
  filter(mode == "SPACE-TCR RMSD") %>%
  transmute(mode = "RMSD only", threshold, retention, purity,
            perc_clusters_gt09, perc_tcrs_gt09)
dmx <- dm %>%
  mutate(mode = str_replace(mode, "SPACE-TCR RMSD\\+dih<", "+ gate < ") %>%
                str_replace("rad", " rad")) %>%
  select(mode, threshold, retention, purity, perc_clusters_gt09, perc_tcrs_gt09)
s3 <- bind_rows(rmsd, dmx) %>%
  mutate(mode = factor(mode, c("RMSD only", "+ gate < 0.5 rad",
                               "+ gate < 0.3 rad", "+ gate < 0.2 rad"))) %>%
  pivot_longer(c(retention, purity, perc_clusters_gt09, perc_tcrs_gt09),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
           c("retention", "purity", "perc_clusters_gt09", "perc_tcrs_gt09"),
           c("Retention", "Global purity", "% clusters\npurity \u2265 0.9",
             "% TCRs in\nhigh-purity clusters")))

expect_close("S3 RMSD-only purity at 0.25 A (Section S7: 100%)",
             rmsd$purity[rmsd$threshold == 0.25], 1.0, 0.002)
# The 0.5 rad gate is the setting Tables 1 and 2 report, so these values must
# match sptce_all_threshold_metrics.csv.
g05 <- dm %>% filter(grepl("0.5rad", mode))
expect_close("S3 0.5 rad gate, % clusters at purity >= 0.9 at 1.25 A (Section S7: 53.8%)",
             g05$perc_clusters_gt09[g05$threshold == 1.25], 0.53846, 0.002)
expect_close("S3 0.5 rad gate, global purity at 1.25 A (Section S7: 84.3%)",
             g05$purity[g05$threshold == 1.25], 0.84317, 0.002)

pS3 <- ggplot(s3, aes(threshold, value, colour = mode, linetype = mode)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.4) +
  facet_wrap(~metric, nrow = 1) +
  scale_colour_manual(values = c("RMSD only" = SPTCE_RED, "+ gate < 0.5 rad" = "#b30000",
                                 "+ gate < 0.3 rad" = "#777777", "+ gate < 0.2 rad" = "#bbbbbb"),
                      name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted"), name = NULL) +
  scale_x_continuous(breaks = c(0.25, 0.75, 1.25)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "RMSD threshold (Å)", y = NULL) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing.x = unit(1.6, "lines"), legend.position = "bottom")
save_plot("figS3_dihedral_sweep.png", pS3, 14, 7)

# =============================================================================
# S4 — what contaminates the GILGFVFTL-dominant clusters
# =============================================================================
cat("\nS4 - GILGFVFTL contamination\n")
ct <- read_csv(T("extensions/focal_epitope_contamination_GILGFVFTL.csv"), show_col_types = FALSE)
tot <- ct %>% group_by(method) %>% summarise(n = sum(n_tcrs), .groups = "drop")
expect_close("S4 focal SPACE-TCR non-GILGFVFTL receptors (caption: 2)",
             tot$n[tot$method == FOCAL], 2, 0.5)
expect_close("S4 clusTCR non-GILGFVFTL receptors (caption: 3)",
             tot$n[tot$method == "clusTCR"], 3, 0.5)

top <- ct %>% group_by(method) %>% slice_max(n_tcrs, n = 5, with_ties = FALSE) %>% ungroup()
# Short names: "SPACE-TCR RMSD (0.75A)" and "SPACE-TCR RMSD+dih (1A)" were
# clipped by the facet strip at this panel width, losing their leading letters.
SHORT4 <- c("SPACE-TCR RMSD (0.75A)"  = "SPACE-TCR 0.75 \u00c5",
            "SPACE-TCR RMSD+dih (1A)" = "SPACE-TCR+dih 1.0 \u00c5",
            "Hamming_Distance"        = "Hamming",
            "Levenshtein_Distance"    = "Levenshtein")
shorten <- function(x) ifelse(x %in% names(SHORT4), SHORT4[x], pretty_method(x))
top <- top %>%
  mutate(method = factor(shorten(method), levels = shorten(tot$method[order(tot$n)])),
         hl = epitope == "NLVPMVATV")

pS4 <- ggplot(top, aes(n_tcrs, reorder_within <- epitope, fill = hl)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n_tcrs), hjust = -0.25, size = 3.2 * FS, colour = "grey20") +
  facet_wrap(~method, scales = "free", ncol = 3) +
  scale_fill_manual(values = c(`TRUE` = SPTCE_RED, `FALSE` = OTHER_GREY),
                    labels = c(`TRUE` = "NLVPMVATV", `FALSE` = "other epitope"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Receptors placed in a GILGFVFTL-dominant cluster", y = NULL,
       caption = "Five most frequent contaminating labels per method; note the free x scale.") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 11 * FS),
        plot.caption = element_text(size = 12 * FS, hjust = 0, colour = "grey35"))
save_plot("figS4_gilgfvftl_contamination.png", pS4, 14, 10)

# =============================================================================
# S5 and S6 — PROMOTED TO THE MAIN TEXT
# =============================================================================
# The sequence-divergent share became Figure 4C and the GILGFVFTL coverage
# graph became Figure 3C (2026-10-11). Both are now built by
# scripts/fig_rebuild_all.R, beside the figures they are panels of.

# =============================================================================
# S5 — paired CDR3a/CDR3b logos for the NLVPMVATV clusters
# =============================================================================
cat("\nS5 - paired NLVPMVATV logos\n")
suppressMessages(library(ggseqlogo))
nl <- read_csv(T("extensions/nlvpmvatv_paired_members.csv"), show_col_types = FALSE)
cat(sprintf("    %d receptors in %d high-purity NLVPMVATV clusters\n",
            nrow(nl), n_distinct(nl$cluster_id)))

top_cl <- nl %>% count(cluster_id, sort = TRUE) %>% slice_head(n = 4) %>% pull(cluster_id)
mk_logo <- function(cid, chain) {
  s <- nl[[chain]][nl$cluster_id == cid]
  s <- s[!is.na(s) & nchar(s) == as.integer(names(sort(table(nchar(s)), decreasing = TRUE))[1])]
  ggseqlogo(s, method = "prob") +
    labs(title = sprintf("%s  (n = %d)",
                         ifelse(chain == "cdr3a", "CDR3α", "CDR3β"), length(s))) +
    # the chemistry legend is identical in all eight panels; it is stated in
    # the caption line instead of being drawn eight times
    theme(plot.title = element_text(size = 13 * FS, face = "bold"),
          axis.text = element_text(size = 11 * FS),
          legend.position = "none")
}
suppressMessages(library(patchwork))
rows <- lapply(top_cl, function(cid) mk_logo(cid, "cdr3a") | mk_logo(cid, "cdr3b"))
pS5_ <- Reduce(`/`, rows) +
  plot_annotation(caption = "Four largest high-purity NLVPMVATV clusters at the focal threshold; one cluster per row, paired α (left) and β (right). Residues are coloured by chemistry.",
                  theme = theme(plot.caption = element_text(size = 12 * FS,
                                                            hjust = 0, colour = "grey35")))
save_plot("figS5_paired_logos.png", pS5_, 13, 12)

# =============================================================================
# S6 — training-set overlap and held-out validation
# =============================================================================
cat("\nS6 - STCRDab overlap and leakage removal\n")
ovp <- read_csv(T("stage7/stcrdab_overlap_per_epitope.csv"), show_col_types = FALSE) %>%
  filter(set == "SPTCE_input", n_in_stcrdab > 0)
sm  <- read_csv(T("stage7/stcrdab_overlap_summary.csv"), show_col_types = FALSE)
expect_close("S6A cohort-wide overlap (caption: 45 of 4,292)", sm$n_overlap[1], 45, 0.5)
expect_close("S6A cohort-wide overlap percent (caption: 1.05%)", sm$pct_overlap[1], 1.048, 0.01)
expect_close("S6A GILGFVFTL overlap (caption: 8 of 454)",
             ovp$n_in_stcrdab[ovp$epitope == "GILGFVFTL"], 8, 0.5)

s8 <- ovp %>% arrange(pct_in_stcrdab) %>%
  mutate(epitope = factor(epitope, levels = epitope),
         six = epitope %in% CFG6)
pS6_a <- ggplot(s8, aes(100 * pct_in_stcrdab, epitope, fill = six)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%d / %d", n_in_stcrdab, n_total)),
            hjust = -0.15, size = 3.2 * FS, colour = "grey20") +
  scale_fill_manual(values = c(`TRUE` = SPTCE_RED, `FALSE` = OTHER_GREY),
                    labels = c(`TRUE` = "benchmark epitope", `FALSE` = "other"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "% of the epitope's unique CDR3\u03b2 also present in STCRDab", y = NULL,
       subtitle = sprintf("Cohort-wide: %d of %s unique modelled sequences (%.2f%%). Only the %d epitopes with any overlap are shown.",
                          sm$n_overlap[1], format(sm$n_total[1], big.mark = ","),
                          sm$pct_overlap[1], nrow(ovp))) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 11 * FS),
        axis.title.x = element_text(size = 13 * FS),
        plot.subtitle = element_text(size = 12 * FS, colour = "grey35"))

# ---- B: the focal metrics with the overlapping sequences removed -----------
lk <- read_tsv(T("stage7/leakage_per_tcr.tsv"), show_col_types = FALSE)
leaked <- unique(lk$cdr3b[lk$in_stcrdab])
foc <- read_csv(T("extensions/cluster_overview_membership.csv"), show_col_types = FALSE) %>%
  filter(method == FOCAL)
met <- function(d) {
  cl <- d %>% filter(!is.na(epitope), epitope != "") %>% group_by(cluster_id) %>%
    summarise(size = n(), pur = max(table(epitope)) / n(), .groups = "drop") %>%
    filter(size >= 2)
  # ASCII names here: R rejects a \u escape inside backticks. The reader-facing
  # label is applied as a factor below.
  tibble(purity = sum(cl$pur * cl$size) / sum(cl$size),
         hi = mean(cl$pur >= 0.9))
}
cmp <- bind_rows(met(foc) %>% mutate(set = "All 4,292 sequences"),
                 met(foc %>% filter(!cdr3b %in% leaked)) %>%
                   mutate(set = "STCRDab-overlapping sequences removed")) %>%
  pivot_longer(-set, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, c("purity", "hi"),
                         c("Global purity", "% clusters purity \u2265 0.9")),
         set = factor(set, c("All 4,292 sequences",
                             "STCRDab-overlapping sequences removed")))
d_pur <- diff(cmp$value[cmp$metric == "Global purity"])
expect_close("S6B change in focal purity, pp (caption: under 1 pp)",
             abs(100 * d_pur), 0.154, 0.02)

pS6_b <- ggplot(cmp, aes(value, metric, fill = set)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.68)) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * value)),
            position = position_dodge(width = 0.68), hjust = -0.15,
            size = 3.5 * FS, colour = "grey20") +
  scale_fill_manual(values = c(SPTCE_RED, OTHER_GREY), name = NULL,
                    guide = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("Removing all %d overlapping sequences moves global purity by %.2f pp and the high-purity share by %.2f pp.",
                          length(leaked), 100 * d_pur,
                          100 * diff(cmp$value[cmp$metric != "Global purity"]))) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        plot.subtitle = element_text(size = 12 * FS, colour = "grey35"))

save_plot("figS6_stcrdab_overlap.png",
          stack_panels(pS6_a, pS6_b, heights = c(2.1, 1), collect = FALSE), 13, 13)

# =============================================================================
# S7 — consensus clustering
# =============================================================================
cat("\nS7 - consensus clustering\n")
cm <- read_csv(T("extensions/consensus_metrics.csv"), show_col_types = FALSE)
expect_close("S7 focal purity on the shared universe (Table S3: 98.7%)",
             cm$purity[cm$method == FOCAL & !cm$is_consensus], 0.9872, 0.002)

cmx <- cm %>% mutate(kind = ifelse(is_consensus, "Consensus", "Parent method"),
                     lab = pretty_method(method))

# (A) the share of clusters at purity >= 0.9, the comparison the caption leads
# with.
a9 <- cmx %>% arrange(perc_clusters_gt09) %>% mutate(lab = factor(lab, levels = lab))
pS7_a <- ggplot(a9, aes(perc_clusters_gt09, lab, fill = kind)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * perc_clusters_gt09)),
            hjust = -0.15, size = 3.4 * FS, colour = "grey20") +
  scale_fill_manual(values = c(Consensus = SPTCE_RED,
                               `Parent method` = OTHER_GREY), name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(x = "% of clusters at purity \u2265 0.9", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        axis.title.x = element_text(size = 13 * FS))

# (B) what that extra purity costs in yield.
#
# No per-point labels here. Every method is already named on panel A's axis,
# and the four consensus settings sit inside a 99.4-99.7% purity by 2.2-4.4%
# retention box -- labelling them individually piled five long names on top of
# one another and told the reader nothing panel A had not. Both panels use
# fill, so one legend serves the pair.
# Both panels key on FILL alone. Panel A is bars, which have no shape, so a
# shape aesthetic here cannot merge with it and patchwork drew two legends,
# the second running off the page.
pS7_b <- ggplot(cmx, aes(retention, purity, fill = kind)) +
  geom_point(size = 4.6, shape = 21, colour = "grey25", stroke = 0.4) +
  scale_fill_manual(values = c(Consensus = SPTCE_RED,
                               `Parent method` = OTHER_GREY), name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Retention (share of the shared universe)",
       y = "Global cluster purity",
       subtitle = "Panel A names each point.") +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.subtitle = element_text(size = 12 * FS, colour = "grey35"))

save_plot("figS7_consensus.png",
          stack_panels(pS7_a, pS7_b, heights = c(1, 1), collect = TRUE), 12, 12)

cat("\nSupplementary figures written to", OUT, "\n")
