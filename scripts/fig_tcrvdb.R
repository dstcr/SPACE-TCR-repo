# =============================================================================
# fig_tcrvdb.R — does structural clustering pick out annotations that validate?
#
#   Rscript scripts/fig_tcrvdb.R [out_dir]
#
# Draws the TCRvdb re-validation computed by scripts/tcrvdb_revalidation.py
# against the functionally screened labels of Messemaker et al. (Immunity
# 59:2775-2788, 2026). Reads only the committed summary tables, so it runs
# without the TCRvdb data, which is request-only and may not be redistributed.
#
# (A) Purity survives. The worry was that labels inherited from VDJdb are ~50%
#     wrong, and that co-reported receptors sharing a WRONG label would cluster
#     together and inflate purity. Purity does not fall on validated labels; it
#     rises. Three bases are shown, the third being the harshest reading, in
#     which non-validating receptors become a class of their own.
#
# (B) The result, drawn as the control demands. The pooled odds ratio of 7.91
#     is confounded: sequence publicness predicts validation by itself, from
#     0.280 at degree 0 to 0.739 at degree 10+ among unclustered receptors, and
#     it also predicts whether a receptor clusters at all. So the figure does
#     NOT show the pooled number. It shows the comparison WITHIN each stratum of
#     LD <= 2 sequence degree, for the structural method and the sequence
#     comparator side by side.
#
#     Each row is one stratum; the two dots are the validation rate among
#     receptors the method clustered and among those it did not; the bar between
#     them is the effect. SPACE-TCR's gaps stay open in three of four strata.
#     ClusTCR's close. That contrast is the claim, and a pooled bar chart would
#     have hidden it -- ClusTCR's pooled OR of 3.71 looks like a real effect
#     until it is stratified, at which point no stratum reaches significance and
#     one has an odds ratio of 0.94.
#
#     Non-significant strata are drawn hollow and greyed. ClusTCR has too few
#     clustered receptors at degree 0 to test at all, and that cell is labelled
#     rather than left blank.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "results/figures/tcrvdb"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

repo <- if (file.exists("scripts/fig_common.R")) "." else ".."
source(file.path(repo, "scripts/fig_common.R"))
suppressMessages(library(patchwork))

T_PUR <- file.path(repo, "results/tables/extensions/tcrvdb_purity_recomputed.csv")
T_STR <- file.path(repo, "results/tables/extensions/tcrvdb_enrichment_stratified.csv")
stopifnot(file.exists(T_PUR), file.exists(T_STR))
pur <- read_csv(T_PUR, show_col_types = FALSE)
str_ <- read_csv(T_STR, show_col_types = FALSE)

M_LAB <- c("SPACE-TCR RMSD (0.75A)" = "SPACE-TCR 0.75 Å",
           "SPACE-TCR RMSD (1A)"    = "SPACE-TCR 1.0 Å",
           "clusTCR"                = "ClusTCR")
B_LAB <- c("as_published"              = "As published",
           "validated_only"            = "Validated only",
           "non-validating_relabelled" = "Non-validating\nrelabelled")

# ---- assertions, so the figure cannot drift from the text ------------------
g <- function(m, b) pur$purity[pur$analysis == "focal_0.75A_cdr3b" &
                                 pur$method == m & pur$basis == b]
expect_close("focal purity, as published",   g("SPACE-TCR RMSD (0.75A)", "as_published"), 0.9963, 0.001)
expect_close("focal purity, validated only", g("SPACE-TCR RMSD (0.75A)", "validated_only"), 1.0000, 0.001)
expect_close("focal purity, relabelled",     g("SPACE-TCR RMSD (0.75A)", "non-validating_relabelled"), 0.9572, 0.001)

ld <- str_ %>% filter(control == "LD<=2 sequence degree")
expect_close("SPACE-TCR OR at degree 0",
             ld$odds_ratio[ld$method == "SPACE-TCR RMSD (0.75A)" & ld$stratum == "0"], 9.65, 0.05)
expect_close("clusTCR OR at degree 3-9",
             ld$odds_ratio[ld$method == "clusTCR" & ld$stratum == "3-9"], 0.94, 0.05)

# =============================== panel A ====================================
pa_d <- pur %>%
  filter(analysis == "focal_0.75A_cdr3b") %>%
  mutate(meth = factor(recode(method, !!!M_LAB), levels = unname(M_LAB)),
         basis = factor(recode(basis, !!!B_LAB), levels = unname(B_LAB)))

pA <- ggplot(pa_d, aes(basis, purity, fill = meth)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", purity)),
            position = position_dodge(width = 0.78), vjust = -0.6,
            size = 3.5 * FS, colour = "grey20") +
  scale_fill_manual(values = c(SPTCE_RED, "#fb9a99", OTHER_GREY), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.9, 1.022), breaks = seq(0.9, 1, 0.025),
                     oob = scales::squish) +
  labs(x = NULL, y = "Cluster purity",
       title = "Purity holds on validated labels",
       subtitle = "Five TCRvdb epitopes. Axis starts at 90%") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "bottom",
        axis.text.x = element_text(size = 12 * FS),
        plot.title = element_text(size = 14.5 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35"))

# =============================== panel B ====================================
LD_LV <- c("0", "1-2", "3-9", "10+")
pb <- ld %>%
  mutate(meth = factor(recode(method, !!!M_LAB), levels = unname(M_LAB[c(1, 3)])),
         stratum = factor(stratum, levels = rev(LD_LV)),
         sig = p_one_sided < 0.05,
         lab = ifelse(sig, sprintf("OR %.1f", odds_ratio),
                      sprintf("OR %.1f  n.s.", odds_ratio)))

# ClusTCR at degree 0 was not testable: too few clustered receptors. Say so in
# the panel rather than leaving a gap the reader has to interpret.
missing <- tibble(meth = factor("ClusTCR", levels = unname(M_LAB[c(1, 3)])),
                  stratum = factor("0", levels = rev(LD_LV)))

long <- pb %>%
  select(meth, stratum, clustered = rate_clustered,
         unclustered = rate_unclustered, sig) %>%
  pivot_longer(c(clustered, unclustered), names_to = "set", values_to = "rate")

pB <- ggplot(pb, aes(y = stratum)) +
  geom_segment(aes(x = rate_unclustered, xend = rate_clustered, yend = stratum,
                   colour = sig), linewidth = 2.4, alpha = 0.55) +
  geom_point(data = long, aes(x = rate, shape = set, colour = sig),
             size = 3.4, fill = "white", stroke = 1.1) +
  geom_text(aes(x = pmax(rate_clustered, rate_unclustered) + 0.035, label = lab,
                colour = sig), hjust = 0, size = 3.5 * FS, fontface = "bold") +
  geom_text(data = missing, aes(x = 0.03, label = "not testable: <5 clustered"),
            hjust = 0, size = 3.3 * FS, colour = "grey55", fontface = "italic") +
  facet_wrap(~meth, ncol = 1) +
  scale_colour_manual(values = c(`TRUE` = SPTCE_RED, `FALSE` = "#b0b0b0"),
                      guide = "none") +
  scale_shape_manual(values = c(clustered = 16, unclustered = 21), name = NULL,
                     labels = c("clustered by the method", "left unclustered")) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.26), breaks = seq(0, 1, 0.25)) +
  labs(x = "Receptors whose annotation validated functionally",
       y = "CDR3β sequence neighbours (Levenshtein ≤ 2)",
       title = "Structure predicts which annotations validate; sequence does not",
       subtitle = paste("Within each neighbourhood stratum, so the result is",
                        "not an artefact of clustering public sequences.\nGrey",
                        "= not significant")) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        legend.position = "bottom",
        strip.text = element_text(size = 13 * FS, face = "bold"),
        axis.title.y = element_text(size = 12.5 * FS),
        plot.title = element_text(size = 14.5 * FS, face = "bold"),
        plot.subtitle = element_text(size = 11.5 * FS, colour = "grey35"))

combo <- (pA / pB) + plot_layout(heights = c(0.82, 1.18)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 19 * FS, face = "bold"),
        plot.tag.position = c(0, 1))

save_plot("fig_tcrvdb.png", combo, 12, 12.5)
cat("\nFigure written to", OUT, "\n")
