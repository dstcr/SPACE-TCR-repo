# =============================================================================
# common_denominator_retention.R — retention for every method on one denominator
#
#   Rscript scripts/common_denominator_retention.R [out_csv]
#
# The nine published methods are not all measured in the same unit: some are
# keyed on the paired receptor, some on the unique CDR3b, and five were run on a
# pooled single-chain list, so their retention figures are not directly
# comparable (Figure 1 caption, Supplementary S13). This script puts every
# method on the same unit and denominator:
#
#   unit        unique CDR3b sequence
#   universe    the cohort's 4,292 unique CDR3b, the sequences SPACE-TCR was
#               given. Every method is restricted to these before cluster sizes
#               are taken, so no method is credited for sequences outside it.
#   numerator   unique CDR3b from that universe in a cluster of size >= 2,
#               size counted in the same unit
#   denominator 4,292
#
# This differs from the 6,163-sequence "shared universe" of Section 2.6 (the
# union of the published methods' outputs), which excludes 166 cohort
# sequences that no published method clustered.
#
# Purity is each method's own; it is computed within clusters and does not
# depend on the denominator.
#
# Inputs, all in this repository: ml/data/published_clusters.tsv.gz (the nine
# methods' cluster assignments, keyed on cdr3b), data/jouannet_labelled_input.tsv,
# and data/cohort_structural/ for SPACE-TCR's clustering.
# =============================================================================
source("sptce_metrics.R")

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args)) args[1] else "results/tables/retention_common_denominator.csv"

coh <- read_tsv("data/jouannet_labelled_input.tsv", show_col_types = FALSE)
pub <- read_tsv("ml/data/published_clusters.tsv.gz", show_col_types = FALSE)

universe <- unique(coh$cdr3_aa_beta)
N <- length(universe)
cat("universe: the cohort's own unique CDR3b =", N, "\n\n")

# ---- SPACE-TCR, from the committed geometry --------------------------------
fs <- list.files("data/cohort_structural/distance_tables", "\\.tsv\\.gz$", full.names = TRUE)
dih_map <- load_dihedral_map("data/cohort_structural/dihedrals.tsv.gz")
labels <- transmute(coh, receptor_id = as.character(id), cdr3b = cdr3_aa_beta,
                    epitope = epitope)

# Distinct CDR3b of the common universe in a cluster of >= 2 distinct CDR3b.
unique_cdr3b_clustered <- function(asg) {
  asg <- filter(asg, cdr3b %in% universe)
  if (!nrow(asg)) return(0L)
  jm_metrics(asg, N, unit = "cdr3b")$n_units_clustered
}
sptce_one <- function(thr, dih_cut = NA_real_)
  unique_cdr3b_clustered(cluster_at(fs, labels, thr, dih_cut,
                                    if (is.na(dih_cut)) NULL else dih_map,
                                    keep_unlabelled = TRUE))

sptce_rows <- bind_rows(
  map_dfr(c(0.5, 0.75, 1.0), ~ tibble(method = sprintf("SPACE-TCR RMSD (%g A)", .x),
                                      n_unique_cdr3b = sptce_one(.x))),
  tibble(method = "SPACE-TCR RMSD+dih (1 A)", n_unique_cdr3b = sptce_one(1.0, 0.5)))

# ---- the nine published methods, restricted to the same universe ------------
pub_rows <- pub %>%
  mutate(epitope = NA_character_) %>%
  split(.$method) %>%
  imap_dfr(~ tibble(method = .y, n_unique_cdr3b = unique_cdr3b_clustered(.x)))

tab <- bind_rows(sptce_rows, pub_rows) %>%
  mutate(retention_common = n_unique_cdr3b / N) %>%
  arrange(desc(retention_common))
print(as.data.frame(tab %>% mutate(retention_common = sprintf("%.2f%%", 100*retention_common))),
      row.names = FALSE)

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write_csv(mutate(tab, universe_n = N), OUT)
cat("\nwrote", OUT, "\n")
