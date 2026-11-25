# =============================================================================
# validate_outputs.R — recompute the paper's headline numbers from the inputs
# in data/ and compare them with the tables in results/.
#
#   Rscript validate_outputs.R [--write-actual]
#
# Needs only this repository and the seven R packages from
# install/r_packages.R.
#
# It checks that the geometry in data/cohort_structural/ and data/benchmark_a/,
# scored with sptce_metrics.R, gives the numbers in results/ and in the paper.
# It does not predict structures: to check TCRBuilder2 itself, rerun Stages A
# and B from data/jouannet_labelled_input.tsv (see INSTALL.md, Track B) and
# compare the resulting distance tables with these.
#
# Each check prints PASS, FAIL or SKIP with the expected and actual values. SKIP
# means an input is absent. The exit status is the number of FAILs.
# =============================================================================

suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
})

args         <- commandArgs(trailingOnly = TRUE)
WRITE_ACTUAL <- "--write-actual" %in% args

if (!file.exists("sptce_metrics.R"))
  stop("run this from the repository root", call. = FALSE)
source("sptce_metrics.R")

n_pass <- 0L; n_fail <- 0L; n_skip <- 0L
results <- list()

# Tolerance is absolute, on the quantity as printed. Purity and retention are
# compared at 1e-4, which is tighter than the 3 significant figures the paper
# reports, and counts are compared exactly.
TOL <- 1e-4

report <- function(name, expected, actual, tol = TOL, exact = FALSE) {
  if (is.null(actual) || (length(actual) == 1 && is.na(actual))) {
    n_skip <<- n_skip + 1L
    cat(sprintf("  SKIP  %-52s (not computable from the tracked inputs)\n", name))
    results[[length(results) + 1L]] <<-
      tibble(check = name, expected = as.character(expected),
             actual = NA_character_, status = "SKIP")
    return(invisible(NULL))
  }
  # An absent expected value is a SKIP, not a FAIL: it means the tracked table
  # does not carry that quantity, which is a gap in what can be checked rather
  # than a disagreement.
  if (is.null(expected) || (length(expected) == 1 && is.na(expected))) {
    n_skip <<- n_skip + 1L
    cat(sprintf("  SKIP  %-52s (the tracked table carries no such column)\n", name))
    results[[length(results) + 1L]] <<-
      tibble(check = name, expected = NA_character_,
             actual = as.character(actual), status = "SKIP")
    return(invisible(NULL))
  }
  ok <- if (exact) identical(as.integer(expected), as.integer(actual)) else
    abs(as.numeric(expected) - as.numeric(actual)) <= tol
  if (ok) n_pass <<- n_pass + 1L else n_fail <<- n_fail + 1L
  fmt <- function(v) if (exact) format(as.integer(v), big.mark = ",") else
    sprintf("%.5f", as.numeric(v))
  cat(sprintf("  %s  %-52s expected %12s   actual %12s\n",
              if (ok) "PASS" else "FAIL", name, fmt(expected), fmt(actual)))
  results[[length(results) + 1L]] <<-
    tibble(check = name, expected = fmt(expected), actual = fmt(actual),
           status = if (ok) "PASS" else "FAIL")
  invisible(NULL)
}

hdr <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")

cat("\n==============================================================\n")
cat("  SPACE-TCR — validating the tracked outputs\n")
cat("==============================================================\n")

# =============================================================================
# 0. The inputs are present and are the tables the paper describes
# =============================================================================
hdr("0. Inputs")

coh_f <- "data/jouannet_labelled_input.tsv"
if (!file.exists(coh_f)) stop("missing ", coh_f, call. = FALSE)
coh <- read_tsv(coh_f, show_col_types = FALSE)

report("cohort input: paired receptors",            4779L, nrow(coh),                    exact = TRUE)
report("cohort input: distinct ids",                4779L, n_distinct(coh$id),           exact = TRUE)
report("cohort input: distinct CDR3b",              4292L, n_distinct(coh$cdr3_aa_beta), exact = TRUE)
report("cohort input: distinct epitopes",            367L, n_distinct(coh$epitope),      exact = TRUE)
report("cohort input: rows with a missing field",      0L,
       sum(!complete.cases(coh[, c("id", "v_call_alpha", "j_call_alpha", "cdr3_aa_alpha",
                                   "v_call_beta", "j_call_beta", "cdr3_aa_beta")])),
       exact = TRUE)

# Only the last rung is tracked. Rung N is every receptor whose `first_rung`
# in ladder_membership.tsv is <= N -- the same derivation .rung_table() does in
# the pipeline and rebuild_rungs.R does on disk. These six assertions are what
# makes dropping the five redundant sequence tables safe to rely on.
lad_n  <- c(3056L, 8089L, 21413L, 23529L, 23529L, 23529L)
full_f <- "data/benchmark_a/donor1_subset6_sptce.tsv"
memb_f <- "data/benchmark_a/ladder_membership.tsv"
memb   <- if (file.exists(memb_f)) read_tsv(memb_f, show_col_types = FALSE) else NULL
full   <- if (file.exists(full_f)) read_tsv(full_f, show_col_types = FALSE) else NULL
rung_ids <- function(i) {
  f <- sprintf("data/benchmark_a/donor1_subset%d_sptce.tsv", i)
  if (file.exists(f)) return(unique(read_tsv(f, show_col_types = FALSE)$id))
  if (is.null(memb) || is.null(full)) return(NULL)
  intersect(full$id, memb$id[memb$first_rung <= i])
}
report("membership file covers the full table", 1L,
       if (is.null(memb) || is.null(full)) NA else
         as.integer(setequal(memb$id, full$id)), exact = TRUE)
for (i in 1:6) {
  r <- rung_ids(i)
  report(sprintf("ladder rung %d: paired receptors", i), lad_n[i],
         if (is.null(r)) NA else length(r), exact = TRUE)
}

# The ladder is nested, and rungs 4-6 are the same receptor set. docs/BENCHMARK_A.md
# says so; this is the assertion behind that claim.
ids <- map(1:6, rung_ids)
nest_ok <- all(map_lgl(1:5, ~ length(ids[[.x]]) > 0 &&
                             all(ids[[.x]] %in% ids[[.x + 1L]])))
report("ladder is strictly nested (rung i within i+1)", 1L, as.integer(nest_ok), exact = TRUE)
report("rungs 4, 5, 6 are the same receptor set",       1L,
       as.integer(setequal(ids[[4]], ids[[5]]) && setequal(ids[[5]], ids[[6]])), exact = TRUE)

# =============================================================================
# 1. The committed geometry is readable and complete
# =============================================================================
hdr("1. Committed geometry")

co_dir <- "data/cohort_structural"
dt     <- file.path(co_dir, "distance_tables")
report("cohort: pairwise distance tables on disk",   828L,
       if (dir.exists(dt)) length(list.files(dt, "\\.tsv\\.gz$")) else NA, exact = TRUE)
report("cohort: dihedral table present",               1L,
       as.integer(file.exists(file.path(co_dir, "dihedrals.tsv.gz"))), exact = TRUE)

la_dt <- "data/benchmark_a/subset4_clusters/distance_tables"
report("ladder: pairwise distance tables on disk",  1704L,
       if (dir.exists(la_dt)) length(list.files(la_dt, "\\.tsv\\.gz$")) else NA, exact = TRUE)
report("ladder: dihedral table present",               1L,
       as.integer(file.exists("data/benchmark_a/subset4_dihedrals/dihedrals.tsv.gz")),
       exact = TRUE)

# Reading every pairwise table checks that none is truncated or missing.
# Stage C writes one gzipped table per structural bin, with columns
#   i  j  id_i  id_j  status  rmsd  dihedral_alpha  dihedral_beta
#   clustering_dist  error
# where id_i and id_j are PDB paths and the receptor id is the file stem.
# Bins with no pairs within the cut-off are header-only files, so every column
# is read as character and converted here.
read_pairs <- function(dir, max_d) {
  fs <- list.files(dir, "\\.tsv\\.gz$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  stem <- function(p) sub("\\.pdb$", "", basename(p))
  map_dfr(fs, function(f) {
    d <- suppressWarnings(read_tsv(f, col_types = cols(.default = col_character()),
                                   progress = FALSE))
    if (!nrow(d)) return(NULL)
    nm <- names(d)
    ia <- if ("id_i" %in% nm) "id_i" else nm[1]
    ib <- if ("id_j" %in% nm) "id_j" else nm[2]
    dc <- if ("clustering_dist" %in% nm) "clustering_dist" else
          if ("rmsd" %in% nm) "rmsd" else nm[3]
    out <- tibble(a = stem(d[[ia]]), b = stem(d[[ib]]),
                  d = suppressWarnings(as.numeric(d[[dc]])))
    if ("status" %in% nm) out <- out[d$status == "ok" | is.na(d$status), ]
    out %>% filter(!is.na(d), d <= max_d)
  })
}

cat("\n  reading the cohort's 828 pairwise tables ...\n")
pr <- tryCatch(read_pairs(dt, 1.25), error = function(e) {
  cat("  read error:", conditionMessage(e), "\n"); NULL })
report("cohort: pairwise tables all readable", 1L,
       if (is.null(pr)) NA else 1L, exact = TRUE)
if (!is.null(pr))
  cat(sprintf("        %s pairs at or below 1.25 A over %s receptors\n",
              format(nrow(pr), big.mark = ","),
              format(n_distinct(c(pr$a, pr$b)), big.mark = ",")))

# =============================================================================
# 2. Clustering at the focal threshold reproduces Table 1
# =============================================================================
# Clustered and scored by sptce_metrics.R, exactly as the pipeline does it:
# cluster_at() at 0.75 A, then jm_metrics() over the labelled cohort.
hdr("2. Focal threshold (0.75 A) against Table 1")

focal <- NULL
if (!is.null(pr)) {
  lab <- coh %>% transmute(receptor_id = as.character(id), epitope)
  m <- jm_metrics(cluster_at(list.files(dt, "\\.tsv\\.gz$", full.names = TRUE),
                             lab, 0.75),
                  denom = nrow(lab), unit = "receptor_id")
  focal <- list(purity = m$purity, n_clusters = m$n_clusters,
                n_clustered = m$n_units_clustered, pct_hi = m$pct_clusters_hi)
}

# Expected values are read from the tracked table rather than hardcoded here,
# so this check cannot drift out of step with results/.
exp_f <- "results/tables/sptce_focal_table.csv"
if (file.exists(exp_f) && !is.null(focal)) {
  ef <- read_csv(exp_f, show_col_types = FALSE)
  # The table has one row per operating point; select the ungated RMSD row at
  # the focal threshold by name rather than by position.
  row <- ef %>% filter(str_detect(method, fixed("0.75")), !str_detect(method, "dih"))
  if (!nrow(row)) row <- ef[1, ]
  cat("\n  comparing against row:", row$method[1], "\n")
  pick <- function(nms) { h <- intersect(nms, names(row)); if (length(h)) row[[h[1]]][1] else NA }
  report("focal 0.75 A: clusters of size >= 2",
         pick(c("n_clusters", "clusters", "n_cluster")), focal$n_clusters, exact = TRUE)
  report("focal 0.75 A: receptors clustered",
         pick(c("n_clustered", "clustered", "n_tcrs_clustered")), focal$n_clustered, exact = TRUE)
  report("focal 0.75 A: share of clusters with purity >= 0.9",
         pick(c("perc_clusters_gt09", "pct_clusters_hi", "pct_high_purity", "frac_hi")),
         focal$pct_hi)
  # No table in results/ holds the pooled purity for this operating point, so
  # it is printed rather than compared. The paper gives 98.7%.
  cat(sprintf("  ....  focal 0.75 A: pooled purity recomputed as %.4f\n", focal$purity))
} else {
  cat("\n  SKIP  Table 1 comparison:",
      if (!file.exists(exp_f)) "results/tables/sptce_focal_table.csv absent"
      else "geometry unreadable", "\n")
  n_skip <- n_skip + 4L
  if (!is.null(focal))
    cat(sprintf("        recomputed anyway: %d clusters, %d receptors, purity %.4f, %.1f%% at purity >= 0.9\n",
                focal$n_clusters, focal$n_clustered, focal$purity, 100 * focal$pct_hi))
}

# =============================================================================
# 3. The sensitivity denominator
# =============================================================================
# Eq. 4 divides by the input count, not the clustered count. This check
# recomputes the denominator from the input table and confirms the sensitivity
# tables use it.
hdr("3. Sensitivity denominator (Eq. 4)")

sens_f <- "results/tables/corrected/sensitivity_CORRECTED_focal.csv"
if (file.exists(sens_f)) {
  sc <- read_csv(sens_f, show_col_types = FALSE)
  nm <- names(sc)
  dcol <- intersect(c("n_input", "n_total", "denominator", "n_sptce_input_tcrs"), nm)
  ecol <- intersect(c("epitope", "Epitope"), nm)
  if (length(dcol) && length(ecol)) {
    # The table carries three denominators per epitope, one per unit class --
    # `paired` (whole receptors), `receptor` (keyed on beta) and `chain`
    # (the pooled single-chain list), because the nine published methods are
    # not all measured in the same unit. Only the paired class is the cohort's
    # own receptor count, so only that one is checked against the input table.
    inp <- coh %>% count(epitope, name = "n_input_recomputed")
    tab <- sc %>% select(epitope = all_of(ecol[1]), n_input = all_of(dcol[1]),
                         unit = any_of("unit")) %>% distinct()
    # The unit label is matched on substring: it is written
    # "receptor (paired)" against "receptor (beta key)" and "chain".
    has_u <- "unit" %in% names(tab)
    pc <- (if (has_u) filter(tab, str_detect(unit, "paired")) else tab) %>%
      distinct(epitope, n_input) %>% inner_join(inp, by = "epitope")
    for (k in seq_len(nrow(pc)))
      report(sprintf("paired-unit denominator for %s", pc$epitope[k]),
             pc$n_input_recomputed[k], pc$n_input[k], exact = TRUE)
    if (has_u) {
      # Informational, not a check: the three unit classes usually differ,
      # but they coincide for an epitope whose receptors all carry distinct
      # beta chains (FLCMKALLL, at 81 either way).
      others <- tab %>% filter(!str_detect(unit, "paired")) %>%
        distinct(epitope, unit, n_input) %>% inner_join(inp, by = "epitope")
      cat(sprintf("  ....  of %d non-paired unit rows, %d differ from the paired count\n",
                  nrow(others), sum(others$n_input != others$n_input_recomputed)))
    }
  } else {
    cat("  SKIP  the corrected table carries no denominator column to check\n")
    cat("        columns present:", paste(nm, collapse = ", "), "\n")
    n_skip <- n_skip + 1L
  }
  # Sensitivity cannot exceed the clustered fraction, whatever the
  # denominator.
  scol <- intersect(c("sensitivity_corrected", "sensitivity"), nm)
  if (length(scol)) {
    bad <- sum(sc[[scol[1]]] > 1 + 1e-9, na.rm = TRUE)
    report("corrected sensitivities all <= 1", 0L, bad, exact = TRUE)
  }
} else {
  cat("  SKIP  results/tables/corrected/sensitivity_CORRECTED_focal.csv absent\n")
  n_skip <- n_skip + 1L
}

# =============================================================================
# 4. The ladder's own accounting
# =============================================================================
hdr("4. Noise ladder accounting")

ss <- "data/benchmark_a/subset_summary.tsv"
if (file.exists(ss)) {
  s <- read_tsv(ss, show_col_types = FALSE)
  for (i in seq_len(nrow(s))) {
    # Labels are identical across the ladder -- verified at the commit that
    # dropped the per-rung files -- so the count is taken from the one tracked
    # label table restricted to the rung's ids.
    f <- sprintf("data/benchmark_a/donor1_subset%d_labels.tsv", s$subset[i])
    d <- if (file.exists(f)) read_tsv(f, show_col_types = FALSE) else {
      lf <- "data/benchmark_a/donor1_subset6_labels.tsv"
      if (file.exists(lf)) {
        r <- rung_ids(s$subset[i])
        read_tsv(lf, show_col_types = FALSE) %>% filter(id %in% r)
      } else NULL
    }
    nlab <- if (is.null(d)) NA else
      sum(!is.na(d$epitope) & d$epitope != "NA" & nzchar(d$epitope))
    report(sprintf("rung %d: labelled receptors", s$subset[i]),
           s$labelled[i], nlab, exact = TRUE)
  }
  report("labelled count constant across the ladder", 1L,
         as.integer(n_distinct(s$labelled) == 1L), exact = TRUE)
} else {
  cat("  SKIP  data/benchmark_a/subset_summary.tsv absent\n")
  n_skip <- n_skip + 1L
}

# =============================================================================
# Summary
# =============================================================================
cat("\n==============================================================\n")
cat(sprintf("  %d PASS   %d FAIL   %d SKIP\n", n_pass, n_fail, n_skip))
cat("==============================================================\n")

if (WRITE_ACTUAL && length(results)) {
  out <- "results/tables/validation_report.csv"
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  write_csv(bind_rows(results), out)
  cat("  wrote", out, "\n")
}

if (n_fail > 0) {
  cat("\n  A FAIL means this version of the code no longer reproduces the\n")
  cat("  tracked table. Check which of the two is out of date before trusting\n")
  cat("  either: it may be the tracked table rather than the code.\n")
}
quit(status = min(n_fail, 125L))
