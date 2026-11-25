# =============================================================================
# SPACE-TCR — metric definitions
# Devinder Kaur Srai, Oxford DPhil 2025
#
#   source("sptce_metrics.R")
#
# The one place in this repository where clusters are built and scored.
# sptce_pipeline.R, validate_outputs.R, the scripts under scripts/ and the R
# tests all source this file; none of them carries its own copy of a metric.
#
# Written directly against the reference implementation in the Jouannet
# repository, so SPACE-TCR and the published methods are scored by the same
# equations on the same units:
#   Functions/purity_functions.R
#   Analysis/benchmark_tools_pooled_database_CD8.Rmd
#
# Everything here is a function of its arguments plus the METRIC constants
# below. Paths, cohorts and operating points live in CFG in sptce_pipeline.R.
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(tibble)
  library(stringdist)
})

# dplyr verbs are bound locally so that an attached seqinr or Bioconductor S4
# stack in the calling session cannot mask them ("Argument 'x' is not a vector:
# list" from seqinr::count).
count      <- dplyr::count
filter     <- dplyr::filter
select     <- dplyr::select
rename     <- dplyr::rename
slice      <- dplyr::slice
first      <- dplyr::first
n_distinct <- dplyr::n_distinct
coalesce   <- dplyr::coalesce
desc       <- dplyr::desc
intersect  <- base::intersect
setdiff    <- base::setdiff
union      <- base::union

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (length(a) == 1 && is.na(a))   return(b)
  a
}

# =============================================================================
# 1. METRIC CONSTANTS — as in the reference implementation
# =============================================================================
METRIC <- list(
  # Reference: Rmd line 229, `th <- 1`, used as `filter(size > th)`. A cluster
  # therefore qualifies at size >= 2.
  min_cluster_size = 2L,
  # Reference: Rmd line 351, `purity >= 0.9`. Not `> 0.9`.
  high_purity_cut  = 0.9,
  # Reference: Rmd line 1547, `size > 3 & perc_epi_in_cluster >= 0.67`. The
  # literal 0.67, not 2/3: 2/3 would admit 4-of-6, 6-of-9 and 8-of-12 clusters
  # that the reference excludes.
  sens_min_size    = 4L,     # size > 3
  sens_min_frac    = 0.67,
  # IMGT positions of CDR3 whose phi/psi angles make up the dihedral vector.
  cdr3_imgt        = 105:117,
  # Median within-cluster CDR3 Levenshtein distance class breaks: <=1, <=3, <=6, >6.
  ld_breaks        = c(1, 3, 6)
)

# =============================================================================
# 2. CLUSTER SCORING
# =============================================================================
# All of these take a long "assignments" table with at least
#   cluster_id, <unit>, epitope
# and are agnostic to how the clusters were produced. `unit` names the column
# that defines receptor identity for this method: a receptor id for SPACE-TCR
# and other paired methods, `cdr3b` for single-chain methods.

# Per-cluster size, dominant count and purity.
#
# A cluster's size is its number of DISTINCT units. Its dominant count is the
# max over labels of the number of distinct units carrying that label. A unit
# with several labels (a CDR3b annotated with two epitopes) counts towards each
# of them, as in the reference: purity_functions.R merges every annotation,
# applies `unique()`, then takes the per-cluster maximum.
#
# `label` is the column scored: `epitope` for purity, a gene column for gene
# homogeneity.
cluster_table <- function(asg, unit = "receptor_id", min_size = METRIC$min_cluster_size,
                          label = "epitope") {
  stopifnot(all(c("cluster_id", unit, label) %in% names(asg)))
  d <- asg %>% filter(!is.na(.data[[unit]]))
  sizes <- d %>% distinct(cluster_id, .data[[unit]]) %>% count(cluster_id, name = "size")
  per_lab <- d %>%
    filter(!is.na(.data[[label]])) %>%
    distinct(cluster_id, .data[[unit]], .data[[label]]) %>%
    count(cluster_id, .lab = .data[[label]], name = "n_lab")
  dom <- if (nrow(per_lab)) per_lab %>%
    group_by(cluster_id) %>%
    summarise(n_dom       = max(n_lab),
              n_labels    = n_distinct(.lab),
              dom_epitope = .lab[which.max(n_lab)],
              .groups = "drop")
  else tibble(cluster_id = character(), n_dom = integer(), n_labels = integer(),
              dom_epitope = character())
  sizes %>%
    left_join(dom, by = "cluster_id") %>%
    mutate(n_dom    = coalesce(n_dom, 0L),
           n_labels = coalesce(n_labels, 0L),
           purity   = ifelse(size > 0, n_dom / size, NA_real_)) %>%
    filter(size >= min_size)
}

# Global purity. Reference: purity_functions.R `purity_function` —
# sum(nb_majoritaire) / sum(size) over clusters above the size threshold.
jm_purity <- function(ct) {
  if (!nrow(ct)) return(NA_real_)
  sum(ct$n_dom) / sum(ct$size)
}

# Polyspecificity: fraction of non-singleton clusters carrying more than one
# epitope label.
jm_polyspecificity <- function(ct) {
  if (!nrow(ct)) return(NA_real_)
  mean(ct$n_labels > 1)
}

# The headline metric row for one clustering.
#
#   retention        distinct units in a qualifying cluster / denom
#   purity           jm_purity()
#   pct_clusters_hi  share of clusters with purity >= 0.9
#   pct_units_hi     distinct units in a high-purity cluster / distinct units clustered
#   pct_members_hi   cluster memberships in high-purity clusters / all memberships;
#                    equal to pct_units_hi unless a unit sits in several clusters
#
# `denom` is supplied by the caller and never inferred from the clustering.
jm_metrics <- function(asg, denom, unit = "receptor_id",
                       min_size = METRIC$min_cluster_size,
                       hi = METRIC$high_purity_cut) {
  ct <- cluster_table(asg, unit = unit, min_size = min_size)
  if (!nrow(ct)) {
    return(tibble(n_clusters = 0L, n_units_clustered = 0L, n_rows_clustered = 0L,
                  retention = 0, purity = NA_real_,
                  pct_clusters_hi = NA_real_, pct_units_hi = NA_real_,
                  pct_members_hi = NA_real_, polyspecificity = NA_real_,
                  n_memberships = 0L, unit = unit, denom = denom))
  }
  kept <- asg %>% filter(cluster_id %in% ct$cluster_id, !is.na(.data[[unit]]))
  n_units <- n_distinct(kept[[unit]])
  hi_ct <- filter(ct, purity >= hi)
  hi_units <- kept %>% filter(cluster_id %in% hi_ct$cluster_id) %>%
    pull(.data[[unit]]) %>% n_distinct()
  tibble(
    n_clusters        = nrow(ct),
    n_units_clustered = n_units,
    n_rows_clustered  = nrow(kept),      # reported separately; never a denominator
    retention         = n_units / denom,
    purity            = jm_purity(ct),
    pct_clusters_hi   = nrow(hi_ct) / nrow(ct),
    pct_units_hi      = if (n_units > 0) hi_units / n_units else NA_real_,
    pct_members_hi    = sum(hi_ct$size) / sum(ct$size),
    polyspecificity   = jm_polyspecificity(ct),
    n_memberships     = sum(ct$size),
    unit              = unit,
    denom             = denom
  )
}

# Epitope-level sensitivity. Reference: Rmd lines 1547 and 1575/1578.
#
#   qualifying cluster c for epitope e:  n_{c,e} equals c's dominant count
#                                        size(c) > 3
#                                        n_{c,e} / size(c) >= 0.67
#   numerator:   distinct units labelled e in those clusters
#   denominator: units labelled e in the INPUT, not in the clustering
#
# The denominator must be supplied: dividing by the clustered count would
# reward a method for clustering less. Returns NA when it is missing.
jm_sensitivity <- function(asg, epi, n_input_e, unit = "receptor_id") {
  jm_sensitivity_detail(asg, epi, n_input_e, unit)$sensitivity
}

# The same, returning the recovered count alongside the ratio.
jm_sensitivity_detail <- function(asg, epi, n_input_e, unit = "receptor_id") {
  if (is.null(n_input_e) || !length(n_input_e) || is.na(n_input_e) || n_input_e <= 0)
    return(list(n_recovered = NA_integer_, sensitivity = NA_real_))
  ct <- cluster_table(asg, unit = unit)
  if (!nrow(ct)) return(list(n_recovered = 0L, sensitivity = 0))
  per_epi <- asg %>%
    filter(!is.na(.data[[unit]]), cluster_id %in% ct$cluster_id,
           !is.na(epitope), epitope == epi) %>%
    distinct(cluster_id, .data[[unit]]) %>%
    count(cluster_id, name = "n_e")
  qual <- ct %>%
    inner_join(per_epi, by = "cluster_id") %>%
    filter(n_e == n_dom,
           size >= METRIC$sens_min_size,
           n_e / size >= METRIC$sens_min_frac)
  if (!nrow(qual)) return(list(n_recovered = 0L, sensitivity = 0))
  recovered <- asg %>%
    filter(cluster_id %in% qual$cluster_id, !is.na(epitope), epitope == epi,
           !is.na(.data[[unit]])) %>%
    pull(.data[[unit]]) %>% n_distinct()
  list(n_recovered = recovered, sensitivity = recovered / n_input_e)
}

# Germline composition. Called homogeneity, not purity: a cluster's V/J makeup
# has no ground truth, only self-consistency. Same counting as cluster_table,
# over the units that carry a gene call.
jm_gene_homogeneity <- function(asg, gene_col, unit = "receptor_id",
                                hi = METRIC$high_purity_cut) {
  if (!gene_col %in% names(asg)) return(tibble())
  d <- asg %>% filter(!is.na(.data[[gene_col]]), !is.na(.data[[unit]]))
  if (!nrow(d)) return(tibble())
  per <- cluster_table(d, unit = unit, label = gene_col)
  if (!nrow(per)) return(tibble())
  tibble(gene = gene_col,
         n_clusters_scored   = nrow(per),
         purity              = jm_purity(per),
         pct_clusters_hi     = mean(per$purity >= hi),
         pct_members_hi      = sum(per$size[per$purity >= hi]) / sum(per$size),
         mean_distinct_genes = mean(per$n_labels),
         n_memberships       = sum(per$size))
}

# Purity of a PUBLISHED method's own clustering, exactly as purity_functions.R
# scores it. Differs from cluster_table() in the two ways the reference does:
#   - `size` is the method's own cluster size (its output file's count), not a
#     recount of labelled units;
#   - a cluster none of whose members carries an annotation is dropped, because
#     the reference inner-joins against the annotation table before counting.
#
#   out          one row per (key, clusters), optionally with `size`
#   annotations  key -> Epitope, one row per annotation (multi-epitope allowed)
#
# Returns one row per scored cluster: clusters, size, nb_maj, purity.
jm_reference_cluster_scores <- function(out, annotations, key) {
  if (!"size" %in% names(out))
    out <- out %>% group_by(clusters) %>% mutate(size = n()) %>% ungroup()
  big <- filter(out, size >= METRIC$min_cluster_size)
  if (!nrow(big)) return(tibble(clusters = character(), size = integer(),
                                nb_maj = integer(), purity = numeric()))
  m <- big %>% inner_join(annotations, by = key) %>% distinct()
  pc <- m %>% count(clusters, Epitope, name = "nb") %>%
    group_by(clusters) %>% summarise(nb_maj = max(nb), .groups = "drop")
  left_join(pc, distinct(big, clusters, size), by = "clusters") %>%
    mutate(purity = nb_maj / size)
}

# The four headline numbers for a published method, from
# jm_reference_cluster_scores(). `pct_tcrs` picks the paper's per-method
# denominator for "% TCRs in high-purity clusters":
#   "size"  -> share of summed cluster sizes          (single-chain, clusTCR, DeepTCR)
#   "units" -> share of distinct annotated keys       (GLIPH2, TCRdist3)
# Returns c(retention, purity, pct_clusters_hi, pct_tcrs_hi).
jm_reference_metrics <- function(out, annotations, key, retention,
                                 pct_tcrs = c("size", "units"),
                                 hi = METRIC$high_purity_cut) {
  pct_tcrs <- match.arg(pct_tcrs)
  s <- jm_reference_cluster_scores(out, annotations, key)
  if (!nrow(s)) return(c(retention, NA, NA, NA))
  if (pct_tcrs == "units") {
    if (!"size" %in% names(out))
      out <- out %>% group_by(clusters) %>% mutate(size = n()) %>% ungroup()
    m <- filter(out, size >= METRIC$min_cluster_size) %>%
      inner_join(annotations, by = key) %>% distinct()
    hi_cl <- s$clusters[s$purity >= hi]
    perc <- (m %>% filter(clusters %in% hi_cl) %>% distinct(.data[[key]]) %>% nrow()) /
            (m %>% distinct(.data[[key]]) %>% nrow())
  } else {
    perc <- sum(s$size[s$purity >= hi]) / sum(s$size)
  }
  c(retention, sum(s$nb_maj) / sum(s$size), mean(s$purity >= hi), perc)
}

# Pairwise CDR3 Levenshtein statistics within one cluster. Clusters larger
# than `cap` are subsampled with a fixed seed so the result is reproducible.
ld_stats <- function(seqs, cap = 300, seed = 42L) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs) > 0]
  if (length(seqs) < 2)
    return(tibble(mean_ld = NA, median_ld = NA, pct_le1 = NA, pct_gt3 = NA))
  if (length(seqs) > cap) { set.seed(seed); seqs <- sample(seqs, cap) }
  v <- stringdistmatrix(seqs, seqs, method = "lv")[
    lower.tri(matrix(0, length(seqs), length(seqs)))]
  tibble(mean_ld = mean(v), median_ld = median(v),
         pct_le1 = mean(v <= 1), pct_gt3 = mean(v > 3))
}

# Median within-cluster CDR3 Levenshtein distance to a divergence class.
LD_CLASSES <- c("Sequence-identical / near-identical", "Sequence-similar",
                "Sequence-diverse", "Sequence-unrelated")
ld_class_of <- function(x, labels = LD_CLASSES) {
  b <- METRIC$ld_breaks
  factor(case_when(x <= b[1] ~ labels[1], x <= b[2] ~ labels[2],
                   x <= b[3] ~ labels[3], !is.na(x) ~ labels[4]),
         levels = labels)
}

# =============================================================================
# 3. SELF-CHECKS
# =============================================================================
check_metrics <- function(m, label = "") {
  bad <- character()
  if (any(m$n_units_clustered > m$denom, na.rm = TRUE))
    bad <- c(bad, "units clustered exceeds the input denominator (modes pooled?)")
  if (any(m$retention > 1, na.rm = TRUE))
    bad <- c(bad, "retention above 1")
  if (any(m$purity > 1 + 1e-9, na.rm = TRUE))
    bad <- c(bad, "purity above 1")
  if (any(m$pct_clusters_hi > 1 + 1e-9, na.rm = TRUE))
    bad <- c(bad, "high-purity cluster share above 1")
  if (length(bad))
    stop("check_metrics failed", if (nzchar(label)) paste0(" [", label, "]") else "",
         ":\n  - ", paste(bad, collapse = "\n  - "), call. = FALSE)
  invisible(TRUE)
}

# Sensitivity cannot exceed the fraction of that epitope's input receptors that
# were clustered at all.
check_sensitivity <- function(sens_df, label = "") {
  bad <- sens_df %>% filter(!is.na(sensitivity), !is.na(clustered_frac),
                            sensitivity > clustered_frac + 1e-9)
  if (nrow(bad)) {
    print(as.data.frame(bad), row.names = FALSE)
    stop("check_sensitivity failed", if (nzchar(label)) paste0(" [", label, "]") else "",
         ": sensitivity exceeds the clustered fraction for the rows above. ",
         "This is the signature of pooled modes or row-counted numerators.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# =============================================================================
# 4. CLUSTERING
# =============================================================================
# Receptor id from a PDB path: "/x/y/tcr12.pdb" -> "tcr12".
clean_id <- function(x) sub("\\.pdb(\\.gz)?$", "", basename(as.character(x)))

# Union-find: component label for each of n nodes given edges (ei, ej).
# Everything is forced to integer: one double index would promote the whole
# parent vector to double and break the closing vapply(..., integer(1)).
connected_components <- function(n, ei, ej) {
  n <- as.integer(n)
  ei <- as.integer(ei); ej <- as.integer(ej)
  p <- seq_len(n)
  root <- function(x) {
    x <- as.integer(x)
    while (p[x] != x) x <- p[x]
    x
  }
  for (k in seq_along(ei)) {
    ri <- root(ei[k]); rj <- root(ej[k])
    if (ri != rj) p[ri] <- rj
  }
  vapply(seq_len(n), root, integer(1))
}

# Mean circular distance between two dihedral-angle vectors, in RADIANS, in
# [0, pi]. The shorter vector is padded with pi: truncating would compare only
# the shared prefix of two loops of different length and could call them
# identical. scripts/dihedral_util.py is the Python twin used on the cluster;
# tests/test_dihedral_parity.R checks the two agree.
circ_dihedral <- function(a, b) {
  if (!length(a) || !length(b)) return(pi)
  n <- max(length(a), length(b))
  a <- c(a, rep(pi, n - length(a))); b <- c(b, rep(pi, n - length(b)))
  d <- abs(a - b) %% (2 * pi)
  mean(pmin(d, 2 * pi - d))
}

# tcr_id -> c(phi..., psi...) over the CDR3 IMGT positions, from a Stage D
# dihedrals table (path, or an already-read data frame).
load_dihedral_map <- function(x) {
  if (is.character(x)) {
    if (!length(x) || is.na(x) || !file.exists(x)) return(NULL)
    x <- read_tsv(x, show_col_types = FALSE)
  }
  dv <- x %>%
    filter(resi %in% METRIC$cdr3_imgt, !is.na(phi_rad), !is.na(psi_rad)) %>%
    arrange(tcr_id, chain, resi) %>%
    group_by(tcr_id) %>%
    summarise(angles = list(c(phi_rad, psi_rad)), .groups = "drop")
  setNames(dv$angles, dv$tcr_id)
}

# SPACE-TCR clusters from a pair table: connected components of the graph
# joining receptors within `rmsd_cut` (and, when gated, within `dih_cut`
# radians of dihedral distance), computed per `bin` (one CDR length bin each).
# Clusters smaller than METRIC$min_cluster_size are dropped.
#
#   pairs     columns bin, ti, tj, rmsd, and dih when gating (NA = no angles;
#             such a pair is never within the gate)
#   labels    one row per receptor, keyed by `receptor_id`; joined onto members
#   keep_ids  if given, only pairs with both ends in it are used, BEFORE
#             components are formed
#   id_fmt    function(bin, comp) giving the cluster id
#
# Returns one row per (receptor_id, cluster_id) with the label columns.
cluster_pairs <- function(pairs, labels, rmsd_cut, dih_cut = NA_real_,
                          keep_ids = NULL, keep_unlabelled = FALSE,
                          id_fmt = function(bin, comp) paste0(bin, "_", comp)) {
  dt <- pairs %>% filter(!is.na(rmsd), rmsd <= rmsd_cut)
  if (!is.null(keep_ids)) dt <- filter(dt, ti %in% keep_ids, tj %in% keep_ids)
  if (!is.na(dih_cut)) dt <- filter(dt, !is.na(dih), dih <= dih_cut)
  if (!nrow(dt)) return(tibble())
  rows <- dt %>% group_split(bin) %>% map_dfr(function(g) {
    ids  <- unique(c(g$ti, g$tj))
    idx  <- setNames(seq_along(ids), ids)
    comp <- connected_components(length(ids), idx[g$ti], idx[g$tj])
    tibble(receptor_id = ids, cluster_id = id_fmt(g$bin[1], comp))
  })
  out <- rows %>%
    distinct(receptor_id, cluster_id) %>%
    group_by(cluster_id) %>% filter(n() >= METRIC$min_cluster_size) %>% ungroup() %>%
    left_join(labels, by = "receptor_id")
  if (!keep_unlabelled) out <- filter(out, !is.na(epitope))
  out
}

# Read Stage C distance tables into a pair table for cluster_pairs(): one bin
# per file, receptor ids from the PDB paths, failed pairs dropped. With a
# dihedral map, `dih` is each pair's circ_dihedral() distance.
read_pair_tables <- function(dist_files, max_rmsd = Inf, dih_map = NULL) {
  map_dfr(dist_files, function(f) {
    dt <- tryCatch(read_tsv(f, show_col_types = FALSE, progress = FALSE),
                   error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    if (all(c("id_a", "id_b") %in% names(dt)) && !"id_i" %in% names(dt))
      dt <- rename(dt, id_i = id_a, id_j = id_b)
    if (!all(c("id_i", "id_j", "rmsd") %in% names(dt))) return(NULL)
    if ("status" %in% names(dt)) dt <- filter(dt, status == "ok")
    dt <- filter(dt, !is.na(rmsd), rmsd <= max_rmsd)
    if (!nrow(dt)) return(NULL)
    out <- tibble(bin = basename(f), ti = clean_id(dt$id_i), tj = clean_id(dt$id_j),
                  rmsd = as.numeric(dt$rmsd))
    if (!is.null(dih_map))
      out$dih <- vapply(seq_len(nrow(out)), function(k) {
        ai <- dih_map[[out$ti[k]]]; aj <- dih_map[[out$tj[k]]]
        if (is.null(ai) || is.null(aj)) NA_real_ else circ_dihedral(ai, aj)
      }, numeric(1))
    out
  })
}

# SPACE-TCR clusters at one operating point, straight from the distance-table
# files: read_pair_tables() then cluster_pairs().
#
#   prefix   NULL -> cluster ids "<file>|<rmsd>|<component>";
#            otherwise "<prefix>_<file>_<component>"
#
# Returns one row per (receptor_id, cluster_id) with the label columns, plus
# `method`, `threshold` and `dih` (NA for the RMSD-only mode).
cluster_at <- function(dist_files, labels, rmsd_cut, dih_cut = NA_real_,
                       dih_map = NULL, keep_unlabelled = FALSE, prefix = NULL) {
  if (!is.na(dih_cut) && is.null(dih_map))
    stop("dih_cut given but no dihedral map loaded", call. = FALSE)
  pairs <- read_pair_tables(dist_files, rmsd_cut, if (is.na(dih_cut)) NULL else dih_map)
  if (!nrow(pairs)) return(tibble())
  fmt <- if (is.null(prefix)) function(bin, comp) sprintf("%s|%s|%s", bin, rmsd_cut, comp)
         else function(bin, comp) paste0(prefix, "_", bin, "_", comp)
  out <- cluster_pairs(pairs, labels, rmsd_cut, dih_cut,
                       keep_unlabelled = keep_unlabelled, id_fmt = fmt)
  if (!nrow(out)) return(out)
  out %>% mutate(method = if (is.na(dih_cut)) "SPACE-TCR RMSD" else "SPACE-TCR RMSD+dih",
                 threshold = rmsd_cut, dih = dih_cut)
}
# End of file.
