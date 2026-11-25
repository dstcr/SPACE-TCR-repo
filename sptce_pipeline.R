# =================================================================
# SPACE-TCR pipeline — full pipeline and every figure
# Devinder Kaur Srai, Oxford DPhil 2025
#
#   source("sptce_pipeline.R")
#   run_benchmark_a()    # paired-only 10X noise-ladder comparison
#   run_benchmark_b()    # Jouannet curated-cohort benchmark + extensions
#   main()               # A + B
#
# Conventions
#   - V/J gene composition is called "consensus", not "purity": purity is scored
#     against a ground truth, while germline composition can only be
#     self-consistent. For the epitope label the same arithmetic IS purity.
#   - High purity is ">= 0.9".
#   - Cohort accounting is emitted per SPACE-TCR mode, never pooled across
#     operating points, with a self-check that fails loudly.
#   - Outputs include sptce_dihedral_mode_metrics.csv, focal-epitope
#     contamination composition, STCRDab leakage on both denominators, a
#     networks-overview check and gene-consensus denominators.
#
# tests/test_pipeline.R checks that every stage function is defined.
# run_stage10_huardb is defined but not called from main(): Benchmark C needs a
# HuARdb download that is not part of this repository.
# =================================================================

suppressPackageStartupMessages({
  library(readr);   library(dplyr);     library(tidyr);  library(purrr)
  library(stringr); library(tibble);    library(forcats)
  library(ggplot2); library(patchwork); library(scales)
  library(stringdist); library(tools);  library(reshape2)
  try(library(igraph),    silent = TRUE)
  try(library(ggraph),    silent = TRUE)
  try(library(tidygraph), silent = TRUE)
  try(library(cowplot),   silent = TRUE)
  try(library(ggseqlogo), silent = TRUE)
  try(library(kableExtra),silent = TRUE)
  try(library(ggforce),   silent = TRUE)   # .plot_networks_overview
  try(library(ggrepel),   silent = TRUE)   # consensus / benchmark-A labels
})

# Every metric -- clustering, purity, retention, sensitivity, the dihedral
# distance -- comes from sptce_metrics.R. Nothing in this file scores clusters
# itself. Run from the repository root.
if (!file.exists("sptce_metrics.R"))
  stop("sptce_metrics.R not found: run from the repository root", call. = FALSE)
source("sptce_metrics.R", local = TRUE)

# Guard against namespace masking in the calling session. Packages such as seqinr
# (k-mer count()) and the Bioconductor S4 stack (S4Vectors/IRanges first(), slice(),
# rename(), ...) export functions that mask the dplyr/base versions used here, e.g.
# "Argument 'x' is not a vector: list" (count) or "unable to find an inherited method
# for function 'first'" (S4 first). Binding the intended versions in this environment
# shadows any attached-package versions for all unqualified calls below, regardless of
# package load order. dplyr verbs -> dplyr; plain vector set-ops -> base.
count      <- dplyr::count
filter     <- dplyr::filter
select     <- dplyr::select
rename     <- dplyr::rename
first      <- dplyr::first
slice      <- dplyr::slice
n_distinct <- dplyr::n_distinct
coalesce   <- dplyr::coalesce
desc       <- dplyr::desc
intersect  <- base::intersect
setdiff    <- base::setdiff

# =================================================================
# 1. NULL-COALESCE  (defined before CONFIG, which uses it)
# =================================================================
# Treats NULL, length-0 AND a scalar NA as missing. CFG lookups can be absent or
# explicitly NA and the extensions rely on the NA behaviour, so this is the
# union of both forms. Defining it here rather than among the utilities matters:
# compute_published_metrics_exact() and the Benchmark A helpers call it.
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (length(a) == 1 && is.na(a))   return(b)
  a
}

# =================================================================
# 2. CONFIG
# =================================================================
# Two kinds of path:
#
#   IN THE REPOSITORY  Stage A inputs and Stage C/D geometry for both
#                      benchmarks (see data/README.md). The defaults point
#                      inside this repository, so it runs without editing.
#
#   EXTERNAL           10X raw tables, the Jouannet repository and the STCRDab
#                      parse. The default is <SPTCE_EXTERNAL>/<filename>, with
#                      SPTCE_EXTERNAL an environment variable defaulting to
#                      ~/sptce_external; individual entries can also be set
#                      after sourcing.
EXT <- Sys.getenv("SPTCE_EXTERNAL", unset = path.expand("~/sptce_external"))
.ext <- function(f) file.path(EXT, f)

CFG <- list(
  # ---- EXTERNAL: 10X donor-1 raw tables (Benchmark A noise fraction) ------
  contig_file           = .ext("vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv"),
  binmat_file           = .ext("vdj_v1_hs_aggregated_donor1_binarized_matrix.csv"),
  master_file           = .ext("donor1_paired_master.tsv"),
  sptce_input           = .ext("donor1_sptce_input.tsv"),

  # ---- TRACKED: Benchmark A ladder inputs + reference structural run ------
  # data/benchmark_a holds donor1_subset{1..6}_sptce.tsv (Stage A input),
  # the matching label tables, and subset4_clusters/ + subset4_dihedrals/,
  # which are the Stage C and Stage D output every rung is derived from.
  subset_dir            = "data/benchmark_a",
  sptce_root            = "data/benchmark_a",

  analysis_dir          = "results/figures/benchmark_a",
  jouannet_metric_table = .ext("metric_table.csv"),
  compare_dir           = "results/tables/sptce_vs_jouannet",

  # ---- TRACKED: Benchmark B cohort input + structural run -----------------
  # data/cohort_structural holds distance_tables/*.pairwise.tsv.gz,
  # dihedrals.tsv.gz, cluster_assignments.tsv and clusters/member_to_bin.tsv
  # for the 4,779-receptor curated cohort. Every threshold in Tables 1 and 2
  # is re-derivable from these without rerunning TCRBuilder2.
  jouannet_sptce_out    = "data/cohort_structural",
  jouannet_labelled     = "data/jouannet_labelled_input.tsv",
  jouannet_out_dir      = "results/tables/cohort",

  # ---- EXTERNAL: the Jouannet benchmark clone ------------------------------
  # third_party/fetch_jouannet.sh clones it; it is needed only for the nine
  # published methods' own cluster output, not for anything SPACE-TCR computes.
  jouannet_repo         = .ext("TCR_Unsupervised_Benchmark-main"),
  jouannet_cell_type    = "CD8",
  
  # ---- EXTERNAL: Stage 7, training-data leakage ---------------------------
  # Needs the parsed IgBLAST output for the STCRDab reference set.
  stcrdab_parsed_csv    = .ext("igblast_tcr_parsed.csv"),
  stcrdab_pdb_table     = .ext("PDB_Database.txt"),  # legacy (unused)
  stcrdab_full_table    = .ext("opig.csv"),          # legacy (unused)
  stcrdab_match_mode    = "cdr3_exact",
  # NULL => <jouannet_out_dir>/stcrdab_overlap, derived below. It was an
  # independent absolute path, so overriding jouannet_out_dir moved the figure
  # but not Stage 7's three CSVs, which then landed in a literal
  # "C:/Users/..." tree created under the working directory on Linux.
  stcrdab_out_dir       = NULL,
  
  # ---- Benchmark A: the Jouannet 10X donor-1 subsets ----------------------
  tenx_root             = NULL,   # <repo>/Data/10X/Donor1
  tenx_out              = NULL,   # <repo>/Data/10X_Output_methods
  # Path to the reference run's dihedral table, if it is not at
  # <sptce_root>/subset<N>_dihedrals/dihedrals.tsv.gz. Benchmark A adds the
  # dihedral-gated modes only when this resolves. It must index the donor-1
  # receptor ids, not Benchmark B's.
  benchA_dih_path       = NULL,
  
  # The ladder has six rungs. Note that nothing in the pipeline reads this
  # value: jouannet_10x_ladder() reads the rungs from the repo directly, and
  # .detect_ref_subset() discovers the structural runs on disk. Kept for
  # reference and for external scripts.
  subsets           = 1:6,
  thresholds        = c(0.25, 0.50, 0.75, 1.00, 1.25),
  focal_threshold   = 0.75,
  # The dihedral gate, in radians, used by run_stage5 for the gated clusters.
  #
  # 0.5 is the value in the paper ("unless otherwise stated this angle was set
  # at a threshold of 0.5-radian or 28 degrees"), so the default run reproduces
  # the published tables.
  #
  # For comparison: at 0.5 and 0.75 A a 0.5 rad gate moves every reported
  # metric by under a percentage point, and on the noise ladder it changes
  # purity under full dilution by <= 0.013; 0.2 rad reduces retention to 5-7.5%
  # at every threshold; 0.3 rad separates visibly while keeping usable
  # coverage. No gated setting reaches the purity of an RMSD-only threshold of
  # comparable retention.
  #
  # Changing it changes every RMSD+dih row in Tables 1 and 2. Output files are
  # prefixed with the value (d0.5_ / d0.3_).
  focal_dihedral    = 0.5,
  dihedral_modes    = c(0.2, 0.3, 0.5),
  dihedral_sweep    = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0),
  rmsd_sweep_fixed  = 1.0,
  focal_epitope     = "GILGFVFTL",
  extra_combos      = list(
    list(thr = 1.00, dih = 0.2),
    list(thr = 1.25, dih = 0.2),
    list(thr = 1.00, dih = 0.3)
  ),
  
  jouannet_epitopes = c("GILGFVFTL","NLVPMVATV","GLCTLVAML",
                        "FLCMKALLL","NYNYLYRLF","LLFGYPVYV"),
  
  jouannet_method_order = c(
    "Hamming_Distance","Levenshtein_Distance","TCRMatch",
    "iSMART","GIANA","clusTCR","GLIPH2","DeepTCR","TCRdist3"),
  jouannet_method_colour = c(
    Hamming_Distance     = "darkorange",  Levenshtein_Distance = "cyan3",
    TCRMatch             = "darkorchid4", iSMART               = "grey",
    GIANA                = "brown1",      clusTCR              = "darkgoldenrod1",
    GLIPH2               = "pink",        DeepTCR              = "chartreuse3",
    TCRdist3             = "black"),
  
  col_method    = "methods",
  col_subset    = "subset",
  col_retention = "retention",
  col_purity    = "purity",
  
  # Metric definitions (cluster size, purity cut, sensitivity rule, CDR3
  # positions) are METRIC in sptce_metrics.R, not here.
  # Per-epitope INPUT counts for Eq. 4's denominator, filled by
  # .set_sens_input_counts() from the labelled input table. NULL means every
  # sensitivity comes back NA, which is deliberate: see compute_sensitivity().
  .sens_input_counts = NULL,

  seed = 42L
)

# Derived config (cannot refer to itself inside the list() above). Resolved
# when called, so that setting CFG$jouannet_out_dir after sourcing takes
# effect.
.stcrdab_dir <- function() {
  d <- CFG$stcrdab_out_dir
  if (is.null(d) || !length(d) || is.na(d) || !nzchar(d))
    d <- file.path(CFG$jouannet_out_dir, "stcrdab_overlap")
  d
}
CFG$tenx_root <- CFG$tenx_root %||% file.path(CFG$jouannet_repo, "Data", "10X", "Donor1")
CFG$tenx_out  <- CFG$tenx_out  %||% file.path(CFG$jouannet_repo, "Data", "10X_Output_methods")
# Raw published-method outputs for the curated cohort, used by the CDR3a
# within-cluster LD analysis to read each paired method's OWN alpha chains
# rather than approximating them from a cdr3b -> modal cdr3a map.
CFG$method_outputs_dir <- file.path(CFG$jouannet_repo, "Data/Output_methods")

set.seed(CFG$seed)

# Output directories are created at load time only when their parent exists,
# so sourcing this file never creates folders for paths that do not apply on
# this machine. The stages that write output create their own directories.
.ensure_dir <- function(d) {
  if (is.null(d) || length(d) != 1 || is.na(d) || !nzchar(d)) return(invisible(FALSE))
  if (dir.exists(d)) return(invisible(TRUE))
  if (!dir.exists(dirname(d))) return(invisible(FALSE))   # not this machine
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  invisible(dir.exists(d))
}
for (d in c(CFG$subset_dir, CFG$analysis_dir, CFG$compare_dir, CFG$jouannet_out_dir))
  .ensure_dir(d)

# =================================================================
# 3. DISPLAY ORDER AND COLOUR PALETTES
# =================================================================
.sptce_focal_rmsd <- sprintf("SPACE-TCR RMSD (%sA)",     CFG$focal_threshold)
.sptce_focal_dih  <- sprintf("SPACE-TCR RMSD+dih (%sA)", CFG$focal_threshold)
# RMSD+dih is displayed in the bubble (fig10) at this looser RMSD cut, NOT 0.75:
.sptce_bubble_dih_thr <- 1.00
.sptce_bubble_dih     <- sprintf("SPACE-TCR RMSD+dih (%sA)", .sptce_bubble_dih_thr)

# The focal SPACE-TCR operating points reported across stage 5 (bubble, metric
# table, LD, polyspecificity, contamination, radar):
#   RMSD      at 0.75 A  (CFG$focal_threshold)
#   RMSD+dih  at 1.0  A  (.sptce_bubble_dih_thr)  -- looser RMSD, dihedral-gated
.sptce_focal_modes <- data.frame(
  mode      = c("SPACE-TCR RMSD",    "SPACE-TCR RMSD+dih"),
  threshold = c(CFG$focal_threshold, .sptce_bubble_dih_thr),
  stringsAsFactors = FALSE
)
filter_sptce_focal <- function(df)
  dplyr::inner_join(df, .sptce_focal_modes, by = c("mode", "threshold"))
# Like filter_sptce_focal but also retains "published"/NA-threshold rows, for
# figures overlaying published methods with the focal SPACE-TCR points.
keep_focal_rows <- function(df) dplyr::filter(df,
                                              is.na(threshold) | mode == "published" |
                                                (mode == "SPACE-TCR RMSD"     & threshold == CFG$focal_threshold) |
                                                (mode == "SPACE-TCR RMSD+dih" & threshold == .sptce_bubble_dih_thr))

.sptce_extras <- vapply(CFG$extra_combos,
                        \(c) sprintf("SPACE-TCR RMSD+dih<%srad (%sA)", c$dih, c$thr),
                        character(1))

METHOD_ORDER <- c(
  .sptce_focal_rmsd, .sptce_focal_dih, .sptce_bubble_dih, .sptce_extras,
  "SPACE-TCR", "SPACE-TCR (RMSD)",
  paste0("SPACE-TCR (RMSD+dih<", CFG$dihedral_modes, "rad)"),
  "clusTCR",
  setdiff(CFG$jouannet_method_order, "clusTCR")
)

order_methods <- function(x) {
  x <- as.character(x)
  factor(x, levels = c(intersect(METHOD_ORDER, x),
                       setdiff(unique(x), METHOD_ORDER)))
}

name_sptce_rmsd      <- \(thr)      sprintf("SPACE-TCR RMSD (%sA)", thr)
name_sptce_focal_dih <- \(thr)      sprintf("SPACE-TCR RMSD+dih (%sA)", thr)
name_sptce_mode      <- \(thr, dih) sprintf("SPACE-TCR RMSD+dih<%srad (%sA)", dih, thr)

# Node-link graph figures (stage-5b networks, figs 22b/23/25-degree, and the
# fig39 focal-epitope overlays) must reflect ONE operating point.
# sptce_focal_clusters bundles BOTH focal modes -- RMSD @ 0.75 A and RMSD+dih @
# 1.0 A -- so plotting it directly mixes 1.0 A dihedral-gated nodes into what
# should be the 0.75 A graph. Restrict to the focal RMSD rows.
# (.plot_networks_overview isolates 0.75 A via sptce_all_thresholds instead and
# is unaffected.)
sptce_graph_nodes <- function(clu)
  dplyr::filter(clu, method == name_sptce_rmsd(CFG$focal_threshold))

METHOD_COLOURS <- c(
  setNames("#e31a1c", .sptce_focal_rmsd),
  setNames("#b30000", .sptce_focal_dih),
  setNames("#b30000", .sptce_bubble_dih),
  setNames(c("#7a0177","#c51b8a","#fb6a4a")[seq_along(.sptce_extras)], .sptce_extras),
  "SPACE-TCR" = "#e31a1c",  "SPACE-TCR (RMSD)" = "#e31a1c",
  "SPACE-TCR (RMSD+dih<0.2rad)" = "#7a0177",
  "SPACE-TCR (RMSD+dih<0.3rad)" = "#c51b8a",
  "SPACE-TCR (RMSD+dih<0.5rad)" = "#b30000",
  "clusTCR"             = "darkgoldenrod1",
  "Hamming_Distance"    = "darkorange",  "Hamming Distance"    = "darkorange",
  "HD"                  = "darkorange",  "LD"                  = "cyan3",
  "Levenshtein_Distance"= "cyan3",       "Levenshtein Distance"= "cyan3",
  "TCRMatch"            = "darkorchid4", "iSMART"              = "grey50",
  "GIANA"               = "brown1",      "GLIPH2"              = "pink3",
  "DeepTCR"             = "chartreuse3", "TCRdist3"            = "black"
)

GROUP_PALETTE <- c(
  "published"                  = "#9c9c9c",
  "SPACE-TCR RMSD"             = "#e31a1c",
  "SPACE-TCR RMSD+dih (focal)" = "#b30000",
  "SPACE-TCR extra"            = "#fb6a4a"
)

PALETTES <- list(
  mode_colour = c(
    "SPACE-TCR RMSD"            = "#e31a1c", "SPACE-TCR RMSD+dih"        = "#b30000",
    "SPACE-TCR RMSD+dih<0.2rad" = "#7a0177", "SPACE-TCR RMSD+dih<0.3rad" = "#c51b8a",
    "SPACE-TCR RMSD+dih<0.5rad" = "#fb6a4a"),
  mode_shape = c(
    "SPACE-TCR RMSD" = 16, "SPACE-TCR RMSD+dih" = 15,
    "SPACE-TCR RMSD+dih<0.2rad" = 17, "SPACE-TCR RMSD+dih<0.3rad" = 15,
    "SPACE-TCR RMSD+dih<0.5rad" = 18),
  mode_line = c(
    "SPACE-TCR RMSD" = "solid", "SPACE-TCR RMSD+dih" = "dashed",
    "SPACE-TCR RMSD+dih<0.2rad" = "solid", "SPACE-TCR RMSD+dih<0.3rad" = "solid",
    "SPACE-TCR RMSD+dih<0.5rad" = "dashed"),
  group_fill = GROUP_PALETTE,
  ld_class = c(
    "LD<=1 (near-identical)"     = "#A50026",
    "LD 2-3 (similar motif)"     = "#FC8D59",
    "LD 4-6 (sequence-diverse)"  = "#4393C3",
    "LD >6 (sequence-unrelated)" = "#2166AC"),
  spec_class = c(
    "Monospecific (1 epitope)"          = "#2166AC",
    "Bispecific (2 epitopes)"           = "#FC8D59",
    "Polyspecific (3-5 epitopes)"       = "#D73027",
    "Highly polyspecific (>5 epitopes)" = "#67000D"),
  size_bin = c(
    "2"     = "#fde725", "3"     = "#90d743", "4-5"   = "#35b779",
    "6-10"  = "#21918c", "11-20" = "#3b528b", ">20"   = "#440154")
)

CONV_LEVELS <- names(PALETTES$ld_class)

# compute_metrics() has always used `purity >= 0.9`; only this label disagreed.
METRIC_LABELS <- c(
  retention          = "Retention",
  purity             = "Purity",
  perc_clusters_gt09 = "% clusters\npurity \u2265 0.9",
  perc_tcrs_gt09     = "% TCRs in\nhigh-purity clusters"
)

# The method is SPACE-TCR everywhere it is rendered.
sptce_label <- function(x) gsub("SPTCE", "SPACE-TCR", as.character(x), fixed = TRUE)

EPITOPE_PALETTE <- c(
  "GILGFVFTL" = "#E69F00", "NLVPMVATV" = "#D55E00",
  "GLCTLVAML" = "#CC79A7", "FLCMKALLL" = "#0072B2",
  "NYNYLYRLF" = "#56B4E9", "LLFGYPVYV" = "#009E73",
  "Other"     = "#1a1a1a"
)

EPITOPE_ORGANISM <- c(
  "GILGFVFTL" = "Influenza A", "NLVPMVATV" = "CMV",
  "GLCTLVAML" = "EBV",         "FLCMKALLL" = "HIV (Gag)",
  "NYNYLYRLF" = "HIV (Nef)",   "LLFGYPVYV" = "HTLV-1"
)

ORGANISM_PALETTE <- c(
  "Influenza A" = "#E69F00", "CMV"        = "#D55E00",
  "EBV"         = "#CC79A7", "HIV (Gag)"  = "#0072B2",
  "HIV (Nef)"   = "#56B4E9", "HTLV-1"     = "#009E73",
  "Other"       = "#9c9c9c", "Unlabelled" = "#dddddd"
)

# Hardcoded published values from Jouannet et al. (2025), kept ONLY as a
# fallback for when the raw Output_methods/ files or the curated database are
# unavailable. When the files are present, run_stage5() recomputes these
# in-script with compute_published_metrics_exact() so that the published methods
# and SPACE-TCR are scored with byte-for-byte the same equations.
PUBLISHED_METRICS_FALLBACK <- tribble(
  ~method,                ~retention, ~purity, ~perc_clusters_gt09, ~perc_tcrs_gt09,
  "clusTCR",              0.09, 0.99, 0.97, 0.98,
  "TCRMatch",             0.09, 0.82, 0.64, 0.63,
  "GLIPH2",               0.22, 0.80, 0.51, 0.48,
  "iSMART",               0.18, 0.75, 0.45, 0.50,
  "Hamming_Distance",     0.27, 0.69, 0.36, 0.34,
  "GIANA",                0.25, 0.65, 0.36, 0.35,
  "DeepTCR",              0.92, 0.63, 0.20, 0.28,
  "Levenshtein_Distance", 0.31, 0.57, 0.34, 0.20,
  "TCRdist3",             0.44, 0.42, 0.07, 0.05
)
PUBLISHED_METRICS <- PUBLISHED_METRICS_FALLBACK
# Paper-exact per-unit cluster membership (one row per pair for paired methods,
# per pooled sequence for single-chain). Set as a side effect by
# compute_published_metrics_exact(); NULL means "fall back to beta-only".
PUBLISHED_CLUSTERS_EXACT <- NULL

# =================================================================
# 4. THEMES
# =================================================================
theme_sptce <- theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    axis.title       = element_text(size = 16, face = "bold"),
    axis.text        = element_text(size = 14, colour = "grey15"),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 13),
    legend.key.size  = unit(1.0, "lines"),
    plot.title       = element_text(size = 17, face = "bold"),
    strip.text       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 13, colour = "grey35"),
    plot.caption     = element_text(size = 14, colour = "grey50"),
    plot.margin      = margin(8, 12, 8, 8)
  )
theme_set(theme_sptce)

# =================================================================
# 5. UTILITY FUNCTIONS
# =================================================================
to_lgl <- function(x) if (is.logical(x)) x else x == "True"
add_allele <- function(x) ifelse(str_detect(x, "\\*"), x, paste0(x, "*01"))

save_plot <- function(name, plot, w, h, dir = CFG$jouannet_out_dir, dpi = 300) {
  ggsave(file.path(dir, name), plot, width = w, height = h, dpi = dpi)
  message("  saved: ", name)
}


# Qualitative palette for colouring clusters in the LD network graphs. The
# first 24 colours are a Kelly/Polychrome-style set chosen so neighbouring
# clusters are easy to tell apart; beyond 24 it rotates HCL hue with
# alternating luminance.
distinct_cluster_palette <- function(n) {
  base <- c(
    "#E6194B","#3CB44B","#4363D8","#F58231","#911EB4","#42D4F4","#F032E6",
    "#BFEF45","#FABED4","#469990","#DCBEFF","#9A6324","#800000","#AAFFC3",
    "#808000","#FFD8B1","#000075","#A9A9A9","#E6BEFF","#1B9E77","#D95F02",
    "#666666","#E7298A","#66A61E"
  )
  if (n <= length(base)) return(base[seq_len(n)])
  extra <- grDevices::hcl(
    h = (seq_len(n - length(base)) * 137.5) %% 360,
    c = 78,
    l = ifelse(seq_len(n - length(base)) %% 2 == 0, 45, 70)
  )
  c(base, extra)
}

# =================================================================
# 6. CORE METRIC FUNCTIONS
# =================================================================
adjusted_rand_index <- function(a, b) {
  ct <- table(a, b); n <- sum(ct)
  if (n < 2) return(NA_real_)
  sc <- function(x) sum(choose(x, 2))
  idx <- sc(as.vector(ct)); ra <- sc(rowSums(ct)); rb <- sc(colSums(ct))
  exp <- ra * rb / choose(n, 2); mx <- (ra + rb) / 2
  if (mx == exp) 0 else (idx - exp) / (mx - exp)
}

normalised_mutual_info <- function(a, b) {
  ct <- table(a, b); n <- sum(ct)
  if (n < 2) return(NA_real_)
  pxy <- ct / n; px <- rowSums(pxy); py <- colSums(pxy)
  H   <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  mi  <- H(px) + H(py) - H(as.vector(pxy))
  if (H(px) == 0 || H(py) == 0) 0 else 2 * mi / (H(px) + H(py))
}


# ---- metric adapters ----------------------------------------------------
# These present sptce_metrics.R results under the column names the figure code
# below reads. No counting happens here: every number comes from cluster_table(),
# jm_metrics() or jm_sensitivity_detail(). Where a frame has one row per
# cluster membership, `.member` (the row number) is the unit, as in the
# reference; SPACE-TCR frames carry one row per receptor.
.as_members <- function(df) mutate(df, .member = row_number())

compute_metrics <- function(df, total_n, min_size = METRIC$min_cluster_size,
                            hi = METRIC$high_purity_cut) {
  m <- jm_metrics(.as_members(df), denom = total_n, unit = ".member",
                  min_size = min_size, hi = hi)
  tibble(n_clusters         = m$n_clusters,
         n_clustered        = m$n_units_clustered,
         retention          = m$retention,
         purity             = m$purity,
         perc_clusters_gt09 = m$pct_clusters_hi,
         perc_tcrs_gt09     = m$pct_members_hi)
}

# Eq. 4 sensitivity (jm_sensitivity). The unit of counting is the receptor:
# an explicit unit vector, else the paired CDR3s, else the id. `n_input` must
# come from the input table; without it the answer is NA.
.sens_unit <- function(df, unit = NULL) {
  unit %||% {
    if (all(c("cdr3a", "cdr3b") %in% names(df))) paste(df$cdr3a, df$cdr3b)
    else if ("receptor_id" %in% names(df))       as.character(df$receptor_id)
    else if ("id" %in% names(df))                as.character(df$id)
    else                                          as.character(seq_len(nrow(df)))
  }
}
compute_sensitivity <- function(df, epi, n_input = NULL, unit = NULL) {
  d <- df %>% mutate(.u = .sens_unit(df, unit)) %>% filter(!is.na(.u), nzchar(.u))
  if (!nrow(d) && !is.null(n_input) && !is.na(n_input) && n_input > 0) return(0)
  jm_sensitivity(d, epi, n_input, unit = ".u")
}

# Used by run_stage5 and by the fig32 held-out comparison.
#
# `input_counts` is a named vector or a two-column table giving, per epitope,
# the number of INPUT receptors carrying it. Without it every sensitivity is NA.
sens_for_df <- function(df, input_counts = CFG$.sens_input_counts) {
  n_of <- function(epi) {
    if (is.null(input_counts)) return(NA_real_)
    if (is.data.frame(input_counts)) {
      hit <- input_counts[[2]][input_counts[[1]] == epi]
      return(if (length(hit)) hit[1] else NA_real_)
    }
    v <- input_counts[[epi]]
    if (is.null(v)) NA_real_ else v
  }
  map_dfr(CFG$jouannet_epitopes, \(epi) tibble(
    epitope          = epi,
    sensitivity      = compute_sensitivity(df, epi, n_input = n_of(epi)),
    n_input          = n_of(epi),
    n_in_any_cluster = sum(df$epitope == epi, na.rm = TRUE)
  ))
}

# Fill CFG$.sens_input_counts from the labelled input table. Called by the
# stages that have the input table in hand, so sens_for_df() has the right
# denominator without every call site having to pass it.
.set_sens_input_counts <- function(input_df, epi_col = "epitope") {
  stopifnot(epi_col %in% names(input_df))
  tab <- table(input_df[[epi_col]])
  CFG$.sens_input_counts <<- as.list(stats::setNames(as.integer(tab), names(tab)))
  invisible(CFG$.sens_input_counts)
}

# =================================================================
# 7. SEQUENCE UTILITY FUNCTIONS
# =================================================================

ld_class <- function(ld) as.character(ld_class_of(ld, CONV_LEVELS))

# Per-cluster dominant epitope and purity (cluster_table over the labelled
# members) plus within-cluster CDR3 Levenshtein statistics.
summarise_clusters <- function(df, seq_col) {
  keys <- intersect(c("method", "threshold", "mode"), names(df))
  df %>%
    group_by(across(all_of(keys))) %>%
    group_modify(function(g, k) {
      ct <- cluster_table(.as_members(filter(g, !is.na(epitope))),
                          unit = ".member", min_size = 1L)
      g %>%
        group_by(cluster_id) %>%
        summarise(n_members = n(), ld = list(ld_stats(.data[[seq_col]])),
                  .groups = "drop") %>%
        left_join(transmute(ct, cluster_id, dom_epitope,
                            epitope_purity = purity), by = "cluster_id") %>%
        select(cluster_id, n_members, dom_epitope, epitope_purity, ld)
    }) %>%
    ungroup() %>%
    unnest(ld)
}

# =================================================================
# 8. DIHEDRAL FUNCTIONS
# =================================================================
circ_var <- function(rad) {
  rad <- rad[!is.na(rad)]
  if (length(rad) < 2) return(NA_real_)
  1 - abs(mean(complex(real = cos(rad), imaginary = sin(rad))))
}


# SPACE-TCR clusters at one operating point (cluster_at), with the column names
# the rest of this file uses: `tcr_id` for the receptor, labels joined by `id`.
build_sptce_clusters <- function(dist_files, labels, rmsd_cut,
                                 dih_cut = NULL, dih_map = NULL,
                                 prefix = "c", keep_unlabelled = FALSE) {
  gated <- !is.null(dih_cut) && !is.null(dih_map)
  out <- cluster_at(dist_files, rename(labels, receptor_id = id), rmsd_cut,
                    dih_cut = if (gated) dih_cut else NA_real_,
                    dih_map = if (gated) dih_map else NULL,
                    keep_unlabelled = keep_unlabelled, prefix = prefix)
  if (!nrow(out)) return(out)
  out %>% select(-method, -threshold, -dih) %>% rename(tcr_id = receptor_id)
}

# =================================================================
# 10. PAPER-EXACT PUBLISHED-METHOD METRICS (curated cohort, in-script)
# =================================================================
# Reproduces, in R, the exact Jouannet et al. calculation the hardcoded table
# was taken from. Verified to reproduce all four metrics for every file-based
# method to the second decimal (clusTCR 0.09/0.99/0.97/0.98, etc.).
#
# Three things make this match the paper, none of which the beta-only
# compute_metrics path did (which is why that returned clusTCR 0.94):
#   * paired methods (clusTCR, GLIPH2, DeepTCR, TCRdist3) are scored on the
#     CDR3b+CDR3a PAIR, single-chain methods on the individual CDR3;
#   * a sequence/pair binding several epitopes is expanded to one row per
#     epitope, so the dominant-epitope count is taken over all annotations
#     while `size` stays the original cluster size;
#   * the high-purity threshold is `purity >= 0.9`.
# Retention denominators are the curated counts (8395 single-chain sequences /
# 4779 pairs) and are overridable.
compute_published_metrics_exact <- function(db_path, method_dir,
                                            n_seq_denom = 8395L,
                                            n_pair_denom = 4779L,
                                            verbose = TRUE) {
  if (!file.exists(db_path)) {
    warning("database not found; using PUBLISHED_METRICS_FALLBACK")
    PUBLISHED_CLUSTERS_EXACT <<- NULL
    return(PUBLISHED_METRICS_FALLBACK)
  }
  if (verbose) cat("  Computing published metrics in-script (paper-exact)...\n")
  
  # ---- 1. curate database exactly as database_processing.R ----
  raw <- suppressWarnings(read_delim(db_path, delim = "\t",
                                     show_col_types = FALSE, guess_max = 200000))
  tenx <- "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf"
  d <- raw %>%
    filter(Verified_score == 2, Identification_score > 4.3, !is.na(Epitope),
           as.character(PubMed_ID) != tenx, Cell_subset == "CD8",
           !is.na(V_alpha), !is.na(J_alpha), !is.na(V_beta), !is.na(J_beta)) %>%
    filter(nchar(CDR3_alpha) >= 6, nchar(CDR3_alpha) <= 23,
           nchar(CDR3_beta)  >= 6, nchar(CDR3_beta)  <= 23) %>%
    mutate(pair = paste0(CDR3_beta, "_", CDR3_alpha)) %>%
    group_by(Epitope) %>% mutate(.nb = n_distinct(pair)) %>% ungroup() %>%
    filter(.nb != 1) %>%
    filter(!pair %in% c("CASSDSRGTEAFF_CASSDSRGTEAFF",
                        "CASSPITGTGAYGYTF_CASSSVNEQYF")) %>%
    distinct() %>%
    mutate(CDR3b_CDR3a = paste0(CDR3_beta, CDR3_alpha))
  
  # multi-epitope maps (one row per annotation)
  pair_epi <- d %>% distinct(CDR3b_CDR3a, Epitope)
  cdr3_epi <- bind_rows(d %>% transmute(CDR3 = CDR3_alpha, Epitope),
                        d %>% transmute(CDR3 = CDR3_beta,  Epitope)) %>% distinct()
  
  n_seq  <- n_seq_denom  %||% n_distinct(c(d$CDR3_alpha, d$CDR3_beta))
  n_pair <- n_pair_denom %||% n_distinct(d$CDR3b_CDR3a)
  valid_pairs <- unique(d$CDR3b_CDR3a)
  
  # ---- 2. scoring: jm_reference_metrics(), as purity_functions.R scores ----
  # `pct_tcrs` picks the paper's per-method denominator for % TCRs:
  #   "size"  -> sum(size)        (single-chain, clusTCR, DeepTCR)
  #   "units" -> distinct pairs   (GLIPH2, TCRdist3)
  metric_onechain <- function(out)
    jm_reference_metrics(out, cdr3_epi, "CDR3",
                         retention = n_distinct(out$CDR3) / n_seq, pct_tcrs = "size")
  metric_pairing <- function(out, retention_units, cdr3_perc = "size")
    jm_reference_metrics(out, pair_epi, "CDR3b_CDR3a",
                         retention = retention_units / n_pair,
                         pct_tcrs = if (cdr3_perc == "pairs") "units" else "size")

  # ---- paper-exact cluster membership for the figures (one row per unit) ----
  mem <- list()
  pair2beta <- setNames(d$CDR3_beta, d$CDR3b_CDR3a)
  pair_epi1 <- pair_epi %>% arrange(CDR3b_CDR3a, Epitope) %>%
    distinct(CDR3b_CDR3a, .keep_all = TRUE)
  cdr3_epi1 <- cdr3_epi %>% arrange(CDR3, Epitope) %>%
    distinct(CDR3, .keep_all = TRUE)
  add_mem_pair <- function(method, df) {
    big <- filter(df, size > 1)
    if (!nrow(big)) return(invisible())
    mem[[method]] <<- big %>%
      transmute(method = method, cluster_id = paste0(method, "_", clusters),
                cdr3b = CDR3b, CDR3b_CDR3a = CDR3b_CDR3a) %>%
      left_join(pair_epi1, by = "CDR3b_CDR3a") %>%
      transmute(method, cluster_id, cdr3b, epitope = Epitope)
  }
  add_mem_one <- function(method, out) {
    o <- out %>% group_by(clusters) %>% mutate(size = n()) %>% ungroup() %>%
      filter(size > 1)
    if (!nrow(o)) return(invisible())
    mem[[method]] <<- o %>%
      transmute(method = method, cluster_id = paste0(method, "_", clusters),
                cdr3b = CDR3) %>%
      left_join(cdr3_epi1, by = c("cdr3b" = "CDR3")) %>%
      transmute(method, cluster_id, cdr3b, epitope = Epitope)
  }
  
  rows <- list()
  rd <- function(p) file.path(method_dir, p)
  
  # ---- 3. parse each method exactly as output_processing.R ----
  if (file.exists(rd("iSMART_output.txt"))) {                 # single-chain
    ism <- read_delim(rd("iSMART_output.txt"), delim = "\t", show_col_types = FALSE)
    names(ism)[1:2] <- c("CDR3", "clusters")
    rows[["iSMART"]] <- c("iSMART", metric_onechain(ism[, c("CDR3","clusters")]))
    add_mem_one("iSMART", ism[, c("CDR3","clusters")])
  }
  if (file.exists(rd("GIANA_output.txt"))) {                   # skip ## header
    gl <- readLines(rd("GIANA_output.txt"), warn = FALSE)
    gl <- gl[!startsWith(gl, "##")]
    sp <- str_split_fixed(gl, "\t", 3)
    gi <- tibble(CDR3 = sp[,1], clusters = sp[,2]) %>% filter(CDR3 != "")
    rows[["GIANA"]] <- c("GIANA", metric_onechain(gi))
    add_mem_one("GIANA", gi)
  }
  if (file.exists(rd("TCRMatch_output.txt"))) {                # edges -> components
    tm <- suppressWarnings(read_delim(rd("TCRMatch_output.txt"), delim = " ",
                                      col_names = c("from","to","score"),
                                      show_col_types = FALSE))
    tm <- tm %>% filter(!is.na(score), score != 1)
    nodes <- unique(c(tm$from, tm$to))
    idx <- setNames(seq_along(nodes), nodes)
    comp <- connected_components(length(nodes), idx[tm$from], idx[tm$to])
    tmc <- tibble(CDR3 = nodes, clusters = as.character(comp))
    rows[["TCRMatch"]] <- c("TCRMatch", metric_onechain(tmc))
    add_mem_one("TCRMatch", tmc)
  }
  if (file.exists(rd("clusTCR_output.csv"))) {                 # paired
    cl <- read_csv(rd("clusTCR_output.csv"), show_col_types = FALSE)
    names(cl)[1:2] <- c("CDR3b_CDR3a", "clusters")
    cs <- read_csv(rd("clusTCR_output_summary.csv"), show_col_types = FALSE)
    names(cs)[1] <- "clusters"
    clm <- cl %>% inner_join(cs[, c("clusters","size")], by = "clusters") %>%
      filter(CDR3b_CDR3a %in% valid_pairs)
    ru <- n_distinct(clm$CDR3b_CDR3a)
    rows[["clusTCR"]] <- c("clusTCR", metric_pairing(clm, ru, "size"))
    clm$CDR3b <- pair2beta[clm$CDR3b_CDR3a]
    add_mem_pair("clusTCR", clm)
  }
  if (file.exists(rd("GLIPH2_output_cluster.csv"))) {          # paired
    g2 <- read_csv(rd("GLIPH2_output_cluster.csv"), show_col_types = FALSE)
    g2 <- g2[, 1:18] %>% na.omit() %>% filter(pattern != "single") %>%
      rename(CDR3b = TcRb, CDR3a = TcRa, clusters = index) %>%
      group_by(clusters) %>% mutate(size = n()) %>% ungroup() %>%
      mutate(CDR3b_CDR3a = paste0(CDR3b, CDR3a))
    ru <- nrow(distinct(g2, CDR3b, CDR3a))
    rows[["GLIPH2"]] <- c("GLIPH2",
                          metric_pairing(g2[, c("CDR3b_CDR3a","clusters","size")], ru, "pairs"))
    add_mem_pair("GLIPH2", g2)
  }
  if (file.exists(rd("DeepTCR_output.csv"))) {                 # paired
    dt <- read_csv(rd("DeepTCR_output.csv"), show_col_types = FALSE) %>%
      rename(CDR3b = Beta_Sequences, CDR3a = Alpha_Sequences, clusters = cluster) %>%
      mutate(CDR3b_CDR3a = paste0(CDR3b, CDR3a))
    ru <- nrow(distinct(dt, CDR3b, CDR3a))
    rows[["DeepTCR"]] <- c("DeepTCR",
                           metric_pairing(dt[, c("CDR3b_CDR3a","clusters","size")], ru, "size"))
    add_mem_pair("DeepTCR", dt)
  }
  if (file.exists(rd("TCRdist3_hclust_output.csv"))) {         # paired, ; delim
    td <- read_delim(rd("TCRdist3_hclust_output.csv"), delim = ";",
                     show_col_types = FALSE) %>%
      group_by(clusters) %>% mutate(size = n()) %>% ungroup() %>%
      filter(size != 1) %>%
      mutate(CDR3b_CDR3a = paste0(CDR3b, CDR3a))
    ru <- nrow(distinct(td, CDR3b, CDR3a))
    rows[["TCRdist3"]] <- c("TCRdist3",
                            metric_pairing(td[, c("CDR3b_CDR3a","clusters","size")], ru, "pairs"))
    add_mem_pair("TCRdist3", td)
  }
  
  # Hamming / Levenshtein: distance==1 graph on the curated unique single-chain
  # sequences. Mirrors output_processing.R, where graph_from_data_frame keeps
  # ONLY sequences in at least one distance-1 edge (isolated ones are dropped),
  # which is what gives e.g. Hamming retention 0.27 rather than ~1.
  uniq_seq <- tibble(CDR3 = unique(c(d$CDR3_alpha, d$CDR3_beta)))
  dist_clusters <- function(method) {
    seqs <- uniq_seq$CDR3
    dm <- stringdistmatrix(seqs, seqs, method = method)
    ei <- which(dm == 1 & upper.tri(dm), arr.ind = TRUE)
    if (!nrow(ei)) return(NULL)
    keep  <- sort(unique(c(ei[,1], ei[,2])))    # edge vertices only
    remap <- match(seq_along(seqs), keep)       # old index -> compact index
    comp  <- connected_components(length(keep), remap[ei[,1]], remap[ei[,2]])
    tibble(CDR3 = seqs[keep], clusters = as.character(comp))
  }
  hd <- dist_clusters("hamming")
  if (!is.null(hd)) {
    rows[["Hamming_Distance"]] <- c("Hamming_Distance", metric_onechain(hd))
    add_mem_one("Hamming_Distance", hd)
  }
  ld <- dist_clusters("lv")
  if (!is.null(ld)) {
    rows[["Levenshtein_Distance"]] <- c("Levenshtein_Distance", metric_onechain(ld))
    add_mem_one("Levenshtein_Distance", ld)
  }
  
  out <- bind_rows(lapply(rows, function(r) tibble(
    method             = r[[1]],
    retention          = round(as.numeric(r[[2]]), 4),
    purity             = round(as.numeric(r[[3]]), 4),
    perc_clusters_gt09 = round(as.numeric(r[[4]]), 4),
    perc_tcrs_gt09     = round(as.numeric(r[[5]]), 4))))
  
  ord <- c("clusTCR","TCRMatch","GLIPH2","iSMART","Hamming_Distance",
           "GIANA","DeepTCR","Levenshtein_Distance","TCRdist3")
  out <- out %>% mutate(method = factor(method, levels = ord)) %>%
    arrange(method) %>% mutate(method = as.character(method))
  PUBLISHED_CLUSTERS_EXACT <<- bind_rows(mem)
  if (verbose) { cat("  Published metrics (in-script):\n"); print(as.data.frame(out)) }
  out
}

# #################################################################
##  PART 1 - CORE REPLACEMENTS
# #################################################################

# =================================================================
# 11. OVERLAY LAYOUT                     [replaces .layout_with_isolates]
# =================================================================
# A plain force-directed layout crowds the informative nodes into a central
# ball: layout_with_fr() at default iterations freezes while still compact,
# normalising to unit radius fixes that in place, and isolates placed outside
# the core push it to under half the linear extent. This layout addresses all
# three.
#
# `spread` is a deliberate distortion and should be described as one: it
# reassigns each node's radius by rank so that density is uniform over the disc.
# It preserves angular order and therefore neighbourhood structure, but not
# distances. These panels are read for WHICH nodes are coloured, not for metric
# topology, so the trade is worth making; set spread = 0 to disable it.

.spread_layout <- function(l, spread = 0.6) {
  n <- nrow(l)
  if (n < 3 || spread <= 0) return(l)
  r <- sqrt(rowSums(l^2))
  mx <- max(r)
  if (!is.finite(mx) || mx <= 0) return(l)
  theta <- atan2(l[, 2], l[, 1])
  # equal-area radii: the k-th of n points sits at sqrt((k - 0.5)/n)
  r_eq  <- sqrt((rank(r, ties.method = "first") - 0.5) / n)
  r_new <- (1 - spread) * (r / mx) + spread * r_eq
  cbind(r_new * cos(theta), r_new * sin(theta))
}

.relax_min_sep <- function(xy, min_sep, iters = 400) {
  n <- nrow(xy)
  if (n < 2 || min_sep <= 0) return(xy)
  # Exactly coincident points have no separating direction and the guard below
  # (d > 1e-9) would leave them stuck on top of each other forever, so break
  # ties with a jitter far smaller than min_sep before relaxing.
  dup <- duplicated(round(xy, 9)) | duplicated(round(xy, 9), fromLast = TRUE)
  if (any(dup))
    xy[dup, ] <- xy[dup, ] + matrix(stats::runif(sum(dup) * 2, -1, 1),
                                    ncol = 2) * min_sep * 0.01
  for (i in seq_len(iters)) {
    dx <- outer(xy[, 1], xy[, 1], "-")
    dy <- outer(xy[, 2], xy[, 2], "-")
    d  <- sqrt(dx * dx + dy * dy)
    diag(d) <- Inf
    close <- d < min_sep
    if (!any(close)) break
    # normalise by the number of crowding neighbours, or nodes in the densest
    # region overshoot and the layout oscillates instead of settling
    push <- ifelse(close & d > 1e-9, (min_sep - d) / d * 0.5, 0)
    ncl  <- pmax(rowSums(close), 1)
    xy[, 1] <- xy[, 1] + rowSums(push * dx) / ncl
    xy[, 2] <- xy[, 2] + rowSums(push * dy) / ncl
  }
  xy
}

.layout_with_isolates <- function(n_nodes, ei, seed = 1,
                                  core_frac  = 0.78,   # share of radius for the connected core
                                  spread     = 0.60,   # 0 = raw FR, 1 = uniform density
                                  min_sep    = 0.045,  # min centre-to-centre, core radii
                                  relax_iter = 400,
                                  ring_pad   = 1.04,   # first isolate ring, x core radius
                                  ring_step  = 0.06) {
  xy  <- matrix(NA_real_, n_nodes, 2)
  deg <- integer(n_nodes)
  if (nrow(ei)) {
    tb <- table(c(ei[, 1], ei[, 2]))
    deg[as.integer(names(tb))] <- as.integer(tb)
  }
  con <- which(deg > 0); iso <- which(deg == 0)
  R <- core_frac
  set.seed(seed)
  
  if (length(con) >= 2 && requireNamespace("igraph", quietly = TRUE)) {
    remap <- match(seq_len(n_nodes), con)
    g <- igraph::graph_from_edgelist(cbind(remap[ei[, 1]], remap[ei[, 2]]),
                                     directed = FALSE)
    g <- igraph::simplify(g)
    if (igraph::vcount(g) < length(con))
      g <- igraph::add_vertices(g, length(con) - igraph::vcount(g))
    l <- tryCatch(
      igraph::layout_with_fr(g, niter = 3000,
                             start.temp = sqrt(igraph::vcount(g))),
      error = function(e) igraph::layout_with_fr(g))
    l <- sweep(l, 2, colMeans(l))
    mx <- max(sqrt(rowSums(l^2)))
    if (is.finite(mx) && mx > 0) l <- l / mx
    l <- .spread_layout(l, spread = spread)
    l <- .relax_min_sep(l, min_sep = min_sep, iters = relax_iter)
    mx <- max(sqrt(rowSums(l^2)))
    if (is.finite(mx) && mx > 0) l <- l / mx
    xy[con, ] <- l * core_frac
  } else if (length(con) == 1) {
    xy[con, ] <- c(0, 0); R <- 0.25
  } else if (length(con) >= 2) {
    a <- 2 * pi * (seq_along(con) - 1) / length(con)
    xy[con, 1] <- cos(a) * core_frac; xy[con, 2] <- sin(a) * core_frac
  } else R <- 0.25
  
  if (length(iso)) {
    # Isolates carry no structure, so pack them into as few tight rings as
    # possible immediately outside the core instead of letting them claim area.
    k <- length(iso)
    per_ring <- max(60L, ceiling(sqrt(k) * 6))
    ring <- ceiling(seq_along(iso) / per_ring)
    pos  <- ave(seq_along(iso), ring, FUN = seq_along)
    cnt  <- ave(seq_along(iso), ring, FUN = length)
    radius <- R * ring_pad + (ring - 1) * ring_step * R
    angle  <- 2 * pi * (pos - 1) / pmax(cnt, 1) + ring * 0.35
    xy[iso, 1] <- radius * cos(angle)
    xy[iso, 2] <- radius * sin(angle)
  }
  xy
}

# =================================================================
# 12. FOCAL-EPITOPE OVERLAYS            [replaces .plot_cluster_overlays]
# =================================================================
# Uses the spread layout above, a legend below the panel (on the right it consumed a third of the plot width),
# fewer named clusters, and smaller points now that the core is spread.
#
# .focal_universe() and .assign_overlay_factor() are defined in PART 2.
.plot_cluster_overlays <- function(s5, focal_epi = CFG$focal_epitope, ld_thr = 2,
                                   universe = "curated_db", max_nodes = 2500,
                                   seed = 1, n_named = 12L, pt_scale = 42) {
  cat("  fig39/39b/39c cluster overlays for", focal_epi,
      "(universe =", universe, ")...\n")
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("  overlay: igraph required; skipping."); return(invisible(NULL))
  }
  uni <- .focal_universe(s5, focal_epi, universe)
  if (is.null(uni) || !nrow(uni)) {
    message("  overlay: no TCRs for ", focal_epi, "; skipping."); return(invisible(NULL))
  }
  
  # A cdr3b can appear in several rows (paired-alpha expansion) and so in more
  # than one cluster. Rule: most frequent cluster, ties to the larger cluster
  # then to the lower id, so the choice is deterministic.
  .dominant_cluster <- function(df) {
    if (is.null(df) || !nrow(df)) return(tibble(cdr3b = character(), cid = character()))
    csize <- df %>% count(cluster_id, name = "cluster_size")
    df %>% count(cdr3b, cluster_id, name = "n_rows") %>%
      left_join(csize, by = "cluster_id") %>%
      arrange(cdr3b, desc(n_rows), desc(cluster_size), cluster_id) %>%
      distinct(cdr3b, .keep_all = TRUE) %>%
      transmute(cdr3b, cid = as.character(cluster_id))
  }
  sptce_map <- .dominant_cluster(sptce_graph_nodes(s5$sptce_clusters)) %>%
    rename(sptce_cid = cid)
  cc_map <- .dominant_cluster(filter(s5$published_clusters,
                                     method %in% c("clusTCR", "ClusTCR"))) %>%
    rename(clustcr_cid = cid)
  nodes <- uni %>%
    left_join(sptce_map, by = "cdr3b") %>%
    left_join(cc_map,    by = "cdr3b") %>%
    mutate(in_sptce = !is.na(sptce_cid), in_clustcr = !is.na(clustcr_cid))
  
  n_uni <- nrow(nodes)
  if (n_uni > max_nodes) {
    warning(sprintf("overlay: %d nodes exceeds max_nodes=%d; the LD matrix is O(n^2).",
                    n_uni, max_nodes), call. = FALSE)
    return(invisible(NULL))
  }
  if (n_uni < 2) { message("  overlay: <2 nodes; skipping."); return(invisible(NULL)) }
  
  dm <- stringdistmatrix(nodes$cdr3b, nodes$cdr3b, method = "lv")
  ei <- which(dm <= ld_thr & upper.tri(dm), arr.ind = TRUE)
  if (!nrow(ei)) ei <- matrix(integer(0), ncol = 2)
  xy <- .layout_with_isolates(n_uni, ei, seed = seed)
  nodes$x <- xy[, 1]; nodes$y <- xy[, 2]
  edge_df <- if (nrow(ei))
    tibble(x = xy[ei[, 1], 1], y = xy[ei[, 1], 2],
           xend = xy[ei[, 2], 1], yend = xy[ei[, 2], 2])
  else tibble(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0))
  
  pt_sz  <- max(1.0, min(3.0, pt_scale / sqrt(n_uni)))
  edge_a <- max(0.05, min(0.20, 180 / max(nrow(edge_df), 1)))
  n_sp   <- sum(nodes$in_sptce);  n_cc <- sum(nodes$in_clustcr)
  n_both <- sum(nodes$in_sptce & nodes$in_clustcr)
  n_none <- sum(!nodes$in_sptce & !nodes$in_clustcr)
  cat(sprintf("    universe %d | SPACE-TCR %d (%.1f%%) | clusTCR %d (%.1f%%) | both %d | neither %d\n",
              n_uni, n_sp, 100 * n_sp / n_uni, n_cc, 100 * n_cc / n_uni, n_both, n_none))
  
  short_sptce <- function(x) gsub("^r[0-9.]+_bin_", "bin",
                                  gsub("\\.pairwise\\.tsv\\.gz", "", as.character(x)))
  base_theme <- theme_void() +
    theme(legend.position   = "bottom",
          legend.title      = element_text(size = 12, face = "bold"),
          legend.text       = element_text(size = 9),
          legend.key.size   = unit(0.4, "cm"),
          plot.title        = element_text(size = 17, face = "bold"),
          plot.subtitle     = element_text(size = 11, colour = "grey40"),
          plot.margin       = margin(6, 6, 6, 6))
  
  panel <- function(fill_factor, pal, title, subtitle, legend_title) {
    # Greys are the last two factor levels, so painting in descending level
    # order puts them down first and the named clusters on top. Sorting by
    # label text would break: each panel uses a different "Other ..." label.
    d <- nodes %>% mutate(.ov = fill_factor) %>% arrange(desc(as.integer(.ov)))
    ggplot() +
      geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
                   colour = "grey65", linewidth = 0.2, alpha = edge_a) +
      geom_point(data = d, aes(x, y, colour = .ov), size = pt_sz, alpha = 0.95) +
      scale_colour_manual(values = pal, name = legend_title, drop = FALSE) +
      coord_fixed() + base_theme +
      guides(colour = guide_legend(nrow = 3, byrow = TRUE,
                                   override.aes = list(size = 3.2))) +
      labs(title = title, subtitle = subtitle)
  }
  sub_line <- function(n_in, who)
    sprintf("All %d %s TCRs in the %s. %s clustered %d (%.1f%%); grey = not clustered. Edges: CDR3b LD <= %d.",
            n_uni, focal_epi, universe, who, n_in, 100 * n_in / n_uni, ld_thr)
  
  cc_ov <- .assign_overlay_factor(nodes$clustcr_cid, n_named,
                                  label_other = "Other clusTCR cluster",
                                  label_none  = "Not clustered")
  save_plot(sprintf("fig39_overlay_clustcr_%s.png", focal_epi),
            panel(cc_ov$factor, cc_ov$palette,
                  sprintf("clusTCR clustering of all %s TCRs", focal_epi),
                  sub_line(n_cc, "clusTCR"), "clusTCR cluster"), 11, 10)
  
  sp_ov <- .assign_overlay_factor(nodes$sptce_cid, n_named,
                                  label_other = "Other SPACE-TCR cluster",
                                  label_none  = "Not clustered", shorten = short_sptce)
  save_plot(sprintf("fig39b_overlay_sptce_%s.png", focal_epi),
            panel(sp_ov$factor, sp_ov$palette,
                  sprintf("SPACE-TCR clustering of all %s TCRs", focal_epi),
                  sub_line(n_sp, "SPACE-TCR"), "SPACE-TCR cluster"), 11, 10)
  
  # fig39c strips cluster identity out, so the only thing compared is WHICH
  # receptors each method captured on an identical layout.
  bin_pal <- c("Clustered" = "#D73027", "Not clustered" = "#DADADA")
  bin_panel <- function(flag, who, n_in) {
    d <- nodes %>%
      mutate(.ov = factor(ifelse(flag, "Clustered", "Not clustered"),
                          levels = names(bin_pal))) %>%
      arrange(desc(.ov == "Not clustered"))
    ggplot() +
      geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
                   colour = "grey70", linewidth = 0.2, alpha = edge_a) +
      geom_point(data = d, aes(x, y, colour = .ov), size = pt_sz, alpha = 0.95) +
      scale_colour_manual(values = bin_pal, name = NULL, drop = FALSE) +
      coord_fixed() + base_theme +
      guides(colour = guide_legend(override.aes = list(size = 3.2))) +
      labs(title = who,
           subtitle = sprintf("%d / %d captured (%.1f%%)", n_in, n_uni, 100 * n_in / n_uni))
  }
  save_plot(sprintf("fig39c_capture_%s.png", focal_epi),
            bin_panel(nodes$in_sptce,
                      sprintf("SPACE-TCR (RMSD %s A)", CFG$focal_threshold), n_sp) +
              bin_panel(nodes$in_clustcr, "clusTCR", n_cc) +
              patchwork::plot_layout(guides = "collect") +
              patchwork::plot_annotation(
                title = sprintf("Which %s TCRs does each method capture?", focal_epi),
                subtitle = sprintf("Identical node set and layout: all %d %s TCRs in the %s. Both: %d. Neither: %d.",
                                   n_uni, focal_epi, universe, n_both, n_none),
                theme = theme(plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
                              plot.subtitle = element_text(size = 12, colour = "grey40", hjust = 0.5))) &
              theme(legend.position = "bottom"),
            16, 9)
  
  deg <- integer(n_uni)
  if (nrow(ei)) { tb <- table(c(ei[, 1], ei[, 2])); deg[as.integer(names(tb))] <- as.integer(tb) }
  memb <- nodes %>%
    transmute(cdr3b, epitope = focal_epi,
              in_curated_db = !is.na(src_db), in_sptce_input = !is.na(src_input),
              in_sptce, sptce_cluster = sptce_cid,
              in_clustcr, clustcr_cluster = clustcr_cid, ld_degree = deg,
              capture = case_when(in_sptce & in_clustcr ~ "both",
                                  in_sptce   ~ "SPACE-TCR only",
                                  in_clustcr ~ "clusTCR only",
                                  TRUE       ~ "neither"))
  write_csv(memb, file.path(CFG$jouannet_out_dir,
                            sprintf("overlay_membership_%s.csv", focal_epi)))
  invisible(list(nodes = memb, n_universe = n_uni, n_sptce = n_sp,
                 n_clustcr = n_cc, n_both = n_both, n_neither = n_none))
}

# =================================================================
# 13. DIHEDRAL-MODE PANEL
# =================================================================
sptce_dihedral_mode_panel <- function(all_sptce, sptce_dih_all_modes, total_n) {
  if (is.null(sptce_dih_all_modes) || !nrow(sptce_dih_all_modes)) {
    message("  dih panel: no dihedral modes; skipping."); return(invisible(NULL))
  }
  # The RMSD-only column is the reference the prune is judged against; without
  # it the panel cannot answer "what does the gate buy".
  combined <- bind_rows(
    all_sptce %>% filter(mode == "SPACE-TCR RMSD") %>%
      mutate(mode = "SPACE-TCR RMSD only"),
    sptce_dih_all_modes)
  
  mode_metrics <- combined %>% group_by(threshold, mode) %>%
    group_modify(~ compute_metrics(.x, total_n) %>%
                   mutate(method = paste0(.y$mode, " (", .y$threshold, "A)"))) %>%
    ungroup()
  write_csv(mode_metrics,
            file.path(CFG$jouannet_out_dir, "sptce_dihedral_mode_metrics.csv"))
  message("  wrote sptce_dihedral_mode_metrics.csv (", nrow(mode_metrics), " rows)")
  
  lv <- c("SPACE-TCR RMSD only",
          paste0("SPACE-TCR RMSD+dih<", sort(CFG$dihedral_modes), "rad"))
  p <- mode_metrics %>%
    pivot_longer(names(METRIC_LABELS), names_to = "metric", values_to = "value") %>%
    mutate(metric      = factor(metric, names(METRIC_LABELS), METRIC_LABELS),
           mode        = factor(mode, levels = lv),
           mode_short  = factor(sub("^SPACE-TCR ", "", as.character(mode)),
                                levels = sub("^SPACE-TCR ", "", lv)),
           threshold_f = factor(threshold)) %>%
    filter(!is.na(mode_short)) %>%
    ggplot(aes(mode_short, value, colour = threshold_f, group = threshold_f,
               shape = threshold_f)) +
    geom_line(linewidth = 0.85, alpha = 0.9) + geom_point(size = 3) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_colour_brewer(palette = "RdYlBu", name = "RMSD threshold (\u00c5)") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8), name = "RMSD threshold (\u00c5)") +
    scale_y_continuous(labels = percent) +
    labs(title = "Effect of the dihedral prune across all RMSD thresholds",
         subtitle = paste0("Each line is one RMSD threshold; x-axis is prune tightness. ",
                           "The 0.5 rad column is the RMSD+dih mode of Tables 1 and 2."),
         x = "SPACE-TCR mode", y = NULL) +
    theme(axis.text.x = element_text(size = 10, angle = 12, hjust = 1))
  save_plot("fig09_dih_modes_panel.png", p, 15, 5)
  invisible(mode_metrics)
}

# =================================================================
# 14. COHORT ACCOUNTING       [replaces run_stage5d_cohort_accounting]
# =================================================================
# Accounts per SPACE-TCR mode: pooling modes would count a receptor clustered at
# two operating points twice and blend their sensitivities. n_jouannet_paired
# and n_sptce_input are taken from separate sources so n_modelled_gap is real.
run_stage5d_cohort_accounting <- function(s5 = NULL) {
  cat("\n=== Stage 5d: cohort accounting (per mode) ===\n")
  if (is.null(s5)) s5 <- run_stage5()
  clu <- s5$sptce_clusters
  if (is.null(clu) || !nrow(clu)) { cat("  No clusters\n"); return(invisible(NULL)) }
  
  db_path <- file.path(CFG$jouannet_repo,
                       "Data/Database/database_pooled_human_2023_03_15.txt")
  curated <- if (file.exists(db_path)) {
    load_jouannet_curated(db_path, cell_type = CFG$jouannet_cell_type,
                          verbose = FALSE) %>% count(epitope, name = "n_curated_cdr3b")
  } else {
    cat("  [warn] curated DB not found; n_curated_cdr3b will be NA\n")
    tibble(epitope = character(), n_curated_cdr3b = integer())
  }
  input_counts <- read_tsv(CFG$jouannet_labelled, show_col_types = FALSE) %>%
    filter(!is.na(epitope)) %>% count(epitope, name = "n_sptce_input_tcrs")
  
  modes <- unique(as.character(clu$method))
  cat("  modes:", paste(modes, collapse = " | "), "\n")
  
  acc <- map_dfr(modes, function(m) {
    sub <- filter(clu, method == m)
    in_cluster <- sub %>% filter(!is.na(epitope)) %>% count(epitope, name = "n_in_any_cluster")
    # Sensitivity (Eq. 4, jm_sensitivity_detail) on the paired receptor, the
    # unit n_sptce_input_tcrs counts. |Te| comes from the INPUT table.
    sub_r <- sub %>%
      mutate(.rid = paste(cdr3a, cdr3b)) %>%
      filter(!is.na(.rid), nzchar(.rid))
    qual <- map_dfr(CFG$jouannet_epitopes, function(epi) {
      n_tot <- input_counts$n_sptce_input_tcrs[input_counts$epitope == epi]
      n_tot <- if (length(n_tot)) n_tot[1] else 0L
      if (!n_tot) return(tibble(epitope = epi, n_in_qualifying = 0L,
                                sensitivity = NA_real_))
      r <- jm_sensitivity_detail(sub_r, epi, n_tot, unit = ".rid")
      tibble(epitope = epi, n_in_qualifying = as.integer(r$n_recovered),
             sensitivity = r$sensitivity)
    })
    tibble(method = m, epitope = CFG$jouannet_epitopes) %>%
      left_join(curated,      by = "epitope") %>%
      left_join(input_counts, by = "epitope") %>%
      left_join(in_cluster,   by = "epitope") %>%
      left_join(qual,         by = "epitope") %>%
      mutate(across(c(n_sptce_input_tcrs, n_in_any_cluster, n_in_qualifying),
                    ~ replace_na(.x, 0L)),
             n_not_clustered = pmax(0L, n_sptce_input_tcrs - n_in_any_cluster),
             sensitivity     = round(sensitivity, 4))
  }) %>%
    select(method, epitope, n_curated_cdr3b, n_sptce_input_tcrs,
           n_in_any_cluster, n_not_clustered, n_in_qualifying, sensitivity)
  
  # A mode cannot cluster more TCRs of an epitope than it was supplied.
  bad <- filter(acc, n_in_any_cluster > n_sptce_input_tcrs)
  if (nrow(bad)) {
    warning("cohort accounting impossible for ", nrow(bad),
            " rows: more clustered than supplied. This is the signature of ",
            "pooling operating points.\n",
            paste(capture.output(print(as.data.frame(bad))), collapse = "\n"),
            call. = FALSE)
  } else cat("  self-check passed (clustered <= input for every row)\n")

  # Eq. 4 cannot recover more receptors than were clustered, so sensitivity
  # cannot exceed the clustered fraction.
  bad2 <- acc %>%
    filter(n_sptce_input_tcrs > 0, !is.na(sensitivity),
           sensitivity > n_in_any_cluster / n_sptce_input_tcrs + 1e-9)
  if (nrow(bad2)) {
    warning("sensitivity exceeds the clustered fraction for ", nrow(bad2),
            " rows, which is impossible under Eq. 4.\n",
            paste(capture.output(print(as.data.frame(bad2))), collapse = "\n"),
            call. = FALSE)
  } else cat("  self-check passed (sensitivity <= clustered fraction)\n")
  
  # n_curated_cdr3b counts unique CDR3b; n_sptce_input_tcrs counts paired TCRs,
  # several of which can share a CDR3b. They are different units, and their
  # difference is not a modelling gap.
  write_csv(acc, file.path(CFG$jouannet_out_dir, "cohort_accounting.csv"))
  cat("\n  Cohort accounting (units: curated = unique CDR3b, input = paired TCRs)\n\n")
  print(as.data.frame(acc), row.names = FALSE); cat("\n")
  
  funnel <- acc %>%
    select(method, epitope, n_sptce_input_tcrs, n_in_any_cluster, n_in_qualifying) %>%
    pivot_longer(-c(method, epitope), names_to = "stage", values_to = "n") %>%
    mutate(stage = factor(stage,
                          c("n_sptce_input_tcrs", "n_in_any_cluster", "n_in_qualifying"),
                          c("1. SPACE-TCR\ninput", "2. In any\ncluster",
                            "3. In qualifying\ncluster")),
           epitope = factor(epitope, CFG$jouannet_epitopes))
  p <- ggplot(funnel, aes(stage, n, group = method, colour = method)) +
    geom_line(linewidth = 0.8, alpha = 0.85) + geom_point(size = 2.6) +
    geom_text(aes(label = n), vjust = -0.9, size = 2.8, show.legend = FALSE) +
    scale_colour_manual(values = METHOD_COLOURS, drop = FALSE, name = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.25))) +
    facet_wrap(~ epitope, scales = "free_y", ncol = 3) +
    labs(title = "SPACE-TCR cohort funnel, per operating point",
         subtitle = "Stage 2 is the sensitivity denominator. Each line is one SPACE-TCR mode.",
         x = NULL, y = "Number of TCRs") +
    theme(axis.text.x = element_text(size = 10), legend.position = "bottom")
  save_plot("fig40_cohort_funnel.png", p, 14, 9)
  cat("Stage 5d complete.\n")
  invisible(acc)
}

# =================================================================
# 15. NETWORKS OVERVIEW CHECK
# =================================================================
# Panel titles that recompute their own n from whatever frame reaches
# layout_tcrs() can disagree with Table 2. Call this before writing any caption
# from that figure.
sptce_networks_overview_audit <- function(s5, min_cluster_size = 2L) {
  audit <- bind_rows(s5$sptce_all_thresholds %||% s5$sptce_clusters,
                     s5$published_clusters) %>%
    filter(!is.na(cdr3b)) %>%
    group_by(method, cluster_id) %>% mutate(csize = n()) %>%
    group_by(method) %>%
    summarise(rows_supplied        = n(),
              clusters_plotted     = n_distinct(cluster_id[csize >= min_cluster_size]),
              rows_plotted         = sum(csize >= min_cluster_size),
              unique_cdr3b_plotted = n_distinct(cdr3b[csize >= min_cluster_size]),
              .groups = "drop") %>%
    mutate(method = as.character(method))
  write_csv(audit, file.path(CFG$jouannet_out_dir, "networks_overview_audit.csv"))
  message("  wrote networks_overview_audit.csv — reconcile against Table 2 ",
          "before writing the caption")
  invisible(audit)
}

# =================================================================
# 16. LEAKAGE ON BOTH DENOMINATORS
# =================================================================
# Overlap against the full cohort, alongside the modelled-input figure.
sptce_leakage_full_cohort <- function() {
  parsed_csv <- CFG$stcrdab_parsed_csv
  db_path <- file.path(CFG$jouannet_repo,
                       "Data/Database/database_pooled_human_2023_03_15.txt")
  if (!file.exists(parsed_csv) || !file.exists(db_path)) {
    message("  leakage (full cohort): inputs missing; skipping."); return(invisible(NULL))
  }
  parsed <- read_csv(parsed_csv, show_col_types = FALSE) %>%
    filter(!is.na(cdr3), nchar(cdr3) >= 6)
  if ("v_gene" %in% names(parsed)) {
    b <- filter(parsed, grepl("TRBV", v_gene)); if (nrow(b)) parsed <- b
  }
  stcrdab_cdr3 <- unique(parsed$cdr3)
  curated <- load_jouannet_curated(db_path, cell_type = CFG$jouannet_cell_type,
                                   verbose = FALSE)
  cohort_cdr3 <- unique(curated$cdr3b)
  out <- tibble(set = c("Curated CD8 cohort", "SPACE-TCR modelled input"),
                n_unique_cdr3b = c(length(cohort_cdr3), NA_integer_),
                n_in_stcrdab   = c(length(intersect(cohort_cdr3, stcrdab_cdr3)), NA_integer_),
                pct            = NA_real_)
  inp <- tryCatch(read_tsv(CFG$jouannet_labelled, show_col_types = FALSE),
                  error = function(e) NULL)
  if (!is.null(inp)) {
    ic <- unique(inp$cdr3_aa_beta); ic <- ic[!is.na(ic)]
    out$n_unique_cdr3b[2] <- length(ic)
    out$n_in_stcrdab[2]   <- length(intersect(ic, stcrdab_cdr3))
  }
  out$pct <- 100 * out$n_in_stcrdab / pmax(out$n_unique_cdr3b, 1)
  dir.create(.stcrdab_dir(), showWarnings = FALSE, recursive = TRUE)
  write_csv(out, file.path(.stcrdab_dir(), "stcrdab_overlap_both_denominators.csv"))
  cat("\n  STCRDab overlap on both denominators:\n\n")
  print(as.data.frame(out), row.names = FALSE); cat("\n")
  invisible(out)
}

# =================================================================
# BENCHMARK A FIGURES
# =================================================================
# fig_paired_thresholds() draws the thresholds panel with the gates in
# BENCH_A_FIG_DIHEDRAL; fig_paired_methods() draws the methods panel RMSD-only
# (dih_keep = NA), as in Figure 9.
#
# Why these gates. score_paired_all() computes every gate in
# CFG$dihedral_modes, but 5 thresholds x 4 modes would be 20 lines. At each
# RMSD threshold:
#   0.5 rad  changes purity under full dilution by <= 0.013 and recovery by
#            <= 1.3 pp, so it lies close to the RMSD-only line.
#   0.2 rad  reduces coverage to ~1% of annotated sequences (0.4-1.1%), so all
#            five thresholds converge.
#   0.3 rad  separates visibly while keeping usable coverage (1.2-4.3%).
#
# At a fixed RMSD the gate raises purity under noise (0.75 A at 89.2% noise:
# 0.304 -> 0.367 with a 0.3 rad gate). Matched on coverage instead,
# tightening the RMSD threshold gives higher purity than adding the gate in 9
# of 11 comparable pairs on this cohort and 9 of 9 on the curated cohort; the
# two exceptions are 0.5 rad at the same RMSD, ahead by ~0.01. For example,
# 0.75 A + dih<0.3 gives 2.0% recovery at 0.367 purity, while 0.5 A alone gives
# 1.8% recovery at 0.624.
fig_paired_methods <- function(res, rmsd = c(0.5, 0.75, 1.0),
                               dih_keep = NA,
                               min_clusters = 5L, size_metric = "mean") {
  sc <- .size_col(res, size_metric)
  d0 <- filter(res, !is.na(purity))
  sp <- d0 %>% filter(str_starts(method, "SPACE-TCR"), param %in% rmsd,
                      is.na(dih) | dih %in% dih_keep)
  pb <- filter(d0, !str_starts(method, "SPACE-TCR"))
  if (!nrow(sp)) {
    message("No SPACE-TCR rows for RMSD ", paste(rmsd, collapse = ", "),
            " / dih ", paste(dih_keep, collapse = ", "))
    return(invisible(NULL))
  }
  
  cov <- sp %>% group_by(param, dih) %>%
    summarise(clusters = sprintf("%d-%d", min(n_clusters), max(n_clusters)),
              recovery_pct = sprintf("%.1f-%.1f", 100 * min(signal_recovery),
                                     100 * max(signal_recovery)), .groups = "drop")
  cat("\n  SPACE-TCR coverage by operating point (a high line on very few\n")
  cat("  clusters is not a strong operating point):\n\n")
  print(as.data.frame(cov), row.names = FALSE); cat("\n")
  
  # The ladder is arithmetic in counts but saturating in percent: the rungs sit
  # at 0, 62, 77, 83, 87 and 89% unannotated, so a linear % axis would spend
  # ~70% of its width on the first step. The rungs are plotted as an ordered
  # factor (evenly spaced), and each tick shows the noise fraction and the size
  # of the repertoire at that rung.
  lv  <- sort(unique(res$pct_noise))
  # score_paired_all() names the per-rung repertoire size `n_universe`; the
  # Python scorer in ml/ called it `n_input`. Accept either, and fall back to a
  # percentage-only label rather than inventing a count.
  ncol_ <- intersect(c("n_universe", "n_input", "n_total"), names(res))[1]
  n_at <- if (!is.na(ncol_)) {
    res %>% group_by(pct_noise) %>%
      summarise(n = max(.data[[ncol_]], na.rm = TRUE), .groups = "drop") %>%
      arrange(pct_noise) %>% pull(n)
  } else rep(NA_integer_, length(lv))
  # Two lines per tick: the noise fraction, then the repertoire size. The "n = "
  # prefix is NOT repeated on every tick -- at six rungs the labels then overlap
  # and the axis becomes less readable than the bare integers it replaced. The
  # prefix lives in the axis title instead.
  rung_lab <- ifelse(is.na(n_at),
                     sprintf("%.0f%%", lv),
                     sprintf("%.0f%%\n%s", lv, format(n_at, big.mark = ",",
                                                      trim = TRUE)))
  rung <- function(x) factor(sprintf("%.0f", x), levels = sprintf("%.0f", lv),
                             labels = rung_lab)
  prep <- function(x) x %>%
    mutate(rung = rung(pct_noise), thin = n_clusters < min_clusters) %>%
    select(any_of(c("method", "param", "dih")), rung, purity,
           size = all_of(sc$col), thin) %>%
    pivot_longer(c(purity, size), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(recode(metric, purity = "Cluster purity", size = sc$lab),
                           c("Cluster purity", sc$lab)))
  
  # Published methods are repeated in every RMSD facet as context, so each panel
  # compares one SPACE-TCR operating point against the whole field.
  thr_lab <- function(p) factor(sprintf("SPACE-TCR RMSD %s \u00c5", p),
                                levels = sprintf("SPACE-TCR RMSD %s \u00c5", sort(rmsd)))
  # prep() keeps param/dih via any_of(), and for published rows both are NA
  # columns that exist. Crossing in a new `param` would then collide with the
  # existing one ("Names must be unique"), so drop them first.
  dp <- prep(pb) %>% select(-any_of(c("param", "dih"))) %>%
    tidyr::crossing(param = rmsd) %>% mutate(facet = thr_lab(param))
  ds <- prep(sp) %>% mutate(facet = thr_lab(param), gated = !is.na(dih),
                            grp = paste0(param, "_", coalesce(as.character(dih), "none")))
  
  p <- ggplot(mapping = aes(rung, value)) +
    geom_line(data = dp, aes(group = method, colour = method), linewidth = 0.6, alpha = 0.85) +
    geom_line(data = ds, aes(group = grp, linetype = gated), colour = "#B2182B", linewidth = 1.25) +
    geom_point(data = ds, aes(group = grp, shape = thin), colour = "#B2182B", size = 2.1) +
    facet_grid(metric ~ facet, scales = "free_y", switch = "y") +
    scale_colour_manual(values = METHOD_COLOURS, drop = FALSE, name = "Published method",
                        guide = guide_legend(order = 2, nrow = 3)) +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"),
                          labels = c(`FALSE` = "SPACE-TCR (RMSD only)",
                                     `TRUE` = sprintf("SPACE-TCR (+ dihedral < %s rad)",
                                                      paste(dih_keep, collapse = ", "))),
                          name = NULL, guide = guide_legend(order = 1, nrow = 2)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                       labels = c(`FALSE` = "5+ clusters",
                                  `TRUE` = sprintf("< %d clusters", min_clusters)),
                       name = NULL, guide = guide_legend(order = 3, nrow = 2)) +
    labs(title = "Purity and cluster size as repertoire noise increases",
         subtitle = paste0("A fixed core of epitope-labelled paired TCRs is diluted into a ",
                           "growing background of unannotated ones, left to right.\nOnly the ",
                           "labelled core can contribute to purity, so a method that absorbs ",
                           "the background is penalised. Same sequences,\nlabels, purity ",
                           "definition and minimum cluster size for every method; published ",
                           "methods repeated in each panel as context.\nSteps are evenly ",
                           "spaced because the ladder is arithmetic in counts but saturates ",
                           "in percent."),
         x = "Unannotated share of the repertoire (second line: n paired TCRs)",
         y = NULL) +
    theme(plot.subtitle = element_text(size = 10.5, colour = "grey35"),
          panel.grid.minor = element_blank(), legend.position = "bottom",
          legend.box = "horizontal", strip.placement = "outside",
          axis.text.x = element_text(size = 9.5, lineheight = 1.05))
  ggsave(file.path(CFG$analysis_dir, "fig_paired_methods.png"), p,
         width = 13.5, height = 8, dpi = 200)
  message("  saved: fig_paired_methods.png")
  invisible(bind_rows(dp, select(ds, -any_of(c("gated", "grp")))))
}

# Figure 9's thresholds panel draws both gates reported in the paper, 0.3 and
# 0.5 rad, so it can be read against Tables 1 and 2 (which use 0.5 for the
# RMSD+dih rows). Pass dih_keep to draw a different set.
BENCH_A_FIG_DIHEDRAL <- c(0.3, 0.5)

fig_paired_thresholds <- function(res, min_clusters = 5L, size_metric = "mean",
                                  dih_keep = BENCH_A_FIG_DIHEDRAL) {
  sc <- .size_col(res, size_metric)
  d <- res %>% filter(str_starts(method, "SPACE-TCR"), !is.na(purity))
  # The argument is named dih_keep, not dih, because inside filter() a bare
  # `dih` refers to the column.
  if (!all(is.na(dih_keep))) d <- filter(d, is.na(dih) | dih %in% dih_keep)
  else                       d <- filter(d, is.na(dih))
  d <- d %>%
    mutate(thin  = n_clusters < min_clusters,
           thr   = factor(param, levels = sort(unique(param))),
           # A factor over the gate, so each gate gets its own line type.
           gated = factor(coalesce(as.character(dih), "none"),
                          levels = c("none", sort(unique(na.omit(as.character(dih)))))),
           grp   = paste0(param, "_", coalesce(as.character(dih), "none")),
           tag   = .sptce_tag(param, dih)) %>%
    select(param, dih, thr, gated, grp, tag, pct_noise, purity,
           size = all_of(sc$col), thin, n_clusters)
  if (!nrow(d)) { message("no SPACE-TCR rows"); return(invisible(NULL)) }
  
  thin_thr <- d %>% filter(thin) %>% count(tag, name = "n_sparse_levels")
  if (nrow(thin_thr)) {
    cat("  operating points with sparse levels (kept, drawn hollow):\n")
    print(as.data.frame(thin_thr), row.names = FALSE)
  }
  
  lt_levels <- levels(d$gated)
  lt_values <- setNames(c("solid", "22", "41", "1343")[seq_along(lt_levels)], lt_levels)
  # Short labels: with two gates drawn the long form ran off the device and the
  # last legend entry was clipped. The subtitle carries the full wording.
  lt_labels <- setNames(ifelse(lt_levels == "none", "RMSD only",
                               sprintf("+ gate %s rad", lt_levels)),
                        lt_levels)
  gate_lab <- paste(setdiff(lt_levels, "none"), collapse = " and ")
  dl <- d %>%
    pivot_longer(c(purity, size), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(recode(metric, purity = "Cluster purity", size = sc$lab),
                           c("Cluster purity", sc$lab)))
  
  p <- ggplot(dl, aes(pct_noise, value, colour = thr, group = grp)) +
    geom_line(aes(linetype = gated), linewidth = 1) +
    geom_point(aes(shape = thin), size = 2.4) +
    scale_linetype_manual(values = lt_values, labels = lt_labels, name = NULL,
                          guide = guide_legend(nrow = 2, order = 2,
                                               keywidth = unit(2.1, "lines"))) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                       labels = c(`FALSE` = "5+ clusters",
                                  `TRUE` = sprintf("< %d clusters", min_clusters)),
                       name = NULL) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_colour_brewer(palette = "RdYlBu", name = "RMSD (\u00c5)") +
    # This panel keeps a linear axis, unlike fig_paired_methods(), to show all
    # thresholds converging once the background dominates.
    scale_x_continuous(breaks = seq(0, 90, 20),
                       labels = function(b) sprintf("%g%%", b)) +
    labs(title = "SPACE-TCR: purity and cluster size across operating points",
         # The subtitle names the gates drawn and summarises the others.
         subtitle = paste0("Same paired-only cohort, diluted into a growing unannotated ",
                           "background. Broken lines are the ", gate_lab,
                           " rad dihedral gates;\nthe 0.2 rad gate is not drawn because it ",
                           "collapses every threshold to about 1% coverage. Hollow points ",
                           "have fewer\nthan ", min_clusters, " clusters. Linear axis, so ",
                           "the convergence at the noisy end is shown to scale."),
         x = "Unannotated share of the repertoire", y = NULL) +
    theme(plot.subtitle = element_text(size = 11, colour = "grey35"),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  # 13.5in, not 12: the legend gained the gate label and the last entry was
  # being clipped at the right edge of the device.
  ggsave(file.path(CFG$analysis_dir, "fig_paired_thresholds.png"), p,
         width = 14, height = 7.4, dpi = 200)
  message("  saved: fig_paired_thresholds.png")
  invisible(dl)
}

# run_benchmark_a() draws the field-comparison panel RMSD-only (dih_fig = NA),
# as in the paper: "we therefore use RMSD-only at 0.75 A for the primary
# comparison and reserve dihedral gating for looser thresholds". To draw a
# gated panel:
#     run_benchmark_a(dih_fig = CFG$focal_dihedral)
# Every gate in CFG$dihedral_modes is written to paired_only_comparison.csv
# whichever is drawn.

# =================================================================
# PART 2 — EXTENSIONS
#
# Two different things are called "consensus" in this codebase and they are
# unrelated. Gene consensus (section 10) is how uniform ONE cluster is in its
# V or J gene usage. Consensus clustering (section 12) is the intersection of
# TWO methods' cluster assignments. Both are standard terms in their own
# literature.
# =================================================================

.ext_dir <- function() {
  d <- file.path(CFG$jouannet_out_dir, "extensions")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}
.ext_save <- function(name, plot, w = 10, h = 7, dpi = 200) {
  path <- file.path(.ext_dir(), name)
  ggsave(path, plot, width = w, height = h, dpi = dpi)
  message("  wrote ", path); invisible(path)
}
.ext_write_csv <- function(x, name) {
  path <- file.path(.ext_dir(), name)
  write_csv(x, path); message("  wrote ", path); invisible(path)
}

# =================================================================
# 17. LABEL CONSENSUS                           [replaces compute_purity_generic]
# =================================================================
# Fraction of a cluster accounted for by its most common value of `label_col`.
# When label_col is the epitope this is purity in the usual sense, because the
# epitope is a ground-truth annotation. When it is a V or J gene call there is
# no ground truth, so the same number is a consensus: a description of how
# uniform the cluster is, not of whether it is correct.
compute_label_consensus <- function(df, label_col = "epitope", min_size = 2,
                                    hi_thresh = 0.9) {
  stopifnot(label_col %in% names(df), "cluster_id" %in% names(df))
  cs <- df %>%
    filter(!is.na(.data[[label_col]]), .data[[label_col]] != "") %>%
    group_by(cluster_id) %>%
    summarise(size     = n(),
              n_dom    = { t <- table(.data[[label_col]]); if (length(t)) max(t) else 0L },
              n_labels = n_distinct(.data[[label_col]]),
              .groups  = "drop") %>%
    filter(size >= min_size) %>%
    mutate(consensus = n_dom / size)
  
  is_truth <- label_col == "epitope"
  if (!nrow(cs))
    return(tibble(label = label_col, measure = if (is_truth) "purity" else "consensus",
                  n_clusters = 0L, n_clustered = 0L, global_consensus = NA_real_,
                  perc_clusters_hi = NA_real_, perc_tcrs_hi = NA_real_,
                  mean_labels_per_cluster = NA_real_))
  
  tibble(label                   = label_col,
         measure                 = if (is_truth) "purity" else "consensus",
         n_clusters              = nrow(cs),
         n_clustered             = sum(cs$size),
         global_consensus        = sum(cs$n_dom) / sum(cs$size),
         perc_clusters_hi        = mean(cs$consensus >= hi_thresh),
         perc_tcrs_hi            = sum(cs$size[cs$consensus >= hi_thresh]) / sum(cs$size),
         mean_labels_per_cluster = mean(cs$n_labels))
}

# =================================================================
# 18. V- AND J-GENE CONSENSUS                    [replaces sptce_gene_purity]
# =================================================================
.strip_allele <- function(x) if (is.null(x)) NULL else sub("\\*.*$", "", x)

.load_gene_labels <- function() {
  for (f in c(CFG$jouannet_labelled, CFG$master_file)) {
    if (!length(f) || !file.exists(f)) next
    lab <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(lab)) next
    nm  <- tolower(names(lab))
    ren <- list(
      v_call_beta   = c("v_call_beta","vb","v_gene_beta","tcrbv","trbv"),
      j_call_beta   = c("j_call_beta","jb","j_gene_beta","tcrbj","trbj"),
      v_call_alpha  = c("v_call_alpha","va","v_gene_alpha","tcrav","trav"),
      j_call_alpha  = c("j_call_alpha","ja","j_gene_alpha","tcraj","traj"),
      cdr3_aa_beta  = c("cdr3_aa_beta","cdr3b","cdr3.beta","beta_cdr3"),
      cdr3_aa_alpha = c("cdr3_aa_alpha","cdr3a","cdr3.alpha","alpha_cdr3"))
    have <- list()
    for (tgt in names(ren)) {
      hit <- which(nm %in% ren[[tgt]])[1]
      if (!is.na(hit)) have[[tgt]] <- names(lab)[hit]
    }
    if (is.null(have$cdr3_aa_beta) ||
        (is.null(have$v_call_beta) && is.null(have$j_call_beta))) next
    out <- lab %>%
      transmute(cdr3b        = .data[[have$cdr3_aa_beta]],
                v_gene_beta  = if (!is.null(have$v_call_beta))  .strip_allele(.data[[have$v_call_beta]])  else NA_character_,
                j_gene_beta  = if (!is.null(have$j_call_beta))  .strip_allele(.data[[have$j_call_beta]])  else NA_character_,
                v_gene_alpha = if (!is.null(have$v_call_alpha)) .strip_allele(.data[[have$v_call_alpha]]) else NA_character_,
                j_gene_alpha = if (!is.null(have$j_call_alpha)) .strip_allele(.data[[have$j_call_alpha]]) else NA_character_) %>%
      filter(!is.na(cdr3b), nchar(cdr3b) >= 6) %>%
      group_by(cdr3b) %>%
      summarise(across(everything(),
                       ~ { t <- table(.[!is.na(.)]); if (length(t)) names(t)[which.max(t)] else NA_character_ }),
                .groups = "drop")
    message("  gene labels: loaded from ", basename(f), " (", nrow(out),
            " unique CDR3b, ", sum(!is.na(out$v_gene_beta)), " with V, ",
            sum(!is.na(out$j_gene_beta)), " with J)")
    return(out)
  }
  message("  gene labels: no V/J columns found; gene consensus skipped")
  NULL
}

sptce_gene_consensus <- function(s5, min_size = 2, hi_thresh = 0.9) {
  gl <- .load_gene_labels()
  if (is.null(gl)) return(invisible(NULL))
  
  combined <- bind_rows(
    s5$sptce_clusters     %>% mutate(method = as.character(method)),
    s5$published_clusters %>% mutate(method = as.character(method))) %>%
    left_join(gl, by = "cdr3b") %>%
    select(method, cluster_id, cdr3b, epitope, v_gene_beta, j_gene_beta)
  
  tbl <- combined %>% group_by(method) %>%
    group_modify(~ bind_rows(
      compute_label_consensus(.x, "epitope",     min_size, hi_thresh),
      compute_label_consensus(.x, "v_gene_beta", min_size, hi_thresh),
      compute_label_consensus(.x, "j_gene_beta", min_size, hi_thresh))) %>%
    ungroup() %>%
    mutate(label = recode(label, epitope = "Epitope (purity)",
                          v_gene_beta = "V gene (consensus)",
                          j_gene_beta = "J gene (consensus)"),
           label = factor(label, c("Epitope (purity)", "V gene (consensus)",
                                   "J gene (consensus)")))
  .ext_write_csv(tbl, "gene_consensus_by_method.csv")
  
  # Gene calls are keyed to CDR3b, so methods that pool CDR3a and CDR3b are
  # scored for gene consensus only over clusters holding >= 2 CDR3b. That
  # denominator is reported alongside.
  denom <- tbl %>% select(method, label, n_clusters) %>%
    pivot_wider(names_from = label, values_from = n_clusters) %>%
    rename(clusters_epitope = `Epitope (purity)`,
           clusters_gene    = `V gene (consensus)`) %>%
    mutate(gene_denominator_frac = round(clusters_gene / clusters_epitope, 3)) %>%
    select(method, clusters_epitope, clusters_gene, gene_denominator_frac)
  .ext_write_csv(denom, "gene_consensus_denominators.csv")
  cat("\n  Gene-consensus denominators (clusters scorable for V/J vs for epitope):\n\n")
  print(as.data.frame(denom), row.names = FALSE); cat("\n")
  
  ord <- tbl %>% filter(label == "Epitope (purity)") %>%
    arrange(desc(replace_na(perc_clusters_hi, -Inf))) %>% pull(method)
  ord <- c(ord, setdiff(unique(tbl$method), ord))
  # Fix the label order (V gene, then J gene) rather than ggplot's alphabetical
  # default.
  tbl <- mutate(tbl, method = factor(sptce_label(method), levels = sptce_label(ord)),
                label  = factor(label, levels = c("Epitope (purity)",
                                                  "V gene (consensus)",
                                                  "J gene (consensus)")))
  
  p <- ggplot(tbl, aes(method, perc_clusters_hi, fill = label)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = percent(perc_clusters_hi, accuracy = 1)),
              position = position_dodge(width = 0.75), vjust = -0.3,
              size = 2.6, na.rm = TRUE) +
    scale_y_continuous(labels = percent, limits = c(0, 1.05),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c("Epitope (purity)"    = "#2C7FB8",
                                 "V gene (consensus)"  = "#7FCDBB",
                                 "J gene (consensus)"  = "#EDF8B1")) +
    labs(title    = "V- and J-gene consensus in clusters",
         subtitle = paste0("% of clusters in which one label accounts for \u2265 ",
                           round(100 * hi_thresh), "% of members. For the epitope ",
                           "label this is purity; for gene calls there is no ground ",
                           "truth, so it measures uniformity only."),
         x = NULL, y = sprintf("%% clusters with \u2265 %d%% consensus", round(100 * hi_thresh)),
         fill = "Label") +
    theme_bw(base_size = 11) +
    # The method names are long and rotated. A rotated label is anchored at its
    # right-hand end, so it overhangs on the left; the left margin is widened
    # to keep the first label on the page.
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.major.x = element_blank(),
          plot.subtitle = element_text(size = 8.5, colour = "grey35"),
          plot.margin = margin(8, 12, 8, 46))
  .ext_save("fig_gene_consensus_by_method.png", p, w = 12, h = 6.6)
  invisible(list(consensus = tbl, denominators = denom))
}

# =================================================================
# 19. PER-POSITION CDR3 ENTROPY
# =================================================================
.shannon_bits <- function(v) {
  v <- v[!is.na(v) & v != ""]
  if (!length(v)) return(NA_real_)
  p <- prop.table(table(v)); -sum(p * log2(p))
}

.cluster_position_entropy <- function(seqs, use_len = NULL) {
  if (!length(seqs)) return(tibble(pos = integer(), entropy = numeric()))
  L <- if (is.null(use_len))
    as.integer(names(sort(table(nchar(seqs)), decreasing = TRUE))[1]) else use_len
  keep <- seqs[nchar(seqs) == L]
  if (length(keep) < 2) return(tibble(pos = seq_len(L), entropy = NA_real_))
  m <- do.call(rbind, strsplit(keep, "", fixed = TRUE))
  tibble(pos = seq_len(L), entropy = apply(m, 2, .shannon_bits),
         n_seq = nrow(m), length = L)
}

sptce_position_entropy_heatmap <- function(s5, min_size = 4, min_purity = 0.9,
                                           top_n = 30, chain = c("beta", "alpha")) {
  chain <- match.arg(chain)
  seq_col <- if (chain == "beta") "cdr3b" else "cdr3a"
  df <- s5$sptce_clusters
  if (!seq_col %in% names(df)) {
    message("  position entropy: no ", seq_col, " column; skipping."); return(invisible(NULL))
  }
  summ <- df %>%
    filter(!is.na(.data[[seq_col]]), nchar(.data[[seq_col]]) >= 6) %>%
    group_by(method, cluster_id) %>%
    summarise(size   = n(),
              dom    = { t <- table(epitope[!is.na(epitope)]); if (length(t)) names(t)[which.max(t)] else NA_character_ },
              purity = { t <- table(epitope[!is.na(epitope)]); if (length(t)) max(t)/sum(t) else NA_real_ },
              seqs   = list(.data[[seq_col]]), .groups = "drop") %>%
    filter(size >= min_size, !is.na(purity), purity >= min_purity)
  if (!nrow(summ)) {
    message("  position entropy: no clusters meet size >= ", min_size,
            " & purity >= ", min_purity, "; skipping."); return(invisible(NULL))
  }
  # Rank-prefix the row label: SPACE-TCR cluster ids share long common prefixes,
  # so truncating alone would collide.
  summ_top <- summ %>% arrange(desc(purity), desc(size)) %>% head(top_n) %>%
    mutate(rank = row_number(),
           row_id = sprintf("%02d. %s | %s | n=%d | pur=%.0f%%", rank,
                            substr(cluster_id, 1, 14), dom, size, 100 * purity))
  ent_long <- pmap_dfr(list(summ_top$row_id, summ_top$seqs), function(row_id, seqs) {
    e <- .cluster_position_entropy(seqs)
    if (!nrow(e)) return(NULL)
    mutate(e, row_id = row_id)
  })
  if (!nrow(ent_long)) { message("  position entropy: nothing computed."); return(invisible(NULL)) }
  ent_long <- ent_long %>%
    mutate(row_id = factor(row_id, levels = rev(summ_top$row_id)),
           pos    = factor(pos, levels = sort(unique(pos))))
  p1 <- ggplot(ent_long, aes(pos, row_id, fill = entropy)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradient(low = "#FFFFCC", high = "#B10026", name = "Shannon\n(bits)",
                        na.value = "grey85",
                        limits = c(0, max(ent_long$entropy, na.rm = TRUE))) +
    labs(title = sprintf("CDR3%s per-position entropy in high-purity clusters",
                         if (chain == "beta") "b" else "a"),
         subtitle = sprintf("Top %d clusters at size >= %d, epitope purity >= %.1f",
                            nrow(summ_top), min_size, min_purity),
         x = "CDR3 position (N -> C)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(family = "mono", size = 7),
          panel.grid = element_blank())
  .ext_save(sprintf("fig_position_entropy_clusters_%s.png", chain), p1,
            w = 10, h = max(4, 0.28 * nrow(summ_top) + 2))
  
  by_len <- pmap_dfr(list(summ$cluster_id, summ$seqs, summ$purity, summ$size),
                     function(cid, seqs, pur, sz) {
                       e <- .cluster_position_entropy(seqs)
                       if (!nrow(e)) return(NULL)
                       mutate(e, cluster_id = cid, purity = pur, cluster_size = sz)
                     })
  ok_lens <- by_len %>% filter(!is.na(entropy)) %>% group_by(length) %>%
    summarise(n_clu = n_distinct(cluster_id), .groups = "drop") %>%
    filter(n_clu >= 2) %>% pull(length)          # avoid one-cluster stripes
  by_len_summ <- by_len %>% filter(!is.na(entropy), length %in% ok_lens) %>%
    group_by(length, pos) %>%
    summarise(mean_entropy = weighted.mean(entropy, cluster_size, na.rm = TRUE),
              n_clusters = n_distinct(cluster_id), .groups = "drop")
  if (nrow(by_len_summ)) {
    p2 <- ggplot(by_len_summ, aes(pos, factor(length), fill = mean_entropy)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      geom_text(aes(label = sprintf("%.2f", mean_entropy)), size = 2.5, colour = "grey20") +
      scale_fill_gradient(low = "#FFFFCC", high = "#B10026", name = "Mean\nbits") +
      scale_x_continuous(breaks = 1:22) +
      labs(title = sprintf("CDR3%s positional variability by length",
                           if (chain == "beta") "b" else "a"),
           subtitle = sprintf("Weighted mean Shannon entropy across %d high-purity clusters",
                              n_distinct(by_len$cluster_id)),
           x = "CDR3 position (N -> C)", y = "CDR3 length (aa)") +
      theme_bw(base_size = 10) + theme(panel.grid = element_blank())
    .ext_save(sprintf("fig_position_entropy_by_length_%s.png", chain), p2, w = 11, h = 5)
    .ext_write_csv(by_len_summ, sprintf("position_entropy_by_length_%s.csv", chain))
  }
  invisible(list(cluster_level = ent_long, by_length = by_len_summ))
}

# =================================================================
# 20. CONSENSUS CLUSTERING (method x method)
# =================================================================
# COUNTING CONVENTION: built and scored on UNIQUE CDR3b. For a strict
# intersection of two clusterings to be comparable, each unit must be counted at
# most once per method, or the Cartesian join between paired-alpha expansions
# and single-chain output inflates the intersection above either component -
# impossible for a set intersection. The denominator is n_distinct(cdr3b), NOT
# s5$total_n (paired-TCR rows), so these numbers differ from the
# stage-5 table.
#
# One consensus per (sptce_mode x baseline); pooling the SPACE-TCR focal modes
# would double-count every cdr3b clustering in both. Formally, two receptors
# share a consensus cluster iff both methods place them together.
.dedupe_method <- function(df, epi_map = NULL) {
  base <- df %>% distinct(cdr3b, cluster_id)
  if (is.null(epi_map)) {
    epi_map <- df %>% filter(!is.na(epitope)) %>% group_by(cdr3b) %>%
      summarise(epitope = { t <- table(epitope); names(sort(t, decreasing = TRUE))[1] },
                .groups = "drop")
  }
  left_join(base, epi_map, by = "cdr3b")
}

# Metrics on unique CDR3b (jm_metrics with unit = cdr3b).
.metrics_unique <- function(df, universe_n, hi = METRIC$high_purity_cut,
                            min_size = METRIC$min_cluster_size) {
  if (!nrow(df))
    return(tibble(n_clusters = 0L, n_clustered = 0L, n_cdr3b_cluster_pairs = 0L,
                  retention = 0, purity = NA_real_,
                  perc_clusters_gt09 = NA_real_, perc_tcrs_gt09 = NA_real_))
  m <- jm_metrics(df, universe_n, unit = "cdr3b", min_size = min_size, hi = hi)
  tibble(n_clusters            = m$n_clusters,
         n_clustered           = m$n_units_clustered,
         n_cdr3b_cluster_pairs = m$n_memberships,
         retention             = m$retention,
         purity                = m$purity,
         perc_clusters_gt09    = m$pct_clusters_hi,
         perc_tcrs_gt09        = m$pct_members_hi)
}

sptce_consensus_clusters <- function(s5, sptce_modes = NULL,
                                     baseline_methods = c("clusTCR", "GLIPH2"),
                                     min_size = 2) {
  if (is.null(sptce_modes)) sptce_modes <- unique(as.character(s5$sptce_clusters$method))
  message("  consensus: sptce modes = ", paste(sptce_modes, collapse = ", "))
  global_epi <- bind_rows(select(s5$sptce_clusters, cdr3b, epitope),
                          select(s5$published_clusters, cdr3b, epitope)) %>%
    filter(!is.na(epitope)) %>% group_by(cdr3b) %>%
    summarise(epitope = { t <- table(epitope); names(sort(t, decreasing = TRUE))[1] },
              .groups = "drop")
  universe <- unique(c(s5$sptce_clusters$cdr3b, s5$published_clusters$cdr3b))
  universe <- universe[!is.na(universe) & nchar(universe) >= 6]
  universe_n <- length(universe)
  message("  consensus: universe of unique cdr3b = ", universe_n)
  
  sptce_by_mode <- map(sptce_modes, ~ s5$sptce_clusters %>% filter(method == .x) %>%
                         .dedupe_method(epi_map = global_epi)) %>% setNames(sptce_modes)
  pub_by_method <- map(baseline_methods, ~ s5$published_clusters %>% filter(method == .x) %>%
                         .dedupe_method(epi_map = global_epi)) %>% setNames(baseline_methods)
  pub_by_method <- compact(pub_by_method[map_lgl(pub_by_method, ~ nrow(.x) > 0)])
  
  .intersect <- function(sp, other, sp_name, other_name) {
    ints <- inner_join(transmute(sp, cdr3b, sptce_cid = cluster_id),
                       transmute(other, cdr3b, other_cid = cluster_id),
                       by = "cdr3b", relationship = "many-to-many") %>% distinct()
    if (!nrow(ints)) {
      message("  consensus (", sp_name, " + ", other_name, "): no shared cdr3b.")
      return(NULL)
    }
    ints %>%
      mutate(cluster_id = sprintf("SPT_%s__%s_%s", sptce_cid, other_name, other_cid)) %>%
      group_by(cluster_id) %>% filter(n_distinct(cdr3b) >= min_size) %>% ungroup() %>%
      left_join(global_epi, by = "cdr3b") %>%
      transmute(method = paste0(sp_name, " + ", other_name), cluster_id, cdr3b, epitope)
  }
  cons_all <- map_dfr(names(sptce_by_mode), function(sm)
    map_dfr(names(pub_by_method), function(om)
      .intersect(sptce_by_mode[[sm]], pub_by_method[[om]], sm, om)))
  if (!nrow(cons_all)) {
    message("  consensus: none at min_size=", min_size, "; skipping."); return(invisible(NULL))
  }
  
  sptce_metrics <- map_dfr(names(sptce_by_mode), ~ .metrics_unique(sptce_by_mode[[.x]], universe_n) %>%
                             mutate(method = .x))
  pub_metrics <- map_dfr(names(pub_by_method), ~ .metrics_unique(pub_by_method[[.x]], universe_n) %>%
                           mutate(method = .x))
  cons_metrics <- cons_all %>% group_by(method) %>%
    group_modify(~ .metrics_unique(.x, universe_n)) %>% ungroup()
  metrics_all <- bind_rows(sptce_metrics, pub_metrics, cons_metrics) %>%
    mutate(universe_n = universe_n, is_consensus = grepl(" \\+ ", method),
           unit = "unique cdr3b")
  
  # SELF-CHECK: an intersection cannot retain more than either component.
  audit <- cons_metrics %>%
    mutate(parts      = strsplit(as.character(method), " \\+ "),
           sptce_part = vapply(parts, function(x) x[1], character(1)),
           other_part = vapply(parts, function(x) x[2], character(1)),
           sptce_ret  = sptce_metrics$retention[match(sptce_part, sptce_metrics$method)],
           other_ret  = pub_metrics$retention[match(other_part, pub_metrics$method)],
           cons_ret   = retention,
           ret_ok     = (cons_ret <= sptce_ret + 1e-9) & (cons_ret <= other_ret + 1e-9)) %>%
    select(method, sptce_ret, other_ret, cons_ret, ret_ok)
  if (!all(audit$ret_ok, na.rm = TRUE))
    warning("Consensus retention exceeds a component's - something is miscounting:\n",
            paste(capture.output(print(as.data.frame(audit))), collapse = "\n"))
  else message("  consensus retention check passed for all ", nrow(audit), " consensuses.")
  
  .ext_write_csv(metrics_all, "consensus_metrics.csv")
  .ext_write_csv(cons_all,    "consensus_clusters_membership.csv")
  .ext_write_csv(audit,       "consensus_retention_selfcheck.csv")
  
  lab <- if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(aes(label = method), size = 3.2, box.padding = 0.4,
                             min.segment.length = 0, max.overlaps = 100, show.legend = FALSE)
  else geom_text(aes(label = method), vjust = -1.2, size = 3.2, show.legend = FALSE)
  
  p <- ggplot(metrics_all, aes(retention, purity, colour = is_consensus, shape = is_consensus)) +
    geom_point(size = 4) + lab +
    scale_x_continuous(labels = percent, expand = expansion(mult = 0.08)) +
    scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
    scale_colour_manual(values = c(`FALSE` = "#3182BD", `TRUE` = "#D95F0E"),
                        labels = c(`FALSE` = "Base method", `TRUE` = "Consensus"), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    labs(title = "Consensus methods: purity vs retention",
         subtitle = sprintf("Unit = unique CDR3b | universe n = %d | intersection retention <= min(component) by construction",
                            universe_n),
         x = "Retention (fraction of universe)", y = "Global cluster purity") +
    theme_bw(base_size = 11)
  .ext_save("fig_consensus_purity_retention.png", p, w = 11, h = 6.5)
  
  p2 <- ggplot(metrics_all, aes(reorder(method, perc_clusters_gt09),
                                perc_clusters_gt09, fill = is_consensus)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = percent(perc_clusters_gt09, accuracy = 1)), hjust = -0.15, size = 3) +
    scale_fill_manual(values = c(`FALSE` = "#7BC1E7", `TRUE` = "#D95F0E"),
                      labels = c(`FALSE` = "Base method", `TRUE` = "Consensus"), name = NULL) +
    scale_y_continuous(labels = percent, limits = c(0, 1.1),
                       expand = expansion(mult = c(0, 0.03))) +
    coord_flip() +
    labs(title = "% clusters at purity >= 0.9",
         subtitle = "Baselines and consensuses on an identical unique-CDR3b denominator",
         x = NULL, y = "% clusters purity >= 0.9") +
    theme_bw(base_size = 11) + theme(panel.grid.major.y = element_blank())
  .ext_save("fig_consensus_perc_pure_clusters.png", p2, w = 10, h = 6.5)
  invisible(list(members = cons_all, metrics = metrics_all, audit = audit))
}

# =================================================================
# 21. COUNTING FUNNEL
# =================================================================
# Shows raw rows against de-duplicated unique cdr3b at every stage, so it is
# visible where multiplicity enters and that the intersection never exceeds
# either component.
sptce_counting_funnel <- function(s5, sptce_modes = NULL,
                                  baseline_methods = c("clusTCR", "GLIPH2"),
                                  min_size = 2) {
  if (is.null(sptce_modes)) sptce_modes <- unique(as.character(s5$sptce_clusters$method))
  all_cdr3b <- unique(c(s5$sptce_clusters$cdr3b, s5$published_clusters$cdr3b))
  all_cdr3b <- all_cdr3b[!is.na(all_cdr3b) & nchar(all_cdr3b) >= 6]
  U <- length(all_cdr3b)
  mk <- function(stage, method, n_rows, uniq)
    tibble(stage_id = stage, method = method, n_rows = as.integer(n_rows),
           n_unique_cdr3b = as.integer(length(uniq)),
           pct_of_universe = 100 * length(uniq) / U)
  rows <- list()
  rows[[1]] <- mk("0. Shared CDR3b universe", "- (all methods)", U, all_cdr3b)
  rows[[2]] <- tibble(stage_id = "1. Input rows (paired-TCR count, s5$total_n)",
                      method = "-", n_rows = as.integer(s5$total_n),
                      n_unique_cdr3b = NA_integer_, pct_of_universe = NA_real_)
  add_stages <- function(df, label) {
    rows[[length(rows) + 1]] <<- mk(sprintf("2. %s - raw rows", label), label,
                                    nrow(df), unique(df$cdr3b))
    dedup <- df %>% distinct(cdr3b, cluster_id)
    rows[[length(rows) + 1]] <<- mk(sprintf("3. %s - unique cdr3b in any cluster", label),
                                    label, nrow(dedup), unique(dedup$cdr3b))
    big <- dedup %>% group_by(cluster_id) %>% filter(n_distinct(cdr3b) >= min_size) %>% ungroup()
    rows[[length(rows) + 1]] <<- mk(sprintf("4. %s - unique cdr3b in cluster size >= %d", label, min_size),
                                    label, nrow(big), unique(big$cdr3b))
  }
  for (sm in sptce_modes) add_stages(filter(s5$sptce_clusters, method == sm), sm)
  for (om in baseline_methods) {
    df <- filter(s5$published_clusters, method == om)
    if (nrow(df)) add_stages(df, om)
  }
  for (sm in sptce_modes) {
    S <- s5$sptce_clusters %>% filter(method == sm) %>% distinct(cdr3b, cluster_id)
    for (om in baseline_methods) {
      O <- s5$published_clusters %>% filter(method == om) %>% distinct(cdr3b, cluster_id)
      if (!nrow(O)) next
      int_cdr3b <- intersect(unique(S$cdr3b), unique(O$cdr3b))
      rows[[length(rows) + 1]] <- mk(sprintf("5. %s AND %s - cdr3b in both", sm, om),
                                     sprintf("%s AND %s", sm, om),
                                     length(int_cdr3b), int_cdr3b)
      ints <- inner_join(transmute(S, cdr3b, sptce_cid = cluster_id),
                         transmute(O, cdr3b, other_cid = cluster_id),
                         by = "cdr3b", relationship = "many-to-many") %>% distinct() %>%
        mutate(cluster_id = sprintf("%s__%s_%s", sptce_cid, om, other_cid))
      big <- ints %>% group_by(cluster_id) %>% filter(n_distinct(cdr3b) >= min_size) %>% ungroup()
      rows[[length(rows) + 1]] <- mk(
        sprintf("6. %s + %s - cdr3b in consensus size >= %d", sm, om, min_size),
        sprintf("%s + %s", sm, om), nrow(big), unique(big$cdr3b))
    }
  }
  funnel <- bind_rows(rows)
  .ext_write_csv(funnel, "counting_funnel.csv")
  
  # Three telescoping counts per method. Gap raw->pairs is paired-alpha or
  # polyspecific row expansion; gap pairs->unique is overlapping clusters.
  method_names <- c(sptce_modes, baseline_methods)
  three <- map_dfr(method_names, function(m) {
    src <- if (m %in% sptce_modes) s5$sptce_clusters else s5$published_clusters
    df <- filter(src, method == m)
    if (!nrow(df)) return(NULL)
    lv <- c("Raw rows in frame", "Unique (cdr3b x cluster) pairs", "Unique cdr3b")
    tibble(method = m, kind = factor(lv, levels = lv),
           n = c(nrow(df), nrow(distinct(df, cdr3b, cluster_id)), n_distinct(df$cdr3b)))
  }) %>% mutate(method = factor(method, levels = method_names))
  p <- ggplot(three, aes(method, n, fill = kind)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(aes(label = formatC(n, big.mark = ",", format = "d")),
              position = position_dodge(width = 0.8), hjust = -0.1, size = 2.8) +
    scale_fill_manual(values = c("Raw rows in frame" = "#FDBB84",
                                 "Unique (cdr3b x cluster) pairs" = "#9ECAE1",
                                 "Unique cdr3b" = "#08519C"), name = NULL) +
    coord_flip() +
    labs(title = "How many CDR3b really? Three counting depths per method",
         subtitle = "Gap raw->pairs = paired-a / polyspecific row expansion; gap pairs->unique = overlapping clusters",
         x = NULL, y = "Count") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
  .ext_save("fig_row_vs_unique_multiplicity.png", p, w = 12, h = 6)
  message("  funnel: universe = ", U, " unique cdr3b")
  invisible(funnel)
}

# =================================================================
# 22. ALL-EPITOPE RANGES
# =================================================================
# The main tables report the 6 focal epitopes; this shows the distribution
# across every epitope each method clustered, to show how representative the
# six are.
sptce_all_epitope_ranges <- function(s5, min_size = 2, focal = NULL) {
  if (is.null(focal)) focal <- CFG$jouannet_epitopes
  per_clu <- bind_rows(
    transmute(s5$sptce_clusters,     method = as.character(method), cluster_id, cdr3b, epitope),
    transmute(s5$published_clusters, method = as.character(method), cluster_id, cdr3b, epitope)
  ) %>% filter(!is.na(epitope)) %>%
    group_by(method, cluster_id) %>%
    summarise(size = n(), dom = { t <- table(epitope); names(t)[which.max(t)] },
              n_dom = { t <- table(epitope); max(t) }, purity = n_dom / size,
              .groups = "drop") %>%
    filter(size >= min_size)
  
  all_epis <- unique(bind_rows(s5$sptce_clusters, s5$published_clusters)$epitope)
  all_epis <- all_epis[!is.na(all_epis)]
  sens_all <- bind_rows(mutate(s5$sptce_clusters,     method = as.character(method)),
                        mutate(s5$published_clusters, method = as.character(method))) %>%
    group_by(method) %>%
    group_modify(~ map_dfr(all_epis, function(epi) tibble(
      epitope = epi, sensitivity = compute_sensitivity(.x, epi),
      n_tcrs = sum(.x$epitope == epi, na.rm = TRUE))) %>% filter(n_tcrs >= 2)) %>%
    ungroup()
  .ext_write_csv(sens_all, "sensitivity_all_epitopes.csv")
  .ext_write_csv(per_clu,  "per_cluster_purity_size.csv")
  
  ord <- sens_all %>% group_by(method) %>%
    summarise(m = median(sensitivity, na.rm = TRUE), .groups = "drop") %>%
    arrange(replace_na(m, -Inf)) %>% pull(method)
  sens_all   <- mutate(sens_all, method = factor(method, levels = ord))
  focal_hits <- filter(sens_all, epitope %in% focal)
  p1 <- ggplot(sens_all, aes(method, sensitivity)) +
    geom_boxplot(fill = "#DEEBF7", outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.35, size = 1, colour = "grey40") +
    geom_point(data = focal_hits, aes(colour = epitope), size = 2.5,
               position = position_jitter(width = 0.15, seed = 42)) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(title = "Per-epitope sensitivity across ALL epitopes",
         subtitle = sprintf("Grey = every epitope with >= 2 TCRs (n=%d); coloured = the 6 focal epitopes",
                            n_distinct(sens_all$epitope)),
         x = NULL, y = "Sensitivity", colour = "Focal epitope") +
    coord_flip() + theme_bw(base_size = 11) + theme(panel.grid.major.y = element_blank())
  .ext_save("fig_sensitivity_all_epitopes.png", p1, w = 11, h = 7)
  
  ordp <- per_clu %>% group_by(method) %>%
    summarise(m = median(purity, na.rm = TRUE), .groups = "drop") %>%
    arrange(replace_na(m, -Inf)) %>% pull(method)
  p2 <- per_clu %>% mutate(method = factor(method, levels = ordp)) %>%
    ggplot(aes(method, purity)) +
    geom_violin(fill = "#C6DBEF", colour = "grey40", scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_y_continuous(labels = percent, limits = c(0, 1.02)) +
    labs(title = "Per-cluster purity distribution",
         subtitle = sprintf("All clusters with size >= %d, every method on one axis", min_size),
         x = NULL, y = "Cluster purity (dominant-epitope fraction)") +
    coord_flip() + theme_bw(base_size = 11) + theme(panel.grid.major.y = element_blank())
  .ext_save("fig_per_cluster_purity_dist.png", p2, w = 10, h = 6.5)
  
  rank_df <- sens_all %>% group_by(method) %>%
    arrange(desc(sensitivity), .by_group = TRUE) %>%
    mutate(rank = row_number(), is_focal = epitope %in% focal) %>% ungroup()
  p3 <- ggplot(rank_df, aes(rank, sensitivity)) +
    geom_line(colour = "grey55", linewidth = 0.4) +
    geom_point(aes(colour = is_focal, size = is_focal, alpha = is_focal)) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#E31A1C"),
                        labels = c(`FALSE` = "other", `TRUE` = "focal 6")) +
    scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2.2), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.5, `TRUE` = 1), guide = "none") +
    scale_y_continuous(labels = percent) +
    facet_wrap(~ method, scales = "free_x", ncol = 3) +
    labs(title = "Sensitivity rank-order per method",
         subtitle = "Red shows where the 6 focal epitopes fall in each method's full distribution",
         x = "Epitope rank (best first)", y = "Sensitivity", colour = NULL) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey92"), legend.position = "top")
  .ext_save("fig_sensitivity_rank_per_method.png", p3, w = 12, h = 8)
  invisible(list(sensitivity = sens_all, per_cluster = per_clu))
}

# =================================================================
# 23. WITHIN-CLUSTER CDR3a LD
# =================================================================
# Companion to fig05. Only methods given paired TRA+TRB input are eligible.
# The alpha is read from each METHOD'S OWN output, not via a cdr3b -> modal
# cdr3a lookup: those clusters are beta-similar but not necessarily
# alpha-similar, so a modal alpha is close to random where a cdr3b appears in
# several paired TCRs, and it inflates the alpha LD.
.load_deeptcr_paired <- function(dir) {
  f <- file.path(dir, "DeepTCR_output.csv"); if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  if (!all(c("Alpha_Sequences", "Beta_Sequences", "cluster") %in% names(d))) return(NULL)
  d %>% transmute(cluster_id = as.character(cluster),
                  cdr3b = as.character(Beta_Sequences),
                  cdr3a = as.character(Alpha_Sequences)) %>%
    filter(!is.na(cdr3a), !is.na(cdr3b), nchar(cdr3a) >= 6, nchar(cdr3b) >= 6)
}
.load_gliph2_paired <- function(dir) {
  f <- file.path(dir, "GLIPH2_output_cluster.csv"); if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  if (!all(c("TcRa", "TcRb", "index") %in% names(d))) return(NULL)
  d %>% transmute(cluster_id = as.character(index), cdr3b = as.character(TcRb),
                  cdr3a = as.character(TcRa)) %>%
    filter(!is.na(cdr3a), !is.na(cdr3b), nchar(cdr3a) >= 6, nchar(cdr3b) >= 6)
}
.load_tcrdist3_paired <- function(dir) {
  f <- file.path(dir, "TCRdist3_hclust_output.csv"); if (!file.exists(f)) return(NULL)
  d <- read_delim(f, delim = ";", show_col_types = FALSE)     # semicolon, not comma
  if (!all(c("CDR3a", "CDR3b", "clusters") %in% names(d))) return(NULL)
  d %>% transmute(cluster_id = as.character(clusters), cdr3b = as.character(CDR3b),
                  cdr3a = as.character(CDR3a)) %>%
    filter(!is.na(cdr3a), !is.na(cdr3b), nchar(cdr3a) >= 6, nchar(cdr3b) >= 6)
}
# Heuristic split of a concatenated BETA||ALPHA junction, used only as a
# fallback: prefer a split whose left half ends F/W and starts CASS/CAS/CSA and
# whose right half starts C, at typical CDR3 lengths.
.split_beta_alpha <- function(s) {
  if (is.na(s) || nchar(s) < 15) return(list(NA_character_, NA_character_))
  n <- nchar(s); best <- NULL; best_score <- -Inf
  for (bk in 8:(n - 8)) {
    b <- substr(s, 1, bk); a <- substr(s, bk + 1, n)
    if (substr(b, nchar(b), nchar(b)) %in% c("F", "W") && substr(a, 1, 1) == "C") {
      sc <- if (bk >= 12 && bk <= 16) 3 else 1
      an <- nchar(a); sc <- sc + if (an >= 10 && an <= 16) 3 else 1
      if (grepl("^CAS[SF]|^CSA", b)) sc <- sc + 2
      if (sc > best_score) { best_score <- sc; best <- list(b, a) }
    }
  }
  if (is.null(best)) list(NA_character_, NA_character_) else best
}
# clusTCR's junction_aa is BETA||ALPHA concatenated. Split by exact lookup
# against a (beta||alpha) index from the labels TSV - no guessing - with the
# heuristic only for rows the lookup misses.
.load_clustcr_paired <- function(dir, ab_pairs) {
  f <- file.path(dir, "clusTCR_output.csv"); if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  if (!all(c("junction_aa", "cluster") %in% names(d))) return(NULL)
  d2 <- if (!is.null(ab_pairs))
    left_join(d, distinct(ab_pairs, concat, cdr3b, cdr3a),
              by = c("junction_aa" = "concat"))
  else mutate(d, cdr3b = NA_character_, cdr3a = NA_character_)
  message(sprintf("  alpha LD: clusTCR label-TSV match hit %d/%d (%.1f%%)",
                  sum(!is.na(d2$cdr3a)), nrow(d2), 100 * sum(!is.na(d2$cdr3a)) / max(1, nrow(d2))))
  need <- which(is.na(d2$cdr3a))
  if (length(need)) {
    sr <- lapply(d2$junction_aa[need], .split_beta_alpha)
    d2$cdr3b[need] <- vapply(sr, `[[`, character(1), 1)
    d2$cdr3a[need] <- vapply(sr, `[[`, character(1), 2)
    message(sprintf("  alpha LD: clusTCR heuristic split hit %d/%d of remaining rows",
                    sum(!is.na(d2$cdr3a[need])), length(need)))
  }
  d2 %>% filter(!is.na(cdr3a), !is.na(cdr3b), nchar(cdr3a) >= 6, nchar(cdr3b) >= 6) %>%
    transmute(cluster_id = as.character(cluster), cdr3b, cdr3a)
}
.find_method_outputs_dir <- function(explicit = NULL) {
  if (!is.null(explicit) && dir.exists(explicit)) return(explicit)
  if (!is.null(CFG$method_outputs_dir) && dir.exists(CFG$method_outputs_dir))
    return(CFG$method_outputs_dir)
  seeds <- unique(c(CFG$jouannet_out_dir, dirname(CFG$jouannet_out_dir %||% ""),
                    dirname(dirname(CFG$jouannet_out_dir %||% "")), ".", ".."))
  seeds <- seeds[nzchar(seeds)]
  for (s in seeds) for (n in c("Output_methods", "OUTPUT_METHODS", "output_methods",
                               "methods_output", "method_outputs", "Methods_Output")) {
    p <- file.path(s, n); if (dir.exists(p)) return(p)
  }
  NULL
}

sptce_within_cluster_ld_alpha <- function(s5, method_outputs_dir = NULL,
                                          paired_methods = c("clusTCR", "GLIPH2",
                                                             "DeepTCR", "TCRdist3"),
                                          min_size = 2, fallback_to_label_lookup = TRUE) {
  lab <- tryCatch(read_tsv(CFG$jouannet_labelled, show_col_types = FALSE),
                  error = function(e) NULL)
  ab_pairs <- NULL
  if (!is.null(lab)) {
    nm <- tolower(names(lab))
    cb <- names(lab)[which(nm %in% c("cdr3_aa_beta", "cdr3b", "cdr3.beta", "beta_cdr3"))[1]]
    ca <- names(lab)[which(nm %in% c("cdr3_aa_alpha", "cdr3a", "cdr3.alpha", "alpha_cdr3"))[1]]
    if (!is.na(cb) && !is.na(ca))
      ab_pairs <- lab %>% transmute(cdr3b = .data[[cb]], cdr3a = .data[[ca]]) %>%
      filter(!is.na(cdr3b), !is.na(cdr3a), nchar(cdr3b) >= 6, nchar(cdr3a) >= 6) %>%
      mutate(concat = paste0(cdr3b, cdr3a))
  }
  mdir <- .find_method_outputs_dir(method_outputs_dir)
  if (!is.null(mdir)) message("  alpha LD: reading paired outputs from ", mdir)
  else message("  alpha LD: no method_outputs_dir found; using label-lookup fallback")
  
  pub_paired <- list()
  if (!is.null(mdir)) {
    loaders <- list(DeepTCR  = function() .load_deeptcr_paired(mdir),
                    GLIPH2   = function() .load_gliph2_paired(mdir),
                    TCRdist3 = function() .load_tcrdist3_paired(mdir),
                    clusTCR  = function() .load_clustcr_paired(mdir, ab_pairs))
    for (m in intersect(paired_methods, names(loaders))) {
      df <- tryCatch(loaders[[m]](), error = function(e) {
        message("  alpha LD: could not load ", m, ": ", conditionMessage(e)); NULL })
      if (!is.null(df) && nrow(df)) {
        pub_paired[[m]] <- df %>% mutate(method = m) %>%
          select(method, cluster_id, cdr3b, cdr3a)
        message(sprintf("  alpha LD: %-10s loaded %d rows (%d unique cdr3b, %d unique cdr3a)",
                        m, nrow(df), n_distinct(df$cdr3b), n_distinct(df$cdr3a)))
      }
    }
  }
  pub_paired_df <- bind_rows(pub_paired)
  
  missing <- setdiff(paired_methods, names(pub_paired))
  fallback_df <- tibble()
  if (length(missing) && fallback_to_label_lookup && !is.null(ab_pairs)) {
    ab_map <- ab_pairs %>% group_by(cdr3b) %>%
      summarise(cdr3a = { t <- table(cdr3a); names(sort(t, decreasing = TRUE))[1] },
                .groups = "drop")
    present <- intersect(missing, unique(as.character(s5$published_clusters$method)))
    if (length(present)) {
      fallback_df <- s5$published_clusters %>% filter(method %in% present) %>%
        mutate(method = as.character(method)) %>% select(-any_of("cdr3a")) %>%
        left_join(ab_map, by = "cdr3b") %>% filter(!is.na(cdr3a)) %>%
        select(method, cluster_id, cdr3b, cdr3a)
      message("  alpha LD: fallback label-lookup used for: ", paste(present, collapse = ", "),
              " (alpha may under-represent the method's true clustering)")
    }
  }
  sptce <- s5$sptce_clusters %>% filter(!is.na(cdr3a), nchar(cdr3a) >= 6) %>%
    transmute(method = as.character(method), cluster_id, cdr3b, cdr3a)
  combined <- bind_rows(sptce, pub_paired_df, fallback_df)
  if (!nrow(combined)) { message("  alpha LD: no usable paired data."); return(invisible(NULL)) }
  
  per_cluster <- combined %>% distinct(method, cluster_id, cdr3a) %>%
    group_by(method, cluster_id) %>% filter(n() >= min_size) %>%
    summarise(size = n(),
              median_ld = { d <- stringdistmatrix(cdr3a, cdr3a, method = "lv"); median(d[upper.tri(d)]) },
              mean_ld   = { d <- stringdistmatrix(cdr3a, cdr3a, method = "lv"); mean(d[upper.tri(d)]) },
              .groups = "drop")
  if (!nrow(per_cluster)) { message("  alpha LD: no cluster meets size >= ", min_size); return(invisible(NULL)) }
  
  sp_ord <- per_cluster %>% filter(grepl("^SPACE-TCR", method)) %>% group_by(method) %>%
    summarise(m = median(median_ld, na.rm = TRUE), .groups = "drop") %>% arrange(m) %>% pull(method)
  pub_ord <- per_cluster %>% filter(!grepl("^SPACE-TCR", method)) %>% group_by(method) %>%
    summarise(m = median(median_ld, na.rm = TRUE), .groups = "drop") %>% arrange(m) %>% pull(method)
  per_cluster <- per_cluster %>%
    mutate(method = factor(method, levels = c(rev(pub_ord), rev(sp_ord))),
           kind = case_when(grepl("SPACE-TCR RMSD\\+dih", method) ~ "SPACE-TCR RMSD+dih (focal)",
                            grepl("SPACE-TCR RMSD", method)       ~ "SPACE-TCR RMSD",
                            TRUE                                  ~ "published"),
           kind = factor(kind, levels = c("published", "SPACE-TCR RMSD",
                                          "SPACE-TCR RMSD+dih (focal)")))
  .ext_write_csv(per_cluster, "within_cluster_ld_alpha.csv")
  
  p <- ggplot(per_cluster, aes(median_ld, method, fill = kind)) +
    geom_violin(scale = "width", trim = FALSE, colour = "grey30", alpha = 0.85, linewidth = 0.3) +
    geom_boxplot(width = 0.18, fill = "white", outlier.size = 0.5, colour = "grey30", linewidth = 0.3) +
    geom_vline(xintercept = c(1, 2, 3), linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    scale_fill_manual(values = c("published" = "grey78", "SPACE-TCR RMSD" = "#F08080",
                                 "SPACE-TCR RMSD+dih (focal)" = "#B22222"), name = NULL) +
    scale_x_continuous(breaks = seq(0, 20, 2), expand = expansion(mult = c(0.01, 0.03))) +
    labs(title = "Within-cluster CDR3a LD",
         subtitle = paste0("Paired-chain methods only. clusTCR/GLIPH2/DeepTCR/TCRdist3 read from ",
                           if (!is.null(mdir)) basename(mdir) else "labels lookup (approx.)"),
         x = "Median within-cluster CDR3a Levenshtein distance", y = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          axis.title.x = element_text(face = "bold"), plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey30"))
  .ext_save("fig_within_cluster_ld_alpha.png", p, w = 10.5, h = 7)
  
  loaded <- names(pub_paired)
  summ <- per_cluster %>% group_by(method) %>%
    summarise(n_clusters = n(), median_of_medians = median(median_ld, na.rm = TRUE),
              q25 = quantile(median_ld, 0.25, na.rm = TRUE),
              q75 = quantile(median_ld, 0.75, na.rm = TRUE),
              max_ld = max(median_ld, na.rm = TRUE), .groups = "drop") %>%
    mutate(source = case_when(grepl("^SPACE-TCR", as.character(method)) ~ "sptce_clusters",
                              as.character(method) %in% loaded ~ "paired output",
                              TRUE ~ "label-lookup (fallback)"))
  .ext_write_csv(summ, "within_cluster_ld_alpha_summary.csv")
  invisible(list(per_cluster = per_cluster, summary = summ, method_outputs_dir = mdir))
}

# =================================================================
# 24. PER-METHOD EPITOPE CAPTURE
# =================================================================
sptce_epitope_capture_stats <- function(s5, hi_purity = 0.9, big_size = 5) {
  combined <- bind_rows(
    transmute(s5$sptce_clusters,     method = as.character(method), cluster_id, cdr3b, epitope),
    transmute(s5$published_clusters, method = as.character(method), cluster_id, cdr3b, epitope)
  ) %>% filter(!is.na(epitope))
  cs <- combined %>% group_by(method, cluster_id) %>%
    summarise(size = n(), dom = { t <- table(epitope); names(t)[which.max(t)] },
              n_dom = { t <- table(epitope); max(t) }, purity = n_dom / size, .groups = "drop")
  stats <- cs %>% group_by(method) %>%
    summarise(total_clusters           = n(),
              total_tcrs_clustered     = sum(size),
              n_distinct_epitopes_seen = n_distinct(dom),
              n_epitopes_hi_purity     = n_distinct(dom[purity >= hi_purity]),
              n_epitopes_hi_purity_big = n_distinct(dom[purity >= hi_purity & size >= big_size]),
              n_hi_purity_clusters     = sum(purity >= hi_purity),
              n_big_hi_purity_clusters = sum(purity >= hi_purity & size >= big_size),
              median_cluster_size      = median(size),
              max_cluster_size         = max(size), .groups = "drop") %>%
    arrange(desc(n_epitopes_hi_purity_big))
  .ext_write_csv(stats, "epitope_capture_stats.csv")
  
  big_label <- sprintf("Purity >= 0.9 & size >= %d", big_size)
  long <- stats %>%
    transmute(method, `Any cluster` = n_distinct_epitopes_seen,
              `Purity >= 0.9` = n_epitopes_hi_purity, `Big&pure` = n_epitopes_hi_purity_big) %>%
    pivot_longer(-method, names_to = "criterion", values_to = "n_epitopes") %>%
    mutate(criterion = recode(criterion, `Big&pure` = big_label),
           criterion = factor(criterion, levels = c("Any cluster", "Purity >= 0.9", big_label)))
  ord <- long %>% group_by(method) %>%
    summarise(m = max(n_epitopes, na.rm = TRUE), .groups = "drop") %>% arrange(m) %>% pull(method)
  p <- long %>% mutate(method = factor(method, levels = ord)) %>%
    ggplot(aes(method, n_epitopes, fill = criterion)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = n_epitopes), position = position_dodge(width = 0.75),
              hjust = -0.15, size = 2.7) +
    scale_fill_manual(values = c("#DEEBF7", "#6BAED6", "#08519C")) + coord_flip() +
    labs(title = "Distinct epitopes captured per method",
         subtitle = sprintf("At three cluster-quality bars (hi-purity = %.1f, big = >= %d)",
                            hi_purity, big_size),
         x = NULL, y = "Number of distinct dominant epitopes", fill = "Criterion") +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
  .ext_save("fig_epitope_capture_stats.png", p, w = 11, h = 6.5)
  invisible(stats)
}

# =================================================================
# 25. FOCAL-EPITOPE CONTAMINATION COMPOSITION
# =================================================================
# fig19_focal_cluster_composition answers "how much contamination"; this answers
# "what is it", which is what separates a clustering error from a genuine
# cross-reactive or variant-peptide relationship.
sptce_focal_epitope_contamination <- function(s5, focal_epi = CFG$focal_epitope,
                                              highlight = "NLVPMVATV",
                                              max_per_panel = 20L) {
  combined <- bind_rows(
    s5$sptce_clusters     %>% transmute(method = as.character(method), cluster_id, cdr3b, epitope),
    s5$published_clusters %>% transmute(method = as.character(method), cluster_id, cdr3b, epitope)
  ) %>% filter(!is.na(epitope))
  
  focal_clusters <- combined %>% group_by(method, cluster_id) %>%
    summarise(dom  = { t <- table(epitope); names(t)[which.max(t)] },
              size = n(), .groups = "drop") %>%
    filter(dom == focal_epi)
  if (!nrow(focal_clusters)) {
    message("  contamination: no ", focal_epi, "-dominant clusters."); return(invisible(NULL))
  }
  contam <- combined %>%
    inner_join(focal_clusters, by = c("method", "cluster_id")) %>%
    filter(epitope != focal_epi) %>%
    count(method, epitope, name = "n_tcrs")
  if (!nrow(contam)) {
    message("  contamination: no non-focal members in any ", focal_epi, " cluster.")
    return(invisible(NULL))
  }
  .ext_write_csv(contam, sprintf("focal_epitope_contamination_%s.csv", focal_epi))
  
  totals <- contam %>% group_by(method) %>%
    summarise(n_contaminating_tcrs    = sum(n_tcrs),
              n_distinct_contaminants = n(),
              top_contaminant         = epitope[which.max(n_tcrs)],
              .groups = "drop") %>%
    arrange(n_contaminating_tcrs)
  .ext_write_csv(totals, sprintf("focal_epitope_contamination_totals_%s.csv", focal_epi))
  cat("\n  Contamination of ", focal_epi, " clusters by method:\n\n", sep = "")
  print(as.data.frame(totals), row.names = FALSE); cat("\n")
  
  # Facets have independent y axes, so one global factor order would be wrong in
  # every panel but one. Key each bar by method+epitope, order that, strip the
  # prefix at render time.
  plot_df <- contam %>% group_by(method) %>%
    slice_max(n_tcrs, n = max_per_panel, with_ties = FALSE) %>% ungroup() %>%
    mutate(is_hl  = epitope == highlight,
           method = order_methods(method),
           ykey   = fct_reorder(paste(method, epitope, sep = "___"), n_tcrs))
  p <- ggplot(plot_df, aes(n_tcrs, ykey, fill = is_hl)) +
    geom_col(width = 0.75) +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_fill_manual(values = c(`FALSE` = "#6BAED6", `TRUE` = "#D73027"), guide = "none") +
    facet_wrap(~ method, scales = "free", ncol = 3,
               labeller = labeller(method = sptce_label)) +
    labs(title = sprintf("Top contaminating epitopes in %s clusters", focal_epi),
         subtitle = sprintf("Non-%s TCRs inside %s-dominant clusters, by their own epitope. %s in red. Note the free x-axis.",
                            focal_epi, focal_epi, highlight),
         x = "Number of non-focal TCRs", y = NULL) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_text(size = 7),
          panel.grid.major.y = element_blank())
  .ext_save(sprintf("fig19_focal_epitope_contam_%s.png", focal_epi), p, w = 14, h = 10)
  invisible(list(contamination = contam, totals = totals))
}

# =================================================================
# 26. EXTENSIONS ORCHESTRATOR       [replaces run_stage11_extensions]
# =================================================================
run_stage11_extensions <- function(s5 = NULL) {
  cat("\n=== Stage 11: paper-response extensions ===\n")
  if (is.null(s5)) {
    message("needs s5 from run_stage5(); running it now.")
    s5 <- run_stage5()
    if (is.null(s5)) { message("run_stage5() returned NULL; aborting."); return(invisible(NULL)) }
  }
  miss <- setdiff(c("sptce_clusters", "published_clusters", "total_n"), names(s5))
  if (length(miss)) {
    message("s5 missing: ", paste(miss, collapse = ", "), " — aborting."); return(invisible(NULL))
  }
  step <- function(label, expr) {
    cat("\n--", label, "----------------------------------\n")
    tryCatch(expr, error = function(e) { message("  ", label, " failed: ", conditionMessage(e)); NULL })
  }
  gc_    <- step("V/J gene consensus",        sptce_gene_consensus(s5))
  he_b   <- step("CDR3b position entropy",    sptce_position_entropy_heatmap(s5, chain = "beta"))
  he_a   <- step("CDR3a position entropy",    sptce_position_entropy_heatmap(s5, chain = "alpha"))
  funnel <- step("counting funnel",           sptce_counting_funnel(s5))
  cons   <- step("consensus clustering",      sptce_consensus_clusters(s5))
  ranges <- step("all-epitope ranges",        sptce_all_epitope_ranges(s5))
  caps   <- step("epitope-capture stats",     sptce_epitope_capture_stats(s5))
  ld_a    <- step("within-cluster CDR3a LD",   sptce_within_cluster_ld_alpha(s5))
  contam <- step("focal-epitope contamination", sptce_focal_epitope_contamination(s5))
  audit  <- step("networks overview audit",   sptce_networks_overview_audit(s5))
  leak   <- step("leakage, both denominators", sptce_leakage_full_cohort())
  
  cat("\n=== Extensions complete. Outputs in: ", .ext_dir(), " ===\n", sep = "")
  invisible(list(gene_consensus = gc_, entropy_beta = he_b, entropy_alpha = he_a,
                 funnel = funnel, consensus = cons, epitope_ranges = ranges,
                 epitope_capture = caps, ld_alpha = ld_a,
                 focal_contamination = contam, networks_audit = audit,
                 leakage_both = leak))
}



# =================================================================
# PREFLIGHT: what structural output exists
# =================================================================
# Benchmark A needs one structural run, not one per rung: it clusters the
# largest available run and scores each rung by restricting that pair list to
# the rung's receptors. .detect_ref_subset() uses the subset*_clusters folder
# covering the most unique TCRs as the reference.
#
#   1. The distance tables are keyed to the subset definitions they were built
#      from. If a _pre_ladder_fix/ folder is present, .paired_betas() reads the
#      reference subset's input from there, so the receptor ids match.
#
#   2. Without a subsetN_dihedrals folder matching the reference run,
#      Benchmark A runs RMSD-only (.find_benchA_dihedrals() finds nothing).
#
#      A dihedral table must come from the same run as the distance tables.
#      Pairs where either receptor lacks angles are excluded from the gated
#      modes, so a table from a different run gives near-zero coverage and
#      empty gated clusters. .load_ref_pairs() reports the coverage and warns
#      below 50%; sptce_preflight(check_dihedral_ids = TRUE) checks it before
#      a long run.
sptce_preflight <- function(root = CFG$sptce_root,
                            check_dihedral_ids = TRUE,
                            n_sample_files = 20L) {
  if (!dir.exists(root)) { message("preflight: root not found: ", root); return(invisible(NULL)) }
  dirs <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  nums <- sort(unique(as.integer(sub("^subset(\\d+)_.*$", "\\1",
                                     grep("^subset\\d+_(clusters|dihedrals)$", dirs, value = TRUE)))))
  if (!length(nums)) { message("preflight: no subset*_clusters directories under ", root); return(invisible(NULL)) }
  
  tab <- purrr::map_dfr(nums, function(i) {
    cdir <- file.path(root, sprintf("subset%d_clusters", i))
    ddir <- file.path(root, sprintf("subset%d_dihedrals", i))
    # read_tsv() reads either, so accept both spellings
    dfil <- file.path(ddir, c("dihedrals.tsv.gz", "dihedrals.tsv"))
    tibble(subset = i,
           clusters_dir = dir.exists(cdir),
           n_distance_files = if (dir.exists(file.path(cdir, "distance_tables")))
             length(list.files(file.path(cdir, "distance_tables"), "\\.tsv\\.gz$")) else 0L,
           dihedrals_dir = dir.exists(ddir),
           dihedral_file = any(file.exists(dfil)))
  })
  cat("\n  Structural runs under ", root, ":\n\n", sep = "")
  print(as.data.frame(tab), row.names = FALSE)
  
  usable <- filter(tab, n_distance_files > 0)
  if (!nrow(usable)) {
    message("\n  No distance tables found - Benchmark A cannot run.")
    return(invisible(tab))
  }
  # file count is a proxy for size; .detect_ref_subset() does the exact count
  ref <- usable$subset[which.max(usable$n_distance_files)]
  cat("\n  Likely reference run: subset ", ref,
      " (most distance tables). One run is enough for the whole ladder.\n", sep = "")
  if (!usable$dihedral_file[usable$subset == ref]) {
    cat("  NOTE: no dihedrals.tsv.gz for subset ", ref,
        " - Benchmark A will run RMSD-only and the\n",
        "        dihedral gate will not be assessed against repertoire noise.\n",
        "        To enable it, generate dihedrals for THAT run, or set\n",
        "        CFG$benchA_dih_path. Do not substitute another rung's table.\n", sep = "")
  }
  pre <- file.path(root, "_pre_ladder_fix")
  cat("  _pre_ladder_fix present: ", dir.exists(pre),
      " (reference-run inputs are read from here when it exists)\n", sep = "")
  
  # Do the dihedral ids actually match the distance tables of the reference run?
  # A mismatch yields empty gated clusters rather than an error, so this is
  # worth a few seconds now rather than a puzzling result after a long run.
  if (isTRUE(check_dihedral_ids) && any(usable$subset == ref) &&
      usable$dihedral_file[usable$subset == ref]) {
    dfil <- file.path(root, sprintf("subset%d_dihedrals", ref), "dihedrals.tsv.gz")
    ddir <- file.path(root, sprintf("subset%d_clusters", ref), "distance_tables")
    fs <- head(list.files(ddir, "\\.tsv\\.gz$", full.names = TRUE), n_sample_files)
    dih_ids <- tryCatch(
      unique(read_tsv(dfil, show_col_types = FALSE, col_select = "tcr_id")$tcr_id),
      error = function(e) { message("  preflight: could not read ", dfil, ": ",
                                    conditionMessage(e)); NULL })
    if (!is.null(dih_ids) && length(fs)) {
      pair_ids <- unique(unlist(lapply(fs, function(f) {
        dt <- tryCatch(read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
        if (is.null(dt)) return(character(0))
        if (all(c("id_a", "id_b") %in% names(dt)) && !"id_i" %in% names(dt))
          dt <- rename(dt, id_i = id_a, id_j = id_b)
        clean_id(unique(c(dt$id_i, dt$id_j)))
      })))
      cov <- mean(pair_ids %in% dih_ids)
      cat(sprintf("  dihedral id coverage on subset %d: %.1f%% of %d sampled ids (table holds %d)\n",
                  ref, 100 * cov, length(pair_ids), length(dih_ids)))
      if (cov < 0.5)
        cat("  WARNING: coverage below 50%. The dihedral table is probably keyed to a\n",
            "           different subset definition than the distance tables; gated modes\n",
            "           would come out near-empty. Check it was built from the same input\n",
            "           (see _pre_ladder_fix/).\n", sep = "")
      else if (cov < 0.95)
        cat("  Note: partial coverage. Pairs with a missing partner are dropped, so gated\n",
            "        retention will be understated by roughly that fraction.\n", sep = "")
    }
  }
  invisible(tab)
}


# =================================================================
# BENCHMARK A: DIHEDRAL DISCOVERY
# =================================================================
#
# Finds dihedrals.tsv.gz or an uncompressed dihedrals.tsv, and reports how many
# of the reference run's TCR ids the table covers; a table from another cohort
# (Benchmark B's, say) gives ~0% coverage.
.find_benchA_dihedrals <- function(ref_subset, dih_path = CFG$benchA_dih_path,
                                   ref_ids = NULL) {
  stems <- c(file.path(CFG$sptce_root, sprintf("subset%d_dihedrals/dihedrals", ref_subset)),
             file.path(CFG$sptce_root, sprintf("subset%d_dihedrals", ref_subset)),
             file.path(CFG$sptce_root, sprintf("subset%d_clusters/dihedrals", ref_subset)),
             file.path(CFG$sptce_root, "dihedrals"))
  cands <- c(dih_path,
             unlist(lapply(stems, function(x) paste0(x, c(".tsv.gz", ".tsv")))))
  cands <- cands[!is.na(cands) & nzchar(cands)]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) {
    cat("  [note] no dihedral table found for the reference run. Benchmark A\n")
    cat("         will run RMSD-only. Searched:\n")
    for (p in cands) cat("           ", p, "\n", sep = "")
    cat("         Set CFG$benchA_dih_path if it lives elsewhere. Do NOT use\n")
    cat("         Benchmark B's file: it indexes a different cohort's TCR ids.\n")
    return(NULL)
  }
  cat("  dihedrals: ", hit, "\n", sep = "")
  dmap <- load_dihedral_map(hit)
  cat("  dihedral table covers ", length(dmap), " TCRs\n", sep = "")
  if (!is.null(ref_ids) && length(ref_ids)) {
    ov <- length(intersect(names(dmap), unique(ref_ids))) / length(unique(ref_ids))
    cat(sprintf("  overlap with reference-run TCR ids: %.1f%%\n", 100 * ov))
    if (ov < 0.5)
      warning(sprintf(paste0("Dihedral table matches only %.1f%% of the reference run's ",
                             "TCR ids. Dihedral-gated modes will look far more restrictive ",
                             "than they are, because an edge is dropped whenever either ",
                             "partner is missing. Check this file belongs to the same run."),
                      100 * ov), call. = FALSE)
  }
  dmap
}

# #################################################################
##  PART 3 - ENTRY POINTS
# #################################################################
#
# HuARdb (Benchmark C) is not run by main(). The disease-level analysis in the
# paper uses the Observed TCR Space (Raybould et al., 2024), a different
# dataset that this pipeline does not produce.

run_benchmark_b <- function(extensions = TRUE) {
  cat("\n=== Benchmark B (Jouannet curated cohort) ===\n")
  run_stage7_training_leakage()   # first: writes the file Stage 5d reads
  s5 <- run_stage5()
  run_stage5b_networks(s5)
  run_stage5c_extras(s5)
  run_stage5d_cohort_accounting(s5)
  run_stage6(s5)
  run_stage6b_cohort_diff(s5)
  if (isTRUE(extensions)) run_stage11_extensions(s5)
  invisible(s5)
}

main <- function(extensions = TRUE) {
  cat("\n=================================================================\n")
  cat("  SPACE-TCR pipeline - started:", format(Sys.time()), "\n")
  cat("=================================================================\n")
  run_stage1()
  run_benchmark_a()
  run_benchmark_b(extensions = extensions)
  cat("\n=================================================================\n")
  cat("  Pipeline complete:", format(Sys.time()), "\n")
  cat("=================================================================\n")
}


# =================================================================
# PART 2 — BENCHMARK A AND SHARED HELPERS
# =================================================================
# Stage functions, Benchmark A, and the figure and overlay helpers they use.
# run_stage10_huardb, .huardb_disease_enrichment and normalise_disease are
# defined but not reachable from main().


# =================================================================
# 9. JOUANNET DATABASE CURATION
# =================================================================
load_jouannet_curated <- function(db_path, cell_type = "CD8", verbose = TRUE) {
  if (verbose) cat("Loading Jouannet database:", db_path, "\n")
  raw <- read_delim(db_path, delim="\t", show_col_types=FALSE, guess_max=100000)
  if (verbose) cat("  Raw rows:", nrow(raw), "\n")
  pick <- function(df, cands) { f <- intersect(cands, names(df)); if(length(f)) f[1] else NA_character_ }
  cols <- list(
    cdr3b    = pick(raw, c("CDR3_beta","cdr3_beta","CDR3b","cdr3b")),
    cdr3a    = pick(raw, c("CDR3_alpha","cdr3_alpha","CDR3a","cdr3a")),
    v_beta   = pick(raw, c("V_beta","v_beta","Vbeta","TRBV")),
    j_beta   = pick(raw, c("J_beta","j_beta","Jbeta","TRBJ")),
    v_alpha  = pick(raw, c("V_alpha","v_alpha","Valpha","TRAV")),
    j_alpha  = pick(raw, c("J_alpha","j_alpha","Jalpha","TRAJ")),
    epitope  = pick(raw, c("Epitope","epitope","Antigen","antigen")),
    vs       = pick(raw, c("Verified_score","verified_score","VS","Verified")),
    ais      = pick(raw, c("Identification_score","identification_score","AIS","Antigen_score")),
    celltype = pick(raw, c("Cell_subset","cell_subset","Cell_type","CellType","T_cell_type")),
    source   = pick(raw, c("PubMed_ID","pubmed_id","PubMedID","Source","source","Database"))
  )
  if (is.na(cols$cdr3b) || is.na(cols$epitope)) stop("Cannot find CDR3b or Epitope column")
  get_col <- function(col, type_fn) if (!is.na(col)) type_fn(raw[[col]]) else NA
  df <- raw %>% transmute(
    cdr3b          = as.character(.data[[cols$cdr3b]]),
    cdr3a          = get_col(cols$cdr3a,  as.character),
    v_beta         = get_col(cols$v_beta,  as.character),
    j_beta         = get_col(cols$j_beta,  as.character),
    v_alpha        = get_col(cols$v_alpha, as.character),
    j_alpha        = get_col(cols$j_alpha, as.character),
    epitope        = as.character(.data[[cols$epitope]]),
    verified_score = get_col(cols$vs,      \(x) suppressWarnings(as.numeric(x))),
    ais            = get_col(cols$ais,     \(x) suppressWarnings(as.numeric(x))),
    cell_subset    = get_col(cols$celltype, as.character),
    pubmed_id      = get_col(cols$source,   as.character)
  )
  n0 <- nrow(df)
  if (!is.na(cols$vs)) df <- filter(df, !is.na(verified_score), verified_score == 2)
  if (!is.na(cols$ais)) df <- filter(df, !is.na(ais), ais > 4.3)
  df <- filter(df, !is.na(epitope), epitope != "", epitope != "NA")
  ten_x <- "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf"
  if (!is.na(cols$source)) df <- filter(df, coalesce(pubmed_id,"") != ten_x)
  if (!is.na(cols$celltype) && cell_type != "both")
    df <- filter(df, cell_subset == cell_type)
  df <- filter(df, !is.na(v_beta), v_beta != "", !is.na(j_beta), j_beta != "",
               !is.na(v_alpha), v_alpha != "", !is.na(j_alpha), j_alpha != "")
  df <- filter(df, !is.na(cdr3b), nchar(cdr3b) >= 6, nchar(cdr3b) <= 23) %>%
    filter(is.na(cdr3a) | cdr3a == "" | (nchar(cdr3a) >= 6 & nchar(cdr3a) <= 23))
  keep_ep <- df %>%
    mutate(pair_id = ifelse(is.na(cdr3a)|cdr3a=="", cdr3b, paste0(cdr3b,"_",cdr3a))) %>%
    distinct(epitope, pair_id) %>% count(epitope) %>% filter(n >= 2) %>% pull(epitope)
  df <- filter(df, epitope %in% keep_ep)
  bl <- c("CASSDSRGTEAFF_CASSDSRGTEAFF","CASSPITGTGAYGYTF_CASSSVNEQYF")
  df <- df %>% mutate(.pair = paste0(cdr3b,"_",cdr3a)) %>%
    filter(!.pair %in% bl) %>% select(-.pair)
  df <- df %>% arrange(desc(verified_score), desc(ais)) %>% distinct(cdr3b, .keep_all=TRUE)
  if (verbose) cat("  Final:", nrow(df), "(started:", n0, ")\n")
  df
}


# =================================================================
# 10. STAGE 1 — 10X contigs -> paired master
# =================================================================
# Still required: Benchmark A matches the Jouannet 10X betas against the paired
# receptors in this master file, and the labels are rebuilt from the same
# dextramer matrix.
run_stage1 <- function() {
  cat("\n=== Stage 1: 10X contigs -> paired master ===\n")
  if (!file.exists(CFG$contig_file)) { message("Contig file missing"); return(invisible(NULL)) }
  contigs <- read_csv(CFG$contig_file, show_col_types=FALSE) %>%
    mutate(across(c(is_cell, high_confidence, productive), to_lgl)) %>%
    filter(is_cell, high_confidence, productive,
           !is.na(cdr3), cdr3 != "None", nchar(cdr3) >= 6, nchar(cdr3) <= 22)
  valid_bc <- contigs %>% filter(chain %in% c("TRA","TRB")) %>%
    count(barcode, chain) %>%
    pivot_wider(names_from=chain, values_from=n, values_fill=0) %>%
    filter(TRB==1, TRA>=1, TRA<=2) %>% pull(barcode)
  contigs <- filter(contigs, barcode %in% valid_bc)
  beta  <- filter(contigs, chain=="TRB") %>%
    select(barcode, v_call_beta=v_gene, j_call_beta=j_gene, cdr3_aa_beta=cdr3,
           reads_beta=reads, umis_beta=umis)
  alpha <- filter(contigs, chain=="TRA") %>%
    arrange(barcode, desc(reads), desc(umis)) %>% group_by(barcode) %>% slice(1) %>% ungroup() %>%
    select(barcode, v_call_alpha=v_gene, j_call_alpha=j_gene, cdr3_aa_alpha=cdr3,
           reads_alpha=reads, umis_alpha=umis)
  paired <- inner_join(alpha, beta, by="barcode") %>%
    filter(!is.na(v_call_alpha), !is.na(j_call_alpha), !is.na(v_call_beta), !is.na(j_call_beta)) %>%
    mutate(across(c(v_call_alpha,j_call_alpha,v_call_beta,j_call_beta), add_allele),
           id = paste0("tcr10x_", row_number()))
  cat("  Paired cells:", nrow(paired), "\n")
  binmat <- read_csv(CFG$binmat_file, show_col_types=FALSE)
  if ("donor" %in% names(binmat)) binmat <- filter(binmat, donor=="donor1")
  meta_cols <- c("barcode","donor","cell_clono_cdr3_aa","cell_clono_cdr3_nt")
  surf_cols <- c("CD3","CD19","CD45RA","CD4","CD8a","CD14","CD45RO","CD279_PD-1",
                 "IgG1","IgG2a","IgG2b","CD127","CD197_CCR7","HLA-DR")
  dex_cols  <- setdiff(names(binmat), c(meta_cols, surf_cols, grep("_binder$", names(binmat), value=TRUE)))
  dex_long <- binmat %>% select(barcode, all_of(dex_cols)) %>%
    pivot_longer(-barcode, names_to="epitope", values_to="dex_umi") %>%
    filter(dex_umi > 10, !str_detect(epitope, "_NC$")) %>%
    group_by(barcode) %>% slice_max(dex_umi, n=1, with_ties=FALSE) %>% ungroup()
  master <- paired %>% left_join(dex_long, by="barcode") %>%
    select(id, barcode, v_call_alpha, j_call_alpha, cdr3_aa_alpha,
           v_call_beta,  j_call_beta,  cdr3_aa_beta,
           reads_alpha, umis_alpha, reads_beta, umis_beta, epitope, dex_umi)
  write_tsv(master, CFG$master_file)
  write_tsv(select(master, id, v_call_alpha, j_call_alpha, cdr3_aa_alpha,
                   v_call_beta, j_call_beta, cdr3_aa_beta), CFG$sptce_input)
  cat("  Written:", CFG$master_file, "\n")
  invisible(master)
}


# #################################################################
##
##   BENCHMARK A — paired-only comparison on the Jouannet 10X subsets
##
##   run_benchmark_a()
##
##   Jouannet ran all nine published methods on the 10X donor-1 noise subsets,
##   so nothing needs re-running: SPACE-TCR and the published methods are scored
##   on the SAME subsets, SAME labels, SAME purity definition and SAME minimum
##   cluster size, restricted to the betas for which a paired structure exists.
##
##   FOUR FACTS FROM THEIR OWN SCRIPTS THAT SHAPE THIS
##     1. The ladder is ARITHMETIC: subsets of 2876 + 4737k, noise
##        0 / 62.2 / 76.7 / 83.2 / 86.8 / 89.2%. (7613/2876 is their SUBSET_2
##        ratio, not their endpoint.)
##     2. The unit is unique CDR3-beta, not paired TCRs; their Gliph2 input
##        even sets CDR3a = NA. SPACE-TCR is collapsed to beta level, a beta
##        counting as clustered if ANY of its receptors clustered.
##     3. Their .csv subset inputs are TAB-delimited despite the extension.
##     4. Labels are multi-epitope: a beta may bind several dextramers, and
##        their purity counts the dominant epitope's rows over the ORIGINAL
##        cluster size, so unlabelled members dilute purity without ever
##        adding to it.
##
##   THE CAVEAT TO STATE IN THE WRITE-UP
##     The published methods CLUSTERED on the full unpaired cohort and are
##     EVALUATED here on the paired subset (~85% of betas). Their outputs are
##     filtered afterwards, so a cluster of five may drop to three. That
##     equalises the evaluation universe, not the clustering input — re-running
##     each tool on paired input is impossible without the binaries. The
##     `kept_frac` column reports how much of each method's output survived.
##     SPACE-TCR is unaffected: it only ever saw paired receptors.
##
# #################################################################

# -----------------------------------------------------------------
# A1. The subsets and their true noise axis
# -----------------------------------------------------------------
jouannet_10x_ladder <- function(root = CFG$tenx_root, quiet = FALSE) {
  subs <- map(1:6, function(i) {
    f <- file.path(root, "subsets", sprintf("subset_%d.txt", i))
    if (!file.exists(f)) return(NULL)
    d <- suppressMessages(read_delim(f, delim = ",", show_col_types = FALSE))
    cn <- names(d)[grepl("cdr3", names(d), ignore.case = TRUE)][1]
    if (is.na(cn)) return(NULL)
    unique(as.character(d[[cn]]))
  })
  names(subs) <- paste0("subset_", 1:6)
  ok <- !map_lgl(subs, is.null)
  if (!any(ok)) {
    message("No subset files under ", file.path(root, "subsets"),
            " — set CFG$tenx_root.")
    return(NULL)
  }
  n1 <- length(subs[[1]])
  tab <- tibble(subset = which(ok), n = map_int(subs[ok], length)) %>%
    mutate(pct_noise = round(100 * (n - n1) / n, 1), increment = c(NA, diff(n)))
  if (!quiet) {
    cat("\n=== Jouannet 10X ladder ===\n\n")
    print(as.data.frame(tab), row.names = FALSE)
    if (n_distinct(na.omit(tab$increment)) == 1)
      cat("\n  Constant increment of ", na.omit(tab$increment)[1],
          " => arithmetic ladder, noise to ", max(tab$pct_noise), "%.\n", sep = "")
  }
  invisible(list(subsets = subs[ok], table = tab,
                 noise = setNames(tab$pct_noise, as.character(tab$subset))))
}


# -----------------------------------------------------------------
# A2. Labels
# -----------------------------------------------------------------
# Their labelled_cdr3b.csv write is commented out in 10X_processing.R, so it is
# usually absent; the reconstruction follows their method (mean dextramer score
# per CDR3b across the cells carrying it, binder if > 10, ALL binders kept).
#
# align_to_subset1: their subset_1 IS their annotated set (2,876 betas, every
# one a binder by construction). The reconstruction returns rather more, because
# 10X_processing.R applies contig filters this pipeline does not. Restricting to
# subset_1 adopts the benchmark's own definition of "labelled" and prints the
# residual disagreement instead of absorbing it.
jouannet_10x_labels <- function(root = CFG$tenx_root, threshold = 10,
                                align_to_subset1 = TRUE) {
  f <- file.path(root, "labelled_cdr3b.csv")
  lab <- NULL
  if (file.exists(f)) {
    cat("\n=== Labels: their labelled_cdr3b.csv ===\n")
    m <- suppressMessages(read_csv(f, show_col_types = FALSE))
    cn <- names(m)[grepl("cdr3", names(m), ignore.case = TRUE)][1]
    if (!is.na(cn))
      lab <- m %>% pivot_longer(-all_of(cn), names_to = "epitope", values_to = "v") %>%
        filter(!is.na(v), v != 0) %>%
        transmute(cdr3b = as.character(.data[[cn]]), epitope)
  }
  if (is.null(lab)) {
    cat("\n=== Labels: rebuilding from the 10x dextramer matrix ===\n")
    if (!file.exists(CFG$binmat_file) || !file.exists(CFG$contig_file)) {
      message("  need CFG$binmat_file and CFG$contig_file"); return(NULL)
    }
    contigs <- read_csv(CFG$contig_file, show_col_types = FALSE) %>%
      mutate(across(any_of(c("is_cell","high_confidence","productive")), to_lgl)) %>%
      filter(is_cell, high_confidence, productive, chain == "TRB",
             !is.na(cdr3), cdr3 != "None", nchar(cdr3) >= 6, nchar(cdr3) <= 22) %>%
      select(barcode, cdr3)
    bm <- read_csv(CFG$binmat_file, show_col_types = FALSE)
    if ("donor" %in% names(bm)) bm <- filter(bm, donor == "donor1")
    meta <- c("barcode","donor","cell_clono_cdr3_aa","cell_clono_cdr3_nt")
    surf <- c("CD3","CD19","CD45RA","CD4","CD8a","CD14","CD45RO","CD279_PD-1",
              "IgG1","IgG2a","IgG2b","CD127","CD197_CCR7","HLA-DR")
    dex <- setdiff(names(bm), c(meta, surf, grep("_binder$", names(bm), value = TRUE)))
    long <- bm %>% select(barcode, all_of(dex)) %>%
      pivot_longer(-barcode, names_to = "epitope", values_to = "score") %>%
      filter(!str_detect(epitope, "_NC$"))
    lab <- contigs %>%
      inner_join(long, by = "barcode", relationship = "many-to-many") %>%
      group_by(cdr3, epitope) %>%
      summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
      filter(mean_score > threshold) %>%
      transmute(cdr3b = cdr3, epitope)
  }
  cat("  ", n_distinct(lab$cdr3b), " annotated CDR3b, ",
      n_distinct(lab$epitope), " epitopes\n", sep = "")

  if (align_to_subset1) {
    lad <- jouannet_10x_ladder(root, quiet = TRUE)
    if (!is.null(lad)) {
      s1 <- lad$subsets[[1]]
      inside  <- intersect(unique(lab$cdr3b), s1)
      outside <- setdiff(unique(lab$cdr3b), s1)
      missing <- setdiff(s1, unique(lab$cdr3b))
      in_noise <- sum(outside %in% unlist(lad$subsets[-1]))
      cat("\n  Alignment to their annotated set (subset_1, n = ", length(s1), "):\n", sep = "")
      cat("    mine inside subset_1   : ", length(inside), "\n", sep = "")
      cat("    mine outside           : ", length(outside), " (", in_noise,
          " of them in subsets 2-6, i.e. their NOISE)\n", sep = "")
      cat("    theirs unlabelled by me: ", length(missing),
          " (dilute purity as if unlabelled, uniformly across all methods)\n", sep = "")
      lab <- filter(lab, cdr3b %in% s1)
      cat("    -> labels restricted to subset_1: ", n_distinct(lab$cdr3b), " betas\n", sep = "")
    }
  }
  lab
}


# -----------------------------------------------------------------
# A3. Parsers, transcribed from their 10X_output_processing.R
# -----------------------------------------------------------------
# Delimiters and filenames matter and are NOT guessable: TCRdist3 is ';',
# HD/LD/iSMART/GIANA are tab, GLIPH2/clusTCR/DeepTCR are ','. TCRMatch emits a
# pairwise edge list rather than a cluster table, and their fread(header=FALSE)
# auto-detects its delimiter. GIANA writes '##' header lines that their
# read.table skips by default.
.j_read <- function(f, ...) tryCatch(
  suppressMessages(read_delim(f, show_col_types = FALSE, ...)), error = \(e) NULL)


.sniff_delim <- function(f) {
  l <- tryCatch(readLines(f, n = 5, warn = FALSE), error = \(e) character(0))
  l <- l[nzchar(l) & !startsWith(l, "#")]
  if (!length(l)) return("\t")
  cnt <- c("\t" = str_count(l[1], "\t"), "," = str_count(l[1], ","),
           ";" = str_count(l[1], ";"), " " = str_count(l[1], " "))
  names(which.max(cnt))
}


# -----------------------------------------------------------------
# A4. Reference structural run: pairs and dihedrals
# -----------------------------------------------------------------
.dist_dir <- function(i) file.path(CFG$sptce_root,
                                   sprintf("subset%d_clusters/distance_tables", i))


.bin_membership <- function(i, max_files = Inf) {
  d <- .dist_dir(i)
  if (!dir.exists(d)) return(NULL)
  fs <- list.files(d, "\\.tsv\\.gz$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  if (is.finite(max_files)) fs <- head(fs, max_files)
  map_dfr(fs, function(f) {
    dt <- tryCatch(read_tsv(f, show_col_types = FALSE), error = \(e) NULL)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    if (all(c("id_a", "id_b") %in% names(dt)) && !"id_i" %in% names(dt))
      dt <- rename(dt, id_i = id_a, id_j = id_b)
    tibble(tcr_id = clean_id(unique(c(dt$id_i, dt$id_j))), bin = basename(f))
  }) %>% distinct(tcr_id, .keep_all = TRUE)
}


.detect_ref_subset <- function() {
  avail <- Filter(\(i) dir.exists(.dist_dir(i)), 1:12)
  if (!length(avail)) return(NULL)
  counts <- map_int(avail, \(i) { m <- .bin_membership(i); if (is.null(m)) 0L else nrow(m) })
  ref <- avail[which.max(counts)]
  cat("  reference run: subset", ref, "covering", max(counts), "unique TCRs\n")
  ref
}


# Load the reference pairs ONCE, with the bin tag and (when a dihedral map is
# supplied) a precomputed pairwise dihedral distance. Precomputing turns every
# (threshold x dihedral cutoff) combination into a cheap column filter instead
# of a fresh circ_dihedral() sweep over the whole pair list.
.load_ref_pairs <- function(ref_subset, rmsd_max = max(CFG$thresholds),
                            dih_map = NULL, with_dih = TRUE) {
  fs <- list.files(.dist_dir(ref_subset), "\\.tsv\\.gz$", full.names = TRUE)
  cat("  loading", length(fs), "distance files (rmsd <=", rmsd_max, ")...\n")
  use_dih <- with_dih && !is.null(dih_map)
  if (use_dih) cat("  precomputing pairwise dihedral distances...\n")
  pairs <- read_pair_tables(fs, rmsd_max, if (use_dih) dih_map else NULL)
  cat("  pairs retained:", nrow(pairs), "\n")
  if (use_dih && nrow(pairs)) {
    cov <- mean(!is.na(pairs$dih))
    cat("  dihedral coverage:", sprintf("%.1f%%", 100 * cov), "of pairs\n")
    if (cov < 0.5)
      warning(sprintf(paste0("Only %.1f%% of pairs have dihedral angles for BOTH ",
                             "partners. Dihedral-gated modes will look far more ",
                             "restrictive than they are, because a missing angle ",
                             "drops the edge. Check that the dihedral table covers ",
                             "the same TCR ids as the distance tables."), 100 * cov),
              call. = FALSE)
  } else if (nrow(pairs)) pairs$dih <- NA_real_
  pairs
}


# Benchmark A clusters from the preloaded pair table (cluster_pairs), restricted
# to `keep_ids` before components are formed, with the column names this file
# uses. Cluster ids are "<prefix>_<bin>_<component>", as build_sptce_clusters().
.cluster_from_pairs <- function(pairs, keep_ids, rmsd_cut, dih_cut, prefix,
                                labels, keep_unlabelled = TRUE) {
  res <- cluster_pairs(pairs, rename(labels, receptor_id = id), rmsd_cut,
                       dih_cut = if (is.null(dih_cut)) NA_real_ else dih_cut,
                       keep_ids = keep_ids, keep_unlabelled = keep_unlabelled,
                       id_fmt = function(bin, comp) paste0(prefix, "_", bin, "_", comp))
  if (!nrow(res)) return(res)
  rename(res, tcr_id = receptor_id)
}


# ---------------------------------------------------------------------------
# rung_table() -- the ladder from one table plus a membership column
# ---------------------------------------------------------------------------
# The ladder is nested: every rung is the last one restricted by receptor id,
# with the same epitope labels. Only the last rung is stored, as
# donor1_subset6_sptce.tsv, with ladder_membership.tsv giving each receptor the
# first rung it appears in; rung N is every receptor whose first_rung is <= N
# (3,056 / 8,089 / 21,413 / 23,529 / 23,529 / 23,529 receptors).
#
# Six per-rung files are also read if present.
.rung_table <- function(rung, dir = CFG$subset_dir) {
  direct <- c(file.path(dir, "_pre_ladder_fix",
                        sprintf("donor1_subset%d_sptce.tsv", rung)),
              file.path(dir, sprintf("donor1_subset%d_sptce.tsv", rung)))
  hit <- which(file.exists(direct))
  if (length(hit)) return(read_tsv(direct[hit[1]], show_col_types = FALSE))

  full <- file.path(dir, "donor1_subset6_sptce.tsv")
  memb <- file.path(dir, "ladder_membership.tsv")
  if (file.exists(full) && file.exists(memb)) {
    f <- read_tsv(full, show_col_types = FALSE)
    m <- read_tsv(memb, show_col_types = FALSE)
    keep <- m$id[m$first_rung <= rung]
    out <- f %>% filter(id %in% keep)
    if (!nrow(out))
      stop(".rung_table(): rung ", rung, " is empty. ", basename(memb),
           " and ", basename(full), " do not share any ids.", call. = FALSE)
    return(out)
  }
  if (!is.null(CFG$master_file) && file.exists(CFG$master_file))
    return(read_tsv(CFG$master_file, show_col_types = FALSE))
  stop(".rung_table(): no ladder table for rung ", rung, " under ", dir,
       ".\nExpected either donor1_subset", rung, "_sptce.tsv, or ",
       "donor1_subset6_sptce.tsv\nbeside ladder_membership.tsv.", call. = FALSE)
}


# Betas for which a paired structure exists.
.paired_betas <- function(ref_subset = NULL) {
  if (is.null(ref_subset)) ref_subset <- .detect_ref_subset()
  if (is.null(ref_subset)) return(NULL)
  # The reference rung's own table, not CFG$master_file: the structural run's
  # TCR ids are keyed to THAT rung's definition and not to any later one. See
  # the rung 4/5/6 note in docs/BENCHMARK_A.md.
  s <- .rung_table(ref_subset) %>%
    filter(!is.na(cdr3_aa_beta), nchar(cdr3_aa_beta) >= 6)
  list(betas = unique(s$cdr3_aa_beta), map = select(s, id, cdr3b = cdr3_aa_beta),
       ref_subset = ref_subset)
}


# Purity on unique CDR3b where a CDR3b may carry several epitope annotations.
.purity_multi <- function(clu, lab, min_size = METRIC$min_cluster_size) {
  d <- clu %>% left_join(lab, by = "cdr3b", relationship = "many-to-many")
  jm_purity(cluster_table(d, unit = "cdr3b", min_size = min_size))
}


.paired_score <- function(clu, lab, n_labelled, min_size = 2) {
  empty <- tibble(n_clusters = 0L, n_clustered = 0L, median_size = NA_real_,
                  mean_size = NA_real_, max_size = NA_integer_,
                  purity = NA_real_, signal_recovery = 0)
  if (is.null(clu) || !nrow(clu)) return(empty)
  clu <- distinct(clu, cdr3b, cluster_id)
  cs <- clu %>% count(cluster_id, name = "size") %>% filter(size >= min_size)
  if (!nrow(cs)) return(empty)
  big <- semi_join(clu, cs, by = "cluster_id")
  tibble(n_clusters      = nrow(cs),
         n_clustered     = n_distinct(big$cdr3b),
         median_size     = median(cs$size),
         mean_size       = round(mean(cs$size), 2),
         max_size        = max(cs$size),
         purity          = .purity_multi(big, lab, min_size),
         # signal recovery is reported alongside purity because purity alone
         # cannot distinguish "clean because selective" from "clean because
         # it barely clustered anything".
         signal_recovery = n_distinct(big$cdr3b[big$cdr3b %in% lab$cdr3b]) / n_labelled)
}



# -----------------------------------------------------------------
# A3b. Published-method parsers
# -----------------------------------------------------------------
# One parser per published method, used by run_benchmark_a() to read each
# method's cluster assignments.

jouannet_parsers <- list(
  HD = function(i, k, out) {
    d <- .j_read(file.path(out, "output_HD", sprintf("subset_%d%s.csv", i, k)), delim = "\t")
    if (is.null(d)) return(NULL)
    tibble(cdr3b = as.character(d$CDR3), cluster_id = as.character(d$clusters)) },
  LD = function(i, k, out) {
    d <- .j_read(file.path(out, "output_LD", sprintf("subset_%d%s.csv", i, k)), delim = "\t")
    if (is.null(d)) return(NULL)
    tibble(cdr3b = as.character(d$CDR3), cluster_id = as.character(d$clusters)) },
  iSMART = function(i, k, out) {
    d <- .j_read(file.path(out, "output_iSMART",
                 sprintf("subset_%d%s_results.txt_ClusteredCDR3s_7.5.txt", i, k)), delim = "\t")
    if (is.null(d) || ncol(d) < 2) return(NULL)
    tibble(cdr3b = as.character(d[[1]]), cluster_id = as.character(d[[2]])) },
  GIANA = function(i, k, out) {
    f <- file.path(out, "output_GIANA", sprintf("subset_%d%s_results.txt", i, k))
    if (!file.exists(f)) return(NULL)
    d <- .j_read(f, delim = "\t", col_names = FALSE, comment = "#")
    if (is.null(d) || ncol(d) < 2) return(NULL)
    tibble(cdr3b = as.character(d[[1]]), cluster_id = as.character(d[[2]])) %>%
      filter(!is.na(cdr3b), nchar(cdr3b) >= 6) },
  TCRMatch = function(i, k, out) {
    f <- file.path(out, "output_TCRMatch", sprintf("output_TCRMatch_sample%d%s.txt", i, k))
    if (!file.exists(f)) return(NULL)
    d <- .j_read(f, delim = .sniff_delim(f), col_names = c("from","to","score"),
                 comment = "#", trim_ws = TRUE)
    if (is.null(d) || ncol(d) < 3) return(NULL)
    d <- d %>% mutate(score = suppressWarnings(as.numeric(score))) %>%
      filter(!is.na(from), !is.na(to), !is.na(score), score != 1, from != to)
    if (!nrow(d)) return(NULL)
    nd <- unique(c(d$from, d$to)); ix <- setNames(seq_along(nd), nd)
    tibble(cdr3b = nd,
           cluster_id = as.character(connected_components(length(nd), ix[d$from], ix[d$to]))) },
  clusTCR = function(i, k, out) {
    d <- .j_read(file.path(out, "output_clusTCR", sprintf("subset%d%s", i, k),
                           sprintf("subset_%d_output.csv", i)), delim = ",")
    if (is.null(d) || ncol(d) < 2) return(NULL)
    tibble(cdr3b = as.character(d[[1]]), cluster_id = as.character(d[[2]])) },
  GLIPH2 = function(i, k, out) {
    d <- .j_read(file.path(out, "output_gliph2", sprintf("output_subset%d%s", i, k),
                           "output_cluster.csv"), delim = ",")
    if (is.null(d)) return(NULL)
    d <- d[, 1:min(18, ncol(d))]
    if ("TcRa" %in% names(d)) d <- select(d, -TcRa)
    d <- na.omit(d)
    if ("pattern" %in% names(d)) d <- filter(d, pattern != "single")
    if (!all(c("TcRb","index") %in% names(d))) return(NULL)
    tibble(cdr3b = as.character(d$TcRb), cluster_id = as.character(d$index)) %>% distinct() },
  DeepTCR = function(i, k, out) {
    d <- .j_read(file.path(out, "output_DeepTCR",
                 sprintf("results_clusteringTRB_subset%d%s.csv", i, k)), delim = ",")
    if (is.null(d)) return(NULL)
    cn <- names(d)[grepl("^cdr3$", names(d), ignore.case = TRUE)][1]
    cl <- names(d)[grepl("^cluster$", names(d), ignore.case = TRUE)][1]
    if (is.na(cn) || is.na(cl)) return(NULL)
    tibble(cdr3b = as.character(d[[cn]]), cluster_id = as.character(d[[cl]])) },
  TCRdist3 = function(i, k, out) {
    d <- .j_read(file.path(out, "output_TCRdist3", sprintf("subset%d%s", i, k),
                           "output_hclust_beta_cdr3.csv"), delim = ";")
    if (is.null(d) || !all(c("CDR3b","clusters") %in% names(d))) return(NULL)
    tibble(cdr3b = as.character(d$CDR3b), cluster_id = as.character(d$clusters)) }
)

DISEASE_METATYPE <- c(
  "Healthy"="Healthy",
  setNames(rep("Healthy",14), c(
    "Control (Boland et al.)", "Control (Buggert et al.)", "Control (Corridoni et al.)",
    "Control (Francis et al.)", "Control (Gao et al.)", "Control (He et al.)",
    "Control (Lee et al.)", "Control (Luoma et al.)", "Control (Notarbartolo et al., 2021)",
    "Control (Park et al.)", "Control (Ramaswamy et al.)", "Control (Suo et al.)",
    "Control (Wang et al., 2022)", "Control (Wen et al.)")),
  "Control (Suo et al.); Control (Park et al.)"="Healthy",
  "COVID-19"="COVID-19","COVID-19 vaccine"="COVID-19","MIS-C"="COVID-19",
  "CMV"="Infection","HIV"="Infection","Plasmodium falciparum infection"="Infection",
  "skin lesion erythema migrans"="Infection",
  setNames(rep("Solid tumor",12), c("Basal cell carinoma tumor","Breast cancer",
    "Clear cell renal cell carcinoma","Endometrial Carcinoma","Esophagus squamous cell carcinoma",
    "Gallbladder carcinoma","Gastric Cancer","HNSCC","Intrahepatic cholangiocarcinoma",
    "low-grade gliomas","Melanoma","Metastatic Cancer")),
  setNames(rep("Solid tumor",8), c("Metastatic colorectal cancer","Nasopharyngeal carcinoma",
    "NSCLC","Ovarian Cancer","Pancreatic Cancer","Prostate Cancer",
    "Squamous cell carcinoma tumor","Human thymic epithelial tumors (TET)")),
  "Large granular lymphocyte leukemia"="Leukemia/lymphoma",
  "Nodal B cell lymphoma"="Leukemia/lymphoma","Multiple Myeloma"="Leukemia/lymphoma",
  "Allergic rhinitis"="Inflammatory","Arthritis"="Inflammatory",
  "Chronic pancreatitis"="Inflammatory","CPI-colitis"="Inflammatory","ICB Arthritis"="Inflammatory",
  "Ulcerative colitis"="Inflammatory","Psoriatic arthritis"="Autoimmune",
  "Rheumatoid arthritis"="Autoimmune",
  "Aging"=NA_character_,"BEAM-T"=NA_character_,"Cord blood"=NA_character_,
  "Fraility"=NA_character_,"Frailty"=NA_character_,"Organoid"=NA_character_,
  "Young"=NA_character_
)

normalise_disease <- function(x) {
  v <- as.character(x)
  mapped <- DISEASE_METATYPE[v]
  unknown <- setdiff(unique(v[!is.na(v)&!v %in% names(DISEASE_METATYPE)]), NA)
  if (length(unknown)) {
    warning(length(unknown)," unrecognised disease labels mapped to 'Other': ",
            paste(head(unknown,5),collapse="; "), call.=FALSE)
    mapped[is.na(mapped) & !is.na(v) & !v %in% names(DISEASE_METATYPE)] <- "Other"
  }
  unname(mapped)
}

run_stage10_huardb <- function() {
  cat("\n=== Stage 10: HuARdb cohort-scale analysis ===\n")
  dir.create(CFG$huardb_out_dir, showWarnings=FALSE, recursive=TRUE)
  dist_dir <- file.path(CFG$huardb_root,"distance_tables")
  if (!dir.exists(dist_dir)) { cat("  No distance_tables/ dir at:", dist_dir, "\n"); return(invisible(NULL)) }
  dist_files <- list.files(dist_dir, "\\.tsv(\\.gz)?$", recursive=TRUE, full.names=TRUE)
  if (!length(dist_files)) { cat("  No distance files\n"); return(invisible(NULL)) }
  cat("  Distance files:", length(dist_files), "\n")
  thr <- CFG$huardb_threshold %||% CFG$focal_threshold
  labels <- NULL
  if (file.exists(CFG$huardb_input)) {
    inp <- tryCatch(read_tsv(CFG$huardb_input, show_col_types=FALSE), error=\(e) NULL)
    if (!is.null(inp)&&nrow(inp)) {
      labels <- mutate(inp, id=clean_id(.[[1]]))
      if (file.exists(CFG$huardb_metadata)) {
        md <- tryCatch(read_tsv(CFG$huardb_metadata,show_col_types=FALSE,guess_max=100000), error=\(e) NULL)
        if (!is.null(md)&&nrow(md)) {
          jk <- intersect(c("v_call_beta","j_call_beta","cdr3_aa_beta"), intersect(names(labels),names(md)))
          if (length(jk)) {
            md_sub <- md %>% select(-any_of(setdiff(intersect(names(md),names(labels)),jk))) %>%
              semi_join(labels, by=jk) %>% group_by(across(all_of(jk))) %>%
              summarise(across(everything(), \(x){v<-unique(x[!is.na(x)]);if(!length(v))NA_character_ else paste(sort(as.character(v)),collapse=";")}), .groups="drop")
            labels <- left_join(labels, md_sub, by=jk, relationship="many-to-one")
          }
        }
      }
    }
  }
  if (is.null(labels)||!nrow(labels)) {
    ids <- unique(unlist(map(head(dist_files,50), \(f) tryCatch({d<-read_tsv(f,show_col_types=FALSE,n_max=5000);unique(c(d$id_i,d$id_j))},error=\(e)character(0)))))
    labels <- tibble(id=clean_id(ids), epitope=NA_character_, disease=NA_character_,
                     hla=NA_character_, donor_id=NA_character_)
  }
  rnm <- function(df, tgt, cands) { if(tgt %in% names(df)) return(df); f<-intersect(cands,names(df)); if(length(f)) df[[tgt]]<-df[[f[1]]] else df[[tgt]]<-NA; df }
  labels <- labels %>%
    rnm("epitope",  c("Epitope","antigen","Antigen")) %>%
    rnm("hla",      c("HLA","hla_allele","HLA_allele","MHC")) %>%
    rnm("donor_id", c("donor","Donor","DonorID","sample_id","Sample_ID","patient_id"))
  src_col <- intersect(c("disease","disease_metatype"), names(labels))[1]
  if (!is.na(src_col)&&!is.null(src_col)) {
    labels <- mutate(labels, disease_detail=.data[[src_col]],
                     disease=vapply(strsplit(as.character(disease_detail),";",fixed=TRUE), function(parts) {
                       parts <- str_trim(parts[nchar(parts)>0]); if(!length(parts)) return(NA_character_)
                       m <- unique(normalise_disease(parts)[!is.na(normalise_disease(parts))])
                       if(!length(m)) NA_character_ else paste(sort(m),collapse=";")
                     }, character(1)))
  }
  clu <- build_sptce_clusters(dist_files, labels, rmsd_cut=thr, prefix="huar", keep_unlabelled=TRUE)
  cat("  Clusters:", n_distinct(clu$cluster_id), "| TCRs:", nrow(clu), "\n")
  write_tsv(clu, file.path(CFG$huardb_out_dir,"clusters.tsv"))
  sz <- clu %>% count(cluster_id, name="size") %>% filter(size>=2) %>% arrange(desc(size))
  write_tsv(sz, file.path(CFG$huardb_out_dir,"cluster_sizes.tsv"))
  cat("  Cluster size summary:\n"); print(summary(sz$size))
  p <- sz %>% count(size, name="n") %>%
    ggplot(aes(size,n)) + geom_col(fill=GROUP_PALETTE["SPACE-TCR RMSD"], width=0.85) +
    scale_x_log10(breaks=c(2,5,10,50,100,500,1000,5000)) +
    scale_y_log10(breaks=c(1,10,100,1000,10000)) +
    labs(title=sprintf("HuARdb cluster sizes (RMSD <= %s A)", thr), x="Cluster size", y="Number of clusters")
  save_plot("fig_huardb_cluster_sizes.png", p, 10, 6, dir=CFG$huardb_out_dir)
  if ("disease" %in% names(clu) && any(!is.na(clu$disease))) {
    enrich <- .huardb_disease_enrichment(clu, min_size=CFG$huardb_min_size)
    if (nrow(enrich)) {
      write_tsv(enrich, file.path(CFG$huardb_out_dir,"disease_enrichment.tsv"))
      p33 <- ggplot(enrich, aes(in_cluster/cluster_size, -log10(p_adj),
                                size=cluster_size, colour=-log10(p_adj)>1.3)) +
        geom_point(alpha=0.7) + geom_hline(yintercept=-log10(0.05),linetype="dashed",colour="grey40") +
        scale_colour_manual(values=c(`TRUE`=unname(GROUP_PALETTE["SPACE-TCR RMSD"]),`FALSE`="#bbbbbb"),guide="none") +
        scale_x_continuous(labels=percent_format()) +
        labs(title="Disease enrichment per cluster",subtitle="BH-adjusted Fisher; dashed = p < 0.05",
             x="Fraction with disease", y=expression(-log[10]~"(adj p)"))
      save_plot("fig_huardb_disease_volcano.png", p33, 10, 8, dir=CFG$huardb_out_dir)
    }
  }
  if ("hla" %in% names(clu) && any(!is.na(clu$hla))) {
    hla_df <- clu %>% filter(!is.na(hla)) %>%
      mutate(hla=strsplit(as.character(hla),";",fixed=TRUE)) %>%
      tidyr::unnest(hla) %>% mutate(hla=str_trim(hla)) %>% filter(nchar(hla)>0) %>%
      group_by(cluster_id) %>%
      summarise(n_hla=n_distinct(hla),
                dom_frac={t<-table(hla);if(length(t))max(t)/sum(t) else NA_real_}, .groups="drop") %>%
      left_join(count(clu,cluster_id,name="size"), by="cluster_id") %>%
      filter(size >= CFG$huardb_min_size)
    if (nrow(hla_df)) {
      write_tsv(hla_df, file.path(CFG$huardb_out_dir,"hla_coherence.tsv"))
      p35 <- ggplot(hla_df, aes(size, dom_frac)) +
        geom_point(alpha=0.45, size=1.1, colour=unname(GROUP_PALETTE["SPACE-TCR RMSD"])) +
        geom_smooth(method="loess", se=TRUE, span=0.9, colour="grey30", linewidth=0.55) +
        scale_x_log10() + scale_y_continuous(labels=percent_format(), limits=c(0,1.01)) +
        labs(title="HLA coherence vs cluster size",
             subtitle="Fraction of cluster sharing the dominant HLA allele",
             x="Cluster size (log10)", y="Dominant-HLA fraction")
      save_plot("fig_huardb_hla_coherence.png", p35, 10, 7, dir=CFG$huardb_out_dir)
    }
  }
  invisible(list(clusters=clu))
}

# -----------------------------------------------------------------
# A6. The comparison
# -----------------------------------------------------------------
# versions = "a" by default: a/b/c are replicate RUNS, and HD/LD/iSMART/GIANA/
# TCRMatch are deterministic (byte-identical files), so averaging them would
# fake precision. Pass c("a","b","c") to see run-to-run spread for the
# stochastic tools (clusTCR, DeepTCR).
#
# dih_cuts: dihedral cutoffs in radians, crossed with every RMSD threshold.
# Set NULL to force RMSD-only; ignored automatically when no dihedral table is
# found. Each SPACE-TCR row carries `param` (RMSD) and `dih` (NA for RMSD-only).
score_paired_all <- function(lab, ladder = NULL, ref_subset = NULL, versions = "a",
                             rmsd_thresholds = CFG$thresholds,
                             dih_cuts = CFG$dihedral_modes,
                             min_size = 2) {
  cat("\n=== Benchmark A: paired-only comparison ===\n")
  if (is.null(ladder)) ladder <- jouannet_10x_ladder(quiet = TRUE)
  if (is.null(ladder)) return(invisible(NULL))
  pr <- .paired_betas(ref_subset)
  if (is.null(pr)) { message("No structural run found"); return(invisible(NULL)) }

  dih_map <- .find_benchA_dihedrals(pr$ref_subset)
  if (is.null(dih_map)) dih_cuts <- NULL
  pairs <- .load_ref_pairs(pr$ref_subset, max(rmsd_thresholds), dih_map,
                           with_dih = !is.null(dih_map))
  if (!nrow(pairs)) { message("No pairs"); return(invisible(NULL)) }

  # SPACE-TCR operating points: RMSD alone, plus RMSD x dihedral if available.
  modes <- c(list(list(dih = NA_real_, tag = "rmsd")),
             if (length(dih_cuts)) map(dih_cuts, \(dc) list(dih = dc, tag = sprintf("d%s", dc))))
  cat("  SPACE-TCR modes: ", length(rmsd_thresholds), " RMSD thresholds x ",
      length(modes), " dihedral setting(s)\n", sep = "")

  rows <- list()
  for (i in ladder$table$subset) {
    uni   <- intersect(ladder$subsets[[i]], pr$betas)   # the shared universe
    n_lab <- sum(uni %in% lab$cdr3b)
    noise <- ladder$noise[[as.character(i)]]
    cat(sprintf("\n  subset %d (%4.1f%% noise): %d of %d betas paired (%.1f%%), %d labelled\n",
                i, noise, length(uni), length(ladder$subsets[[i]]),
                100 * length(uni) / length(ladder$subsets[[i]]), n_lab))

    # ---- published methods, filtered to the shared universe ----
    for (m in names(jouannet_parsers)) for (k in versions) {
      a <- tryCatch(jouannet_parsers[[m]](i, k, CFG$tenx_out), error = \(e) NULL)
      if (is.null(a) || !nrow(a)) next
      before <- n_distinct(a$cdr3b)
      a2 <- filter(a, cdr3b %in% uni)
      if (!nrow(a2)) next
      s <- .paired_score(a2, lab, n_lab, min_size)
      rows[[length(rows) + 1]] <- s %>%
        mutate(method = m, param = NA_real_, dih = NA_real_,
               subset = i, pct_noise = noise,
               n_universe = length(uni), n_labelled = n_lab,
               kept_frac = round(n_distinct(a2$cdr3b) / before, 3), .before = 1)
      cat(sprintf("    %-9s purity %.3f  mean size %5.2f  clusters %5d  (kept %.0f%% of its betas)\n",
                  m, s$purity, s$mean_size, s$n_clusters,
                  100 * n_distinct(a2$cdr3b) / before))
    }

    # ---- SPACE-TCR, paired by construction ----
    ids <- pr$map$id[pr$map$cdr3b %in% uni]
    lab_stub <- tibble(id = pr$map$id, epitope = NA_character_)
    for (thr in rmsd_thresholds) for (mc in modes) {
      cl <- .cluster_from_pairs(
        pairs, ids, rmsd_cut = thr,
        dih_cut = if (is.na(mc$dih)) NULL else mc$dih,
        prefix  = sprintf("pp_s%d_%s_r%s", i, mc$tag, thr),
        labels  = lab_stub, keep_unlabelled = TRUE)
      if (is.null(cl) || !nrow(cl)) next
      b <- cl %>% left_join(pr$map, by = c("tcr_id" = "id")) %>%
        filter(!is.na(cdr3b)) %>% distinct(cdr3b, cluster_id)
      s <- .paired_score(b, lab, n_lab, min_size)
      nm <- if (is.na(mc$dih)) "SPACE-TCR" else "SPACE-TCR+dih"
      rows[[length(rows) + 1]] <- s %>%
        mutate(method = nm, param = thr, dih = mc$dih,
               subset = i, pct_noise = noise,
               n_universe = length(uni), n_labelled = n_lab,
               kept_frac = 1, .before = 1)
      cat(sprintf("    %-13s purity %.3f  mean size %5.2f  clusters %5d  (RMSD %s A%s)\n",
                  nm, s$purity, s$mean_size, s$n_clusters, thr,
                  if (is.na(mc$dih)) "" else sprintf(", dih<%s rad", mc$dih)))
    }
  }
  if (!length(rows)) { message("nothing scored"); return(invisible(NULL)) }
  res <- bind_rows(rows)
  write_csv(res, file.path(CFG$analysis_dir, "paired_only_comparison.csv"))
  cat("\n  wrote paired_only_comparison.csv\n")
  invisible(res)
}


# -----------------------------------------------------------------
# A7. Figures
# -----------------------------------------------------------------
# Mean, not median. Median cluster size is exactly 2 in ~93% of method/noise
# cells, so every visible movement is one cluster tipping 2 -> 3 on a
# compressed axis. Mean separates the methods and shows the effect that
# matters: LD's mean cluster size roughly doubles as noise rises (it is
# absorbing junk) while SPACE-TCR at tight thresholds stays flat.
.size_col <- function(res, size_metric) {
  size_metric <- match.arg(size_metric, c("mean", "median"))
  if (size_metric == "median") {
    d <- filter(res, !is.na(median_size))
    if (nrow(d) && mean(d$median_size == 2) > 0.7)
      warning("Median cluster size is 2 in >70% of cells and cannot separate ",
              "the methods. Use size_metric = 'mean'.", call. = FALSE)
    list(col = "median_size", lab = "Median cluster size")
  } else list(col = "mean_size", lab = "Mean cluster size")
}


# Label for a SPACE-TCR operating point.
.sptce_tag <- function(param, dih)
  ifelse(is.na(dih), sprintf("%s Å", param), sprintf("%s Å + dih<%s", param, dih))


# -----------------------------------------------------------------
# A8. ENTRYPOINT
# -----------------------------------------------------------------
run_benchmark_a <- function(versions = "a", rmsd_fig = c(0.5, 0.75, 1.0),
                            dih_fig = NA, size_metric = "mean",
                            min_clusters = 5L, dih_cuts = CFG$dihedral_modes) {
  lad <- jouannet_10x_ladder()
  if (is.null(lad)) return(invisible(NULL))
  lab <- jouannet_10x_labels()
  if (is.null(lab)) return(invisible(NULL))
  res <- score_paired_all(lab, ladder = lad, versions = versions, dih_cuts = dih_cuts)
  if (is.null(res)) return(invisible(NULL))
  fig_paired_methods(res, rmsd = rmsd_fig, dih_keep = dih_fig,
                     min_clusters = min_clusters, size_metric = size_metric)
  fig_paired_thresholds(res, min_clusters = min_clusters, size_metric = size_metric)
  cat("\nBenchmark A complete. Outputs in: ", CFG$analysis_dir, "\n", sep = "")
  invisible(res)
}


# #################################################################
##
##   BENCHMARK B — Jouannet curated-cohort benchmark (stages 5-7)
##
##   run_benchmark_b()
##
# #################################################################

# =================================================================
# 11. STAGE 5 — Jouannet-cohort benchmark
# =================================================================
run_stage5 <- function() {
  cat("\n=== Stage 5: Jouannet-cohort benchmark ===\n")
  if (!file.exists(CFG$jouannet_labelled)) {
    message("Labelled input missing: ", CFG$jouannet_labelled); return(invisible(NULL))
  }
  labels    <- read_tsv(CFG$jouannet_labelled, show_col_types=FALSE) %>%
    select(id, cdr3_aa_beta, epitope, any_of("cdr3_aa_alpha"))
  total_n   <- nrow(labels)
  # Eq. 4's denominator is the INPUT count per epitope. The input table is
  # available here, so the counts are set once for sens_for_df() downstream;
  # without them every sensitivity is NA.
  .set_sens_input_counts(labels)
  cat("  Sensitivity denominators set for", length(CFG$.sens_input_counts),
      "epitopes from the input table\n")
  dist_dir  <- file.path(CFG$jouannet_sptce_out, "distance_tables")
  dist_files <- list.files(dist_dir, "\\.tsv\\.gz$", full.names=TRUE)
  cat("  Input TCRs:", total_n, " | Distance files:", length(dist_files), "\n")
  if (!length(dist_files)) { message("No distance files"); return(invisible(NULL)) }

  cat("  Building SPACE-TCR RMSD clusters...\n")
  sptce_rmsd <- map_dfr(CFG$thresholds, function(thr) {
    cat("   ", thr, "A\n")
    build_sptce_clusters(dist_files, labels, rmsd_cut=thr, prefix=paste0("r",thr)) %>%
      mutate(threshold=thr, mode="SPACE-TCR RMSD")
  })

  dih_map <- NULL; sptce_dih <- NULL; sptce_dih_all_modes <- NULL
  dih_file <- file.path(CFG$jouannet_sptce_out, "dihedrals.tsv.gz")
  if (file.exists(dih_file)) {
    cat("  Loading dihedrals...\n")
    dih_map <- load_dihedral_map(dih_file)
    sptce_dih <- map_dfr(CFG$thresholds, function(thr) {
      build_sptce_clusters(dist_files, labels, rmsd_cut=thr,
                           dih_cut=CFG$focal_dihedral, dih_map=dih_map,
                           prefix=paste0("d",CFG$focal_dihedral,"_",thr)) %>%
        mutate(threshold=thr, mode="SPACE-TCR RMSD+dih")
    })
    sptce_dih_all_modes <- map_dfr(CFG$dihedral_modes, function(dc) {
      map_dfr(CFG$thresholds, function(thr) {
        build_sptce_clusters(dist_files, labels, rmsd_cut=thr,
                             dih_cut=dc, dih_map=dih_map,
                             prefix=paste0("dm",dc,"_",thr)) %>%
          mutate(threshold=thr, mode=paste0("SPACE-TCR RMSD+dih<",dc,"rad"))
      })
    })
  }
  all_sptce <- bind_rows(sptce_rmsd, sptce_dih)
  if (!"cdr3_aa_alpha" %in% names(all_sptce)) all_sptce$cdr3_aa_alpha <- NA_character_

  sptce_metrics <- all_sptce %>% group_by(threshold, mode) %>%
    group_modify(~ compute_metrics(.x, total_n) %>%
                   mutate(method=paste0(.y$mode," (",.y$threshold,"A)"))) %>%
    ungroup()
  write_csv(sptce_metrics, file.path(CFG$jouannet_out_dir,"sptce_all_threshold_metrics.csv"))

  cat("  Parsing published method outputs...\n")
  db_path    <- file.path(CFG$jouannet_repo,"Data/Database/database_pooled_human_2023_03_15.txt")
  method_dir <- file.path(CFG$jouannet_repo,"Data/Output_methods")

  # Compute the published-method metrics in-script with the paper-exact engine
  # (paired units + multi-epitope expansion + purity>=0.9) and use them for all
  # downstream figures, instead of the hardcoded fallback table. SPACE-TCR is
  # scored by compute_metrics() with the identical purity>=0.9 rule, so the two
  # sides of every comparison come from the same equations.
  PUBLISHED_METRICS <<- tryCatch(
    compute_published_metrics_exact(db_path, method_dir),
    error = function(e) {
      warning("compute_published_metrics_exact failed (", conditionMessage(e),
              "); falling back to hardcoded table")
      PUBLISHED_CLUSTERS_EXACT <<- NULL
      PUBLISHED_METRICS_FALLBACK
    })

  db_join    <- load_jouannet_curated(db_path, cell_type=CFG$jouannet_cell_type,
                                      verbose=FALSE) %>% select(cdr3b, epitope)
  finalise <- function(df) {
    df %>% mutate(across(c(cluster_id,cdr3b), as.character)) %>%
      filter(nchar(cdr3b)>=6) %>%
      inner_join(db_join, by="cdr3b") %>%
      group_by(cluster_id) %>% filter(n()>=2) %>% ungroup() %>%
      select(method, cluster_id, cdr3b, epitope)
  }
  members <- list()
  read_method <- function(name, path, parse_fn) {
    if (!file.exists(path)) return(NULL)
    tryCatch(parse_fn(path), error=\(e) { warning(name,": ",conditionMessage(e)); NULL })
  }
  members$clusTCR <- read_method("clusTCR",
                                 file.path(method_dir,"clusTCR_output.csv"),
                                 \(f) read_csv(f, show_col_types=FALSE) %>%
                                   mutate(cdr3b=vapply(junction_aa, \(x){m<-regmatches(x,regexpr("^C[A-Z]+?F",x));if(!length(m))NA_character_ else m},character(1))) %>%
                                   filter(!is.na(cdr3b)) %>% rename(cluster_id=cluster) %>% mutate(method="clusTCR") %>% finalise())
  members$GIANA <- read_method("GIANA",
                               file.path(method_dir,"GIANA_output.txt"),
                               \(f) readLines(f, warn=FALSE) %>% .[!str_starts(.,"##")] %>%
                                 str_split_fixed("\t",2) %>% as_tibble(.name_repair="minimal") %>%
                                 setNames(c("cdr3b","cluster_id")) %>%
                                 filter(nchar(cdr3b)>=6, str_detect(cdr3b,"^C[A-Z]+F$")) %>%
                                 mutate(method="GIANA") %>% finalise())
  members$iSMART <- read_method("iSMART",
                                file.path(method_dir,"iSMART_output.txt"),
                                \(f) { d<-read_delim(f,delim="\t",show_col_types=FALSE)
                                rename(d, cdr3b=names(d)[str_detect(names(d),"CDR3")][1],
                                       cluster_id=names(d)[str_detect(tolower(names(d)),"group|cluster")][1]) %>%
                                  mutate(method="iSMART") %>% finalise() })
  members$GLIPH2 <- read_method("GLIPH2",
                                file.path(method_dir,"GLIPH2_output_cluster.csv"),
                                \(f) read_csv(f,show_col_types=FALSE) %>%
                                  rename(cdr3b=TcRb, cluster_id=index) %>% mutate(method="GLIPH2") %>% finalise())
  members$TCRMatch <- read_method("TCRMatch",
                                  file.path(method_dir,"TCRMatch_output.txt"),
                                  \(f) { raw<-readLines(f,warn=FALSE); raw<-raw[!str_starts(raw,"#")&nchar(raw)>0]
                                  dl<-raw[1]; delim<-if(str_count(dl,"\t")>=1)"\t" else if(str_count(dl,",")>=2)"," else " "
                                  mat<-str_split_fixed(raw,delim,3)
                                  if(str_detect(mat[1,1],regex("query|cdr3",ignore_case=TRUE))) mat<-mat[-1,,drop=FALSE]
                                  edges<-tibble(from=mat[,1],to=mat[,2]) %>%
                                    filter(nchar(from)>=6,nchar(to)>=6,from!=to,str_detect(from,"^C"),str_detect(to,"^C"))
                                  if(!nrow(edges)) return(NULL)
                                  nodes<-unique(c(edges$from,edges$to)); idx<-setNames(seq_along(nodes),nodes)
                                  comp<-connected_components(length(nodes),idx[edges$from],idx[edges$to])
                                  tibble(cdr3b=nodes,cluster_id=as.character(comp)) %>% mutate(method="TCRMatch") %>% finalise()})
  members$DeepTCR <- read_method("DeepTCR",
                                 file.path(method_dir,"DeepTCR_output.csv"),
                                 \(f) read_csv(f,show_col_types=FALSE) %>%
                                   rename(cdr3b=Beta_Sequences, cluster_id=cluster) %>%
                                   mutate(method="DeepTCR") %>% finalise())
  cat("  Building LD<=1 clusters (slow)...\n")
  ld_raw <- db_join %>% filter(nchar(cdr3b)>=6) %>% mutate(len=nchar(cdr3b)) %>%
    group_by(len) %>% group_modify(function(g, k) {
      seqs<-g$cdr3b; n<-length(seqs); if(n<2) return(tibble())
      dm<-stringdistmatrix(seqs,seqs,method="lv")
      adj<-which(dm<=1&dm>0,arr.ind=TRUE); if(!nrow(adj)) return(tibble())
      g$cluster_id<-paste0("len",k$len,"_",connected_components(n,adj[,1],adj[,2])); g
    }) %>% ungroup() %>%
    group_by(cluster_id) %>% filter(n()>=2) %>% ungroup()
  members$Levenshtein_Distance <- mutate(ld_raw, method="Levenshtein_Distance") %>%
    select(method,cluster_id,cdr3b,epitope)
  members$Hamming_Distance <- mutate(ld_raw, method="Hamming_Distance",
                                     cluster_id=paste0("hd_",cluster_id)) %>%
    select(method,cluster_id,cdr3b,epitope)
  all_jouannet <- bind_rows(compact(members))
  # Prefer the paper-exact membership (paired methods keep pairs; single-chain
  # keep pooled alpha+beta) so the network / LD / sensitivity / polyspecificity
  # figures match the bubble and metric table. Falls back to beta-only if the
  # in-script engine did not run.
  pub_members <- PUBLISHED_CLUSTERS_EXACT %||% all_jouannet
  if (!is.null(PUBLISHED_CLUSTERS_EXACT))
    cat("  Published figures use paper-exact membership (", nrow(pub_members),
        " units across ", n_distinct(pub_members$method), " methods)\n", sep="")

  # 5.6 Sensitivity — uses the GLOBAL sens_for_df (defined in section 5)
  sptce_sens <- all_sptce %>% filter_sptce_focal() %>%
    group_by(mode, threshold) %>%
    group_modify(\(df, key) sens_for_df(df) %>%
                   mutate(method = paste0(key$mode, " (", key$threshold, "A)"))) %>%
    ungroup() %>% select(-mode, -threshold)
  pub_sens <- pub_members %>%
    group_by(method) %>%
    group_modify(\(df, key) sens_for_df(df)) %>%
    ungroup()
  sensitivity_all <- bind_rows(sptce_sens, pub_sens)
  write_csv(sensitivity_all, file.path(CFG$jouannet_out_dir,"sensitivity_all_methods.csv"))

  pub_ld <- summarise_clusters(pub_members, "cdr3b") %>%
    mutate(threshold=NA_real_, mode="published")
  sptce_ld <- all_sptce %>%
    mutate(method=paste0(mode," (",threshold,"A)")) %>%
    summarise_clusters("cdr3_aa_beta")
  all_ld <- bind_rows(pub_ld, sptce_ld)
  write_csv(all_ld, file.path(CFG$jouannet_out_dir,"ld_all_clusters.csv"))

  ld_summary <- all_ld %>% filter(!is.na(median_ld)) %>%
    group_by(method, threshold, mode) %>%
    summarise(n_clusters=n(), median_ld=median(median_ld,na.rm=TRUE),
              pct_le1=mean(median_ld<=1,na.rm=TRUE), pct_2_3=mean(median_ld>1&median_ld<=3,na.rm=TRUE),
              pct_4_6=mean(median_ld>3&median_ld<=6,na.rm=TRUE), pct_gt6=mean(median_ld>6,na.rm=TRUE),
              .groups="drop") %>% arrange(median_ld)
  write_csv(ld_summary, file.path(CFG$jouannet_out_dir,"ld_summary_all_methods.csv"))

  published_clusters <- select(pub_members, method, cluster_id, cdr3b, epitope)
  sptce_all_thresholds <- all_sptce %>%
    mutate(method=case_when(mode=="SPACE-TCR RMSD"~name_sptce_rmsd(threshold),
                            mode=="SPACE-TCR RMSD+dih"~name_sptce_focal_dih(threshold),
                            TRUE~paste0(mode," (",threshold,"A)")),
           cdr3b=cdr3_aa_beta, cdr3a=cdr3_aa_alpha) %>%
    select(method, threshold, mode, cluster_id, cdr3b, cdr3a, epitope)
  sptce_focal_clusters <- filter_sptce_focal(sptce_all_thresholds) %>%
    select(method, cluster_id, cdr3b, cdr3a, epitope)

  # ---- focal SPACE-TCR rows for Table 4.2 (RMSD@0.75 + RMSD+dih@1.0) ----
  sptce_focal_poly <- sptce_focal_clusters %>%
    group_by(method, cluster_id) %>%
    summarise(is_poly = n_distinct(epitope) > 1, .groups = "drop_last") %>%
    summarise(perc_polyspecific = mean(is_poly), .groups = "drop")
  sptce_focal_table <- sptce_metrics %>%
    filter_sptce_focal() %>%
    transmute(method, n_clusters, n_clustered, perc_clusters_gt09) %>%
    left_join(filter_sptce_focal(ld_summary) %>% transmute(method, median_ld),
              by = "method") %>%
    left_join(sptce_focal_poly, by = "method")
  write_csv(sptce_focal_table, file.path(CFG$jouannet_out_dir, "sptce_focal_table.csv"))
  cat("  Focal SPACE-TCR table rows (Table 4.2):\n")
  print(as.data.frame(sptce_focal_table))

  write_csv(PUBLISHED_METRICS,
            file.path(CFG$jouannet_out_dir, "published_method_metrics.csv"))

  # ---- combined Table 4.2 (SPACE-TCR focal + published) in ONE csv ----
  pub_counts <- pub_members %>% group_by(method) %>%
    summarise(n_clusters = n_distinct(cluster_id), n_clustered = n(), .groups = "drop")
  pub_poly <- pub_members %>% group_by(method, cluster_id) %>%
    summarise(np = n_distinct(epitope) > 1, .groups = "drop_last") %>%
    summarise(perc_polyspecific = mean(np), .groups = "drop")
  pub_ld_med <- ld_summary %>% filter(mode == "published") %>% select(method, median_ld)
  published_table <- PUBLISHED_METRICS %>% select(method, perc_clusters_gt09) %>%
    left_join(pub_counts, by = "method") %>%
    left_join(pub_ld_med, by = "method") %>%
    left_join(pub_poly,   by = "method") %>%
    transmute(block = "published", method, n_clusters, n_clustered,
              perc_clusters_gt09, median_ld, perc_polyspecific)
  benchmark_table <- bind_rows(
    sptce_focal_table %>% mutate(block = "SPACE-TCR") %>%
      select(block, method, n_clusters, n_clustered,
             perc_clusters_gt09, median_ld, perc_polyspecific),
    published_table)
  write_csv(benchmark_table,
            file.path(CFG$jouannet_out_dir, "benchmark_comparison_table.csv"))
  cat("  Combined benchmark table (Table 4.2):\n")
  print(as.data.frame(benchmark_table))

  # ---- CDR3b sequence-divergence classes of high-purity SPACE-TCR clusters ----
  conv <- all_ld %>%
    filter(mode == "SPACE-TCR RMSD",
           abs(threshold - CFG$focal_threshold) < 1e-6,
           n_members >= 2, !is.na(median_ld), epitope_purity >= 0.9) %>%
    mutate(class = as.character(ld_class_of(median_ld))) %>%
    count(class, name = "n_clusters") %>%
    mutate(class = factor(class, levels = LD_CLASSES)) %>%
    tidyr::complete(class, fill = list(n_clusters = 0)) %>%
    arrange(class) %>%
    mutate(share = n_clusters / sum(n_clusters))
  write_csv(conv, file.path(CFG$jouannet_out_dir, "sptce_convergence_classes.csv"))
  cat("  Convergence classes (high-purity SPACE-TCR @ focal):\n")
  print(as.data.frame(conv))

  cat("  Generating figures...\n")
  .s5_figures(sptce_metrics, sptce_focal_clusters, sptce_all_thresholds,
              published_clusters, all_sptce, all_ld, sensitivity_all,
              total_n, dih_map)
  .s5_figures_extra(sptce_metrics, sptce_focal_clusters, sptce_all_thresholds,
                    published_clusters, all_sptce, all_ld, sensitivity_all,
                    total_n, dih_map, sptce_dih, sptce_dih_all_modes, members)

  cat("Stage 5 complete.\n")
  invisible(list(
    metrics              = .s5_bubble_data(sptce_focal_clusters, published_clusters,
                                           all_sptce, total_n, 2),
    sensitivity          = sensitivity_all,
    sptce_clusters       = sptce_focal_clusters,
    sptce_all_thresholds = sptce_all_thresholds,
    published_clusters   = published_clusters,
    all_sptce            = all_sptce,
    all_ld               = all_ld,
    sptce_metrics        = sptce_metrics,
    sptce_dih            = sptce_dih,
    sptce_dih_all_modes  = sptce_dih_all_modes,
    members              = members,
    dih_map              = dih_map,
    total_n              = total_n
  ))
}


.s5_bubble_data <- function(sptce_focal, pub_clu, all_sptce, total_n, min_size) {
  pub <- PUBLISHED_METRICS %>% mutate(n_clusters=9999L, colour_group="published")
  sptce_at <- all_sptce %>%
    filter_sptce_focal() %>%
    group_by(mode, threshold) %>%
    group_modify(~compute_metrics(.x,total_n,min_size=min_size) %>%
                   mutate(method=paste0(.y$mode," (",.y$threshold,"A)"))) %>% ungroup() %>%
    mutate(colour_group=case_when(mode=="SPACE-TCR RMSD"~"SPACE-TCR RMSD",
                                  mode=="SPACE-TCR RMSD+dih"~"SPACE-TCR RMSD+dih (focal)")) %>%
    select(method,retention,purity,perc_clusters_gt09,perc_tcrs_gt09,n_clusters,colour_group)
  bind_rows(pub, sptce_at)
}


# ── .s5_figures: figs 01,02,03,05,10,11,14,17,18 ───────────────
.s5_figures <- function(sptce_metrics, sptce_focal, sptce_all, pub_clu,
                        all_sptce, all_ld, sensitivity_all, total_n, dih_map) {
  out <- CFG$jouannet_out_dir
  p1 <- ggplot() +
    geom_point(data=PUBLISHED_METRICS, aes(retention,purity), colour="grey55", size=2.5) +
    geom_text(data=PUBLISHED_METRICS, aes(retention,purity,label=method), colour="grey40", size=2.4, vjust=-0.9) +
    geom_line(data=filter(sptce_metrics, !is.na(mode)), aes(retention,purity,colour=mode,linetype=mode), linewidth=1.1) +
    geom_point(data=filter(sptce_metrics, !is.na(mode)), aes(retention,purity,colour=mode,shape=mode), size=3.5) +
    geom_text(data=sptce_metrics, aes(retention,purity,label=threshold,colour=mode), size=2.4, vjust=2.5, show.legend=FALSE) +
    scale_colour_manual(values=PALETTES$mode_colour) + scale_linetype_manual(values=PALETTES$mode_line) +
    scale_shape_manual(values=PALETTES$mode_shape) +
    scale_x_continuous(labels=percent_format(), limits=c(0,1)) +
    scale_y_continuous(labels=percent_format(), limits=c(0.35,1.02)) +
    labs(title="SPACE-TCR purity-retention tradeoff", x="Retention", y="Global cluster purity",
         colour=NULL, linetype=NULL, shape=NULL)
  save_plot("fig01_threshold_sweep_arc.png", p1, 10, 7)

  sptce_long <- sptce_metrics %>%
    pivot_longer(names(METRIC_LABELS), names_to="metric", values_to="value") %>%
    mutate(metric=factor(metric, names(METRIC_LABELS), METRIC_LABELS))
  pub_long <- PUBLISHED_METRICS %>%
    pivot_longer(names(METRIC_LABELS), names_to="metric", values_to="value") %>%
    mutate(metric=factor(metric, names(METRIC_LABELS), METRIC_LABELS))
  p2 <- ggplot() +
    geom_hline(data=pub_long, aes(yintercept=value, group=method), colour="grey83", linewidth=0.35) +
    geom_line(data=sptce_long, aes(threshold,value,colour=mode,linetype=mode), linewidth=1.1) +
    geom_point(data=sptce_long, aes(threshold,value,colour=mode,shape=mode), size=2.8) +
    geom_vline(xintercept=CFG$focal_threshold, linetype="dotted", colour="grey40", linewidth=0.4) +
    facet_wrap(~metric, nrow=2, scales="free_y") +
    scale_colour_manual(values=PALETTES$mode_colour) + scale_linetype_manual(values=PALETTES$mode_line) +
    scale_shape_manual(values=PALETTES$mode_shape) + scale_x_continuous(breaks=CFG$thresholds) +
    scale_y_continuous(labels=percent_format()) +
    labs(title="SPACE-TCR metrics across RMSD thresholds", x="RMSD threshold (A)", y=NULL, colour=NULL, linetype=NULL, shape=NULL)
  save_plot("fig02_metric_panel.png", p2, 11, 8)

  p3 <- sensitivity_all %>%
    mutate(method=fct_rev(order_methods(method)),
           epitope=factor(epitope, CFG$jouannet_epitopes),
           label=ifelse(is.na(sensitivity),"",sprintf("%.2f",sensitivity))) %>%
    ggplot(aes(epitope,method,fill=sensitivity)) +
    geom_tile(colour="white", linewidth=0.6) +
    geom_text(aes(label=label, colour=sensitivity>0.45), size=4.4) +
    scale_fill_gradient2(low="#F7F7F7", mid="#FEE08B", high="#D73027",
                         midpoint=0.35, limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1),
                         na.value="grey93", name="Sensitivity",
                         guide=guide_colourbar(barwidth=18, barheight=0.8,
                                               title.position="top", title.hjust=0.5)) +
    scale_colour_manual(values=c(`FALSE`="grey20",`TRUE`="white"), guide="none") +
    scale_x_discrete(guide=guide_axis(angle=30)) +
    labs(title="Sensitivity: all methods x 6 epitopes", x="Epitope", y=NULL) +
    theme(panel.grid=element_blank(), legend.position="bottom",
          plot.margin=margin(6, 14, 6, 6))
  save_plot("fig03_sensitivity_heatmap.png", p3, 11, 8)

  focal_ld <- keep_focal_rows(all_ld) %>% filter(!is.na(median_ld)) %>%
    mutate(fill_group=case_when(
      mode=="published"     ~ "published",
      mode=="SPACE-TCR RMSD"    ~ "SPACE-TCR RMSD",
      mode=="SPACE-TCR RMSD+dih"~ "SPACE-TCR RMSD+dih (focal)",
      TRUE                  ~ "SPACE-TCR extra"),
      method=order_methods(method))
  p5 <- ggplot(focal_ld, aes(median_ld, fct_rev(method), fill=fill_group)) +
    geom_violin(alpha=0.6, scale="width", linewidth=0.3, draw_quantiles=c(0.25,0.5,0.75)) +
    geom_boxplot(width=0.12, fill="white", alpha=0.85, outlier.size=0.4, linewidth=0.35) +
    geom_vline(xintercept=c(1,3,6), linetype="dashed", colour="grey30", linewidth=0.45) +
    scale_fill_manual(values=GROUP_PALETTE, name=NULL) +
    scale_x_continuous(breaks=seq(0,16,2), limits=c(0,NA)) +
    labs(title="Within-cluster CDR3b LD", x="Median within-cluster CDR3b Levenshtein distance", y=NULL) +
    theme(axis.text.y=element_text(size=12))
  save_plot("fig05_ld_violin.png", p5, 13, 12)

  bubble_d <- .s5_bubble_data(sptce_focal, pub_clu, all_sptce, total_n, 2)
  p10 <- bubble_d %>%
    pivot_longer(c(retention,purity,perc_clusters_gt09,perc_tcrs_gt09),
                 names_to="metric", values_to="value") %>%
    mutate(metric=factor(metric, names(METRIC_LABELS), METRIC_LABELS),
           method=order_methods(method)) %>%
    ggplot(aes(metric, fct_rev(method), size=value, fill=value)) +
    geom_point(shape=21, colour="grey30", stroke=0.3, alpha=0.95) +
    geom_text(aes(label=sprintf("%.2f",value), colour=value>0.6),
              fontface="bold", size=4.0, show.legend=FALSE) +
    scale_colour_manual(values=c(`TRUE`="white", `FALSE`="grey15"), guide="none") +
    scale_size_area(max_size=16, guide="none") +
    scale_fill_gradient2(low="#F7F7F7", mid="#FEE08B", high="#D73027", midpoint=0.5, limits=c(0,1),
                         name="Value",
                         guide=guide_colourbar(barwidth=12, barheight=0.6,
                                               title.position="left", title.vjust=0.9)) +
    labs(title="Performance metrics (cluster size >= 2)", x=NULL, y=NULL) +
    theme(axis.text.x=element_text(size=14, face="bold"),
          axis.text.y=element_text(size=13),
          plot.title=element_text(size=17, face="bold"),
          panel.grid.major.y=element_line(colour="grey92", linewidth=0.3),
          panel.grid.major.x=element_blank(),
          legend.position="bottom")
  save_plot("fig10_bubble_metrics.png", p10, 11.5, 11)

  size_data <- bind_rows(
    filter(all_sptce, threshold==CFG$focal_threshold, mode=="SPACE-TCR RMSD") %>%
      group_by(cluster_id) %>% summarise(size=n(),.groups="drop") %>%
      mutate(method=name_sptce_rmsd(CFG$focal_threshold), colour_group="SPACE-TCR RMSD"),
    pub_clu %>% group_by(method,cluster_id) %>% summarise(size=n(),.groups="drop") %>%
      mutate(colour_group="published")
  ) %>% filter(size>=2) %>%
    mutate(size_bin=factor(case_when(size==2~"2",size==3~"3",size<=5~"4-5",
                                     size<=10~"6-10",size<=20~"11-20",TRUE~">20"),
                           levels=c("2","3","4-5","6-10","11-20",">20")))
  write_csv(size_data %>% group_by(method,colour_group) %>%
              summarise(n_clusters=n(),median_size=median(size),mean_size=round(mean(size),1),.groups="drop"),
            file.path(out,"cluster_size_summary.csv"))
  p11 <- size_data %>% count(method,colour_group,size_bin) %>%
    group_by(method) %>% mutate(pct=n/sum(n)) %>% ungroup() %>%
    mutate(method=order_methods(method)) %>%
    ggplot(aes(fct_rev(method),pct,fill=size_bin)) +
    geom_col(width=0.75, alpha=0.95, colour="white", linewidth=0.15) +
    coord_flip() +
    scale_fill_manual(values=PALETTES$size_bin, name="Cluster size",
                      guide=guide_legend(reverse=TRUE)) +
    scale_y_continuous(labels=percent_format(), expand=expansion(mult=c(0,0.01))) +
    labs(title="Cluster size distribution", x=NULL, y="Fraction of clusters") +
    theme(axis.text.y=element_text(size=13), plot.margin=margin(6, 14, 6, 6))
  save_plot("fig11_cluster_size_stacked.png", p11, 13, 9)

  poly_data <- bind_rows(
    filter(all_sptce, threshold==CFG$focal_threshold, mode=="SPACE-TCR RMSD") %>%
      select(cluster_id,epitope) %>%
      mutate(method=name_sptce_rmsd(CFG$focal_threshold), colour_group="SPACE-TCR RMSD"),
    pub_clu %>% select(method,cluster_id,epitope) %>% mutate(colour_group="published")
  ) %>% filter(!is.na(epitope))
  poly_stats <- poly_data %>%
    group_by(method,colour_group,cluster_id) %>%
    summarise(n_members=n(), n_epi=n_distinct(epitope), .groups="drop") %>%
    filter(n_members>=2) %>%
    mutate(spec_class=factor(case_when(
      n_epi==1~"Monospecific (1 epitope)", n_epi==2~"Bispecific (2 epitopes)",
      n_epi<=5~"Polyspecific (3-5 epitopes)", TRUE~"Highly polyspecific (>5 epitopes)"),
      levels=names(PALETTES$spec_class)))
  write_csv(poly_stats %>% group_by(method,colour_group) %>%
              summarise(n_clusters=n(), pct_monospecific=mean(n_epi==1),
                        pct_polyspecific=mean(n_epi>=2), .groups="drop"),
            file.path(out,"polyspecific_summary.csv"))
  p14 <- poly_stats %>% count(method,colour_group,spec_class) %>%
    group_by(method) %>% mutate(pct=n/sum(n)) %>% ungroup() %>%
    mutate(method=order_methods(method)) %>%
    ggplot(aes(fct_rev(method),pct,fill=spec_class)) +
    geom_col(width=0.75, alpha=0.92) + coord_flip() +
    scale_fill_manual(values=PALETTES$spec_class, name="Cluster type",
                      guide=guide_legend(nrow=2, byrow=TRUE)) +
    scale_y_continuous(labels=percent_format(), expand=expansion(mult=c(0,0.01))) +
    labs(title="Polyspecific clusters", x=NULL, y="Fraction of clusters") +
    theme(axis.text.y=element_text(size=13), legend.text=element_text(size=12),
          plot.margin=margin(6, 14, 6, 6))
  save_plot("fig14_polyspecific_stacked.png", p14, 14, 9)

  size_thrs <- c(2,3,4,5,6,8,10,15,20)
  purity_at_size <- function(df) {
    cs <- df %>% group_by(cluster_id) %>%
      summarise(size=n(), dom={t<-table(epitope[!is.na(epitope)]);if(length(t))max(t) else 0L}, .groups="drop")
    map_dfr(size_thrs, \(s) { k<-filter(cs,size>=s); if(!nrow(k)) return(tibble(min_size=s,purity=NA_real_))
    tibble(min_size=s, purity=sum(k$dom)/sum(k$size)) })
  }
  size_curves <- bind_rows(
    pub_clu %>% group_by(method) %>% group_modify(~purity_at_size(.x)) %>% ungroup(),
    sptce_focal %>% group_by(method) %>% group_modify(~purity_at_size(.x)) %>% ungroup()
  ) %>% mutate(method=order_methods(method))
  p17 <- ggplot(size_curves, aes(min_size,purity,colour=method,group=method)) +
    geom_line(linewidth=0.7) + geom_point(size=1.8) +
    scale_colour_manual(values=METHOD_COLOURS, drop=FALSE, name=NULL) +
    scale_x_continuous(breaks=size_thrs) + scale_y_continuous(labels=percent_format(), limits=c(0.4,1.01)) +
    labs(title="Purity vs minimum cluster size", x="Minimum cluster size", y="Purity") +
    guides(colour=guide_legend(ncol=3))
  save_plot("fig17_purity_vs_size_threshold.png", p17, 12, 8)

  focal_sens <- bind_rows(
    pub_clu %>% mutate(threshold=NA_real_),
    sptce_all
  ) %>% group_by(method, threshold) %>%
    group_modify(\(df, key) tibble(sensitivity=compute_sensitivity(df, CFG$focal_epitope))) %>%
    ungroup() %>%
    mutate(method=order_methods(method))
  p18 <- ggplot(focal_sens, aes(method,sensitivity,fill=method)) +
    geom_col(width=0.7) +
    scale_fill_manual(values=METHOD_COLOURS, drop=FALSE, guide="none") +
    scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
    labs(title=sprintf("Focal-epitope sensitivity (%s)", CFG$focal_epitope), x=NULL, y="Sensitivity") +
    theme(axis.text.x=element_text(angle=35,hjust=1,size=12))
  save_plot("fig18_focal_sensitivity.png", p18, 11, 7)
  invisible(NULL)
}


# =================================================================
# 12. STAGE 5b — cluster network visualisations
# =================================================================
run_stage5b_networks <- function(s5 = NULL, n_top = 12, ld_thr = 2) {
  cat("\n=== Stage 5b: cluster networks ===\n")
  if (is.null(s5)) s5 <- run_stage5()
  clu <- sptce_graph_nodes(s5$sptce_clusters)   # 0.75 A RMSD operating point only
  if (is.null(clu) || !nrow(clu)) { cat("  No clusters to plot\n"); return(invisible(NULL)) }
  if (!requireNamespace("igraph", quietly=TRUE) ||
      !requireNamespace("ggraph", quietly=TRUE) ||
      !requireNamespace("tidygraph", quietly=TRUE)) {
    cat("  igraph/ggraph/tidygraph not available; skipping\n"); return(invisible(NULL))
  }
  make_tcr_graph <- function(df, seq_col="cdr3b") {
    rows <- filter(df, !is.na(.data[[seq_col]]))
    if (nrow(rows) < 2) return(NULL)
    seqs <- rows[[seq_col]]
    dm   <- stringdistmatrix(seqs, seqs, method="lv")
    ei   <- which(dm <= ld_thr & upper.tri(dm), arr.ind=TRUE)
    if (!nrow(ei)) g <- igraph::make_empty_graph(n=nrow(rows), directed=FALSE)
    else {
      g <- igraph::graph_from_data_frame(
        d=data.frame(from=ei[,1], to=ei[,2], weight=1/(1+dm[ei])),
        vertices=data.frame(name=seq_len(nrow(rows))), directed=FALSE)
      g <- igraph::simplify(g, edge.attr.comb=list(weight="max"))
      bad <- is.na(igraph::E(g)$weight) | igraph::E(g)$weight <= 0
      if (any(bad)) g <- igraph::delete_edges(g, igraph::E(g)[bad])
    }
    igraph::V(g)$organism <- coalesce(EPITOPE_ORGANISM[rows$epitope],
                                      ifelse(is.na(rows$epitope),"Unlabelled","Other"))
    igraph::V(g)$cluster_id <- rows$cluster_id
    g
  }
  plot_net <- function(g, title="") {
    if (is.null(g)) return(ggplot()+theme_void()+labs(title=title))
    tg <- tidygraph::as_tbl_graph(g)
    ggraph::ggraph(tg, layout="fr") +
      ggraph::geom_edge_link(alpha=0.25, edge_width=0.3, colour="grey40") +
      ggraph::geom_node_point(aes(colour=organism), size=2.6, alpha=0.95) +
      scale_colour_manual(values=ORGANISM_PALETTE, name="Organism") +
      ggraph::theme_graph(base_family="sans") +
      labs(title=title) +
      theme(legend.position="right", plot.title=element_text(size=14,face="bold"))
  }
  top_ids <- clu %>% count(cluster_id, sort=TRUE) %>% head(n_top) %>% pull(cluster_id)
  plots <- map(top_ids, function(cid) {
    sub <- filter(clu, cluster_id==cid)
    dom <- { t<-table(sub$epitope[!is.na(sub$epitope)]); if(length(t)) names(t)[which.max(t)] else "unlab" }
    plot_net(make_tcr_graph(sub), sprintf("Cluster %s (n=%d, dom=%s)", substr(cid,1,12), nrow(sub), dom))
  })
  if (length(plots))
    save_plot("fig22_sptce_cluster_networks.png",
              patchwork::wrap_plots(plots, ncol=3, guides="collect") & theme(legend.position="right"),
              16, 14)
  focal_net <- filter(clu, cluster_id %in% top_ids) %>% make_tcr_graph()
  if (!is.null(focal_net))
    save_plot("fig22b_networks_by_organism.png", plot_net(focal_net,"Top SPACE-TCR clusters by epitope organism"), 12, 9)
  compact_df <- bind_rows(s5$published_clusters, clu) %>%
    filter(!is.na(cdr3b)) %>% group_by(method, cluster_id) %>% filter(n()>=3) %>%
    summarise(mean_ld=ld_stats(cdr3b)$mean_ld, size=n(), .groups="drop") %>%
    filter(!is.na(mean_ld)) %>% mutate(method=order_methods(method))
  if (nrow(compact_df)) {
    p24 <- ggplot(compact_df, aes(fct_rev(method), mean_ld, fill=method)) +
      geom_boxplot(outlier.size=0.5, linewidth=0.35) + coord_flip() +
      scale_fill_manual(values=METHOD_COLOURS, drop=FALSE, guide="none") +
      labs(title="Within-cluster LD per method", subtitle="Lower = tighter sequence convergence",
           x=NULL, y="Mean within-cluster LD")
    save_plot("fig24_cluster_compactness.png", p24, 11, 8)
  }
  .plot_cluster_overlays(s5, ld_thr = ld_thr)
  .plot_networks_overview(s5, ld_thr)
  .s5b_extra_figs(s5, ld_thr)
  cat("Stage 5b complete.\n")
  invisible(NULL)
}


# -----------------------------------------------------------------
# 12b. Cluster overlays (fig39 / 39b / 39c)
# -----------------------------------------------------------------
# The node set is EVERY TCR in the chosen universe carrying the focal
# epitope. Unclustered TCRs appear as pale grey, so each panel shows both the
# partition structure AND the retention of the method it is coloured by.
#
# One precomputed layout is used for both panels, so each node sits in the
# same place in both.
#
# "All TCRs in the database" is ambiguous and the choice matters, so it is an
# explicit argument:
#   "curated_db"  the Jouannet curated database (what the published methods
#                 were scored against). Includes TCRs SPACE-TCR never modelled,
#                 so SPACE-TCR's grey fraction mixes clustering failure with the
#                 modelling gap that stage 5d reports as n_modelled_gap.
#   "sptce_input" only the TCRs handed to SPACE-TCR: grey is then purely
#                 "modelled but not clustered".
#   "union"       both; widest, inherits the caveat above.
.focal_universe <- function(s5, focal_epi,
                            universe = c("curated_db", "sptce_input", "union")) {
  universe <- match.arg(universe)
  out <- list()
  if (universe %in% c("curated_db", "union")) {
    db_path <- file.path(CFG$jouannet_repo,
                         "Data/Database/database_pooled_human_2023_03_15.txt")
    if (file.exists(db_path)) {
      db <- tryCatch(load_jouannet_curated(db_path, cell_type = CFG$jouannet_cell_type,
                                           verbose = FALSE),
                     error = \(e) { message("  overlay: curated DB load failed: ",
                                            conditionMessage(e)); NULL })
      if (!is.null(db))
        out$curated_db <- db %>% filter(epitope == focal_epi) %>%
          distinct(cdr3b) %>% mutate(src_db = TRUE)
    } else message("  overlay: curated database not found at ", db_path)
  }
  if (universe %in% c("sptce_input", "union")) {
    if (file.exists(CFG$jouannet_labelled)) {
      inp <- read_tsv(CFG$jouannet_labelled, show_col_types = FALSE)
      if (all(c("cdr3_aa_beta", "epitope") %in% names(inp)))
        out$sptce_input <- inp %>% filter(epitope == focal_epi) %>%
          transmute(cdr3b = cdr3_aa_beta) %>% distinct() %>% mutate(src_input = TRUE)
    }
  }
  if (!length(out)) return(NULL)
  u <- reduce(out, full_join, by = "cdr3b") %>%
    filter(!is.na(cdr3b), nchar(cdr3b) >= 6) %>% distinct(cdr3b, .keep_all = TRUE)
  if (!"src_db"    %in% names(u)) u$src_db    <- NA
  if (!"src_input" %in% names(u)) u$src_input <- NA
  u
}


# Top n_named clusters get distinct hues; smaller multi-member clusters share a
# mid grey; unclustered nodes a pale grey.
#
# NOTE on the mid grey: a node whose cluster has only ONE member *within this
# epitope* lands in "Other cluster", not "Not clustered". The method did cluster
# it, but its cluster-mates carry a different epitope and are not on this graph.
# Mid grey = "clustered, structure not shown here"; pale grey = "not clustered".
# Only pale grey counts against retention.
.assign_overlay_factor <- function(cluster_id, n_named = 24L,
                                   label_other = "Other cluster",
                                   label_none  = "Not clustered",
                                   shorten = identity) {
  cid <- as.character(cluster_id)
  sizes <- sort(table(cid[!is.na(cid)]), decreasing = TRUE)
  multi <- names(sizes)[sizes >= 2]
  top   <- head(multi, n_named)
  pal   <- c(setNames(distinct_cluster_palette(length(top)), shorten(top)),
             setNames("#9E9E9E", label_other), setNames("#E8E8E8", label_none))
  f <- factor(dplyr::case_when(
    is.na(cid)     ~ label_none,
    cid %in% top   ~ shorten(cid),
    TRUE           ~ label_other), levels = names(pal))
  list(factor = f, palette = pal, n_named = length(top))
}


# -----------------------------------------------------------------
# 12c. Paired CDR3 alpha / beta sequence logos
# -----------------------------------------------------------------
# For the top SPACE-TCR clusters, the CDR3-alpha and CDR3-beta motifs side by
# side (one cluster per row: alpha | beta) so alpha/beta co-conservation can be
# read directly. Auto-run from run_stage5c_extras().
sptce_pairing_logos <- function(s5, epi_focus = CFG$focal_epitope, top_n = 8,
                                min_size = 5, min_purity = 0.8, file = NULL) {
  if (is.null(s5) || is.null(s5$sptce_clusters) || !nrow(s5$sptce_clusters)) {
    message("  pairing logos: no SPACE-TCR clusters; skipping."); return(invisible(NULL))
  }
  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    message("  pairing logos: ggseqlogo not installed; skipping."); return(invisible(NULL))
  }
  clu <- s5$sptce_clusters
  if (!"cdr3a" %in% names(clu) || all(is.na(clu$cdr3a) | clu$cdr3a == "")) {
    message("  pairing logos: no CDR3-alpha in clusters - add 'cdr3_aa_alpha' to ",
            "CFG$jouannet_labelled and re-run stage 5. Skipping.")
    return(invisible(NULL))
  }
  paired <- clu %>% filter(!is.na(cdr3a), cdr3a != "", !is.na(cdr3b), cdr3b != "", !is.na(epitope))
  if (!nrow(paired)) { message("  pairing logos: no paired rows; skipping."); return(invisible(NULL)) }
  summ <- paired %>% group_by(cluster_id) %>%
    summarise(size = n(), dom = names(sort(table(epitope), decreasing = TRUE))[1],
              dom_n = max(table(epitope)), .groups = "drop") %>%
    mutate(purity = dom_n / size)
  if (!is.null(epi_focus)) {
    summ <- filter(summ, dom == epi_focus, size >= min_size, purity >= min_purity)
    lab_title <- sprintf("Paired CDR3 alpha / beta motifs - monospecific %s clusters", epi_focus)
  } else {
    summ <- filter(summ, size >= min_size)
    lab_title <- "Paired CDR3 alpha / beta motifs - largest SPACE-TCR clusters"
  }
  summ <- summ %>% arrange(desc(size)) %>% head(top_n)
  if (!nrow(summ)) { message("  pairing logos: no qualifying clusters."); return(invisible(NULL)) }
  modal_set <- function(x) {
    x <- x[!is.na(x) & nchar(x) > 0]
    if (!length(x)) return(character(0))
    ml <- as.integer(names(sort(table(nchar(x)), decreasing = TRUE))[1])
    x[nchar(x) == ml]
  }
  seqs_list <- list()
  for (i in seq_len(nrow(summ))) {
    cid <- summ$cluster_id[i]
    sub <- filter(paired, cluster_id == cid)
    if (!is.null(epi_focus)) sub <- filter(sub, epitope == epi_focus)
    a_seqs <- modal_set(sub$cdr3a); b_seqs <- modal_set(sub$cdr3b)
    if (length(a_seqs) < 2 || length(b_seqs) < 2) next
    tag <- sprintf("%s  (%s, size %d, purity %.0f%%)", substr(cid, 1, 10),
                   summ$dom[i], summ$size[i], 100 * summ$purity[i])
    seqs_list[[paste0(tag, "\nalpha (n=", length(a_seqs), ")")]] <- a_seqs
    seqs_list[[paste0(tag, "\nbeta (n=",  length(b_seqs), ")")]] <- b_seqs
  }
  if (!length(seqs_list)) {
    message("  pairing logos: clusters too length-heterogeneous; skipping."); return(invisible(NULL))
  }
  p <- ggseqlogo::ggseqlogo(seqs_list, ncol = 2, method = "bits") +
    labs(title = lab_title,
         subtitle = "One cluster per row - left: CDR3 alpha, right: CDR3 beta (dominant-length sequences)") +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(size = 8))
  if (is.null(file))
    file <- if (is.null(epi_focus)) "fig37_pairing_logos_top.png"
            else sprintf("fig37_pairing_logos_%s.png", epi_focus)
  save_plot(file, p, 11, max(4, 2.3 * (length(seqs_list) / 2)))
  invisible(p)
}


# =================================================================
# 13. STAGE 5c — extra analyses
# =================================================================
run_stage5c_extras <- function(s5 = NULL) {
  cat("\n=== Stage 5c: extra analyses ===\n")
  if (is.null(s5)) s5 <- run_stage5()
  clu <- s5$sptce_clusters
  sp <- bind_rows(
    s5$published_clusters %>% group_by(method,cluster_id) %>%
      summarise(size=n(), purity={t<-table(epitope[!is.na(epitope)]);if(length(t))max(t)/sum(t) else NA_real_},.groups="drop"),
    clu %>% group_by(method,cluster_id) %>%
      summarise(size=n(), purity={t<-table(epitope[!is.na(epitope)]);if(length(t))max(t)/sum(t) else NA_real_},.groups="drop")
  ) %>% filter(size>=2, !is.na(purity)) %>% mutate(method=order_methods(method))
  p20b <- ggplot(sp, aes(size,purity,colour=method,fill=method)) +
    geom_smooth(method="loess",span=0.9,se=TRUE,linewidth=0.7,alpha=0.10) +
    scale_colour_manual(values=METHOD_COLOURS,drop=FALSE,name=NULL) +
    scale_fill_manual(values=METHOD_COLOURS,drop=FALSE,guide="none") +
    scale_x_log10(breaks=c(2,3,5,10,20,50,100,500)) +
    scale_y_continuous(labels=percent_format(),limits=c(0,1.01)) +
    labs(title="Cluster size vs purity - all methods overlaid", x="Cluster size (log10)",y="Cluster purity") +
    guides(colour=guide_legend(ncol=3, override.aes=list(fill=NA,linewidth=1)))
  save_plot("fig20b_size_purity_overlay.png", p20b, 12, 8)
  focal_epi <- CFG$focal_epitope
  per_clu <- clu %>% filter(!is.na(epitope)) %>% group_by(cluster_id) %>%
    summarise(size=n(), n_focal=sum(epitope==focal_epi), purity=n_focal/size,
              dom=names(sort(table(epitope),decreasing=TRUE))[1], .groups="drop") %>%
    filter(dom==focal_epi, size>=5, purity>=0.8) %>% arrange(desc(size)) %>% head(12)
  if (nrow(per_clu) && requireNamespace("ggseqlogo", quietly=TRUE)) {
    seqs_list <- map(per_clu$cluster_id, function(cid) {
      s <- filter(clu, cluster_id==cid, epitope==focal_epi) %>% pull(cdr3b)
      s <- s[!is.na(s)&nchar(s)>0]
      ml <- as.integer(names(sort(table(nchar(s)),decreasing=TRUE))[1])
      s[nchar(s)==ml]
    })
    names(seqs_list) <- sprintf("%s\n(n=%d, pur=%.0f%%)", substr(per_clu$cluster_id,1,12),
                                per_clu$size, 100*per_clu$purity)
    seqs_list <- seqs_list[lengths(seqs_list)>=2]
    if (length(seqs_list)) {
      p36c <- ggseqlogo::ggseqlogo(seqs_list, ncol=3, method="bits") +
        labs(title=sprintf("SPACE-TCR CDR3b motifs - monospecific %s clusters", focal_epi),
             subtitle=sprintf("%d qualifying clusters (size >= 5, purity >= 80%%)", nrow(per_clu)))
      save_plot(sprintf("fig36c_sptce_motifs_%s.png", focal_epi), p36c, 15, 12)
    }
  }
  sptce_pairing_logos(s5)
  .s5c_extra_figs(s5)
  cat("Stage 5c complete.\n")
  invisible(NULL)
}


# =================================================================
# 15. STAGE 6 — fair-cohort overlap diagnostic
# =================================================================
run_stage6 <- function(s5 = NULL) {
  cat("\n=== Stage 6: cohort overlap diagnostic ===\n")
  if (is.null(s5)) s5 <- run_stage5()
  sptce_cdr3 <- unique(s5$sptce_clusters$cdr3b)
  pub_cdr3   <- unique(s5$published_clusters$cdr3b)
  common     <- intersect(sptce_cdr3, pub_cdr3)
  cat("  SPACE-TCR:", length(sptce_cdr3), "| Published:", length(pub_cdr3),
      "| Common:", length(common), "\n")
  fair_clu <- filter(s5$sptce_clusters, cdr3b %in% common)
  fair_pub <- filter(s5$published_clusters, cdr3b %in% common)
  fair_sens <- bind_rows(fair_pub, fair_clu) %>% group_by(method) %>%
    group_modify(\(df, key) tibble(sensitivity=compute_sensitivity(df, CFG$focal_epitope))) %>%
    ungroup() %>% mutate(method=order_methods(method))
  p18f <- ggplot(fair_sens, aes(method,sensitivity,fill=method)) +
    geom_col(width=0.7) +
    scale_fill_manual(values=METHOD_COLOURS, drop=FALSE, guide="none") +
    scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
    labs(title=sprintf("Fair-cohort sensitivity (%s)", CFG$focal_epitope),
         subtitle=sprintf("Restricted to %d TCRs in both SPACE-TCR and published outputs", length(common)),
         x=NULL, y="Sensitivity") +
    theme(axis.text.x=element_text(angle=35,hjust=1,size=13))
  save_plot("fig18_fair_cohort_recovery.png", p18f, 11, 7)
  write_csv(tibble(metric=c("n_sptce","n_published","n_common"),
                   n=c(length(sptce_cdr3),length(pub_cdr3),length(common))),
            file.path(CFG$jouannet_out_dir,"cohort_overlap_summary.csv"))
  cat("Stage 6 complete.\n")
  invisible(list(common_cohort=common, fair_sensitivity=fair_sens))
}


# =================================================================
# 16. STAGE 7 — TCRBuilder2 training-data leakage check
# =================================================================
run_stage7_training_leakage <- function() {
  cat("\n=== Stage 7: TCRBuilder2 training-data leakage (parsed IgBLAST) ===\n")
  parsed_csv <- CFG$stcrdab_parsed_csv %||% .ext("igblast_tcr_parsed.csv")
  if (!file.exists(parsed_csv)) {
    cat("  Parsed CSV not found at:", parsed_csv, "; skipping\n"); return(invisible(NULL))
  }
  dir.create(.stcrdab_dir(), showWarnings = FALSE, recursive = TRUE)
  parsed <- read_csv(parsed_csv, show_col_types = FALSE)
  cat("  Parsed rows:", nrow(parsed), "\n")
  parsed <- parsed %>% filter(!is.na(cdr3), nchar(cdr3) >= 6)
  if ("v_gene" %in% names(parsed)) {
    beta <- parsed %>% filter(str_detect(v_gene, "TRBV"))
    if (nrow(beta) > 0) {
      cat("  Beta-chain CDR3s (TRBV):", nrow(beta), "of", nrow(parsed), "total\n")
      parsed <- beta
    } else cat("  No TRBV rows found; using all parsed CDR3s\n")
  }
  stcrdab_cdr3 <- unique(parsed$cdr3)
  stcrdab_cdr3 <- stcrdab_cdr3[!is.na(stcrdab_cdr3) & nchar(stcrdab_cdr3) >= 6]
  if (!file.exists(CFG$jouannet_labelled)) { cat("  jouannet_labelled missing\n"); return(invisible(NULL)) }
  sptce_input <- read_tsv(CFG$jouannet_labelled, show_col_types = FALSE)
  all_cdr3b   <- unique(c(sptce_input$cdr3_aa_beta %||% character(0)))
  all_cdr3b   <- all_cdr3b[!is.na(all_cdr3b)]
  leaked <- intersect(all_cdr3b, stcrdab_cdr3)
  cat(sprintf("  SPACE-TCR CDR3beta: %d | PDB CDR3beta: %d | Exact-match overlap: %d (%.2f%%)\n",
              length(all_cdr3b), length(stcrdab_cdr3), length(leaked),
              100 * length(leaked) / max(1,length(all_cdr3b))))
  epi_lookup <- if ("epitope" %in% names(sptce_input))
    setNames(sptce_input$epitope, sptce_input$cdr3_aa_beta)
  else setNames(rep(NA_character_, length(all_cdr3b)), all_cdr3b)
  leak_df <- tibble(cdr3b = all_cdr3b, in_stcrdab = all_cdr3b %in% stcrdab_cdr3,
                    epitope = epi_lookup[all_cdr3b])
  write_tsv(leak_df, file.path(.stcrdab_dir(), "leakage_per_tcr.tsv"))

  n_pdb_ids      <- n_distinct(str_extract(parsed$sequence_id, "^[^_]+"))
  n_unique_cdr3b <- length(stcrdab_cdr3)
  n_cohort       <- length(all_cdr3b)
  pct_overall    <- length(leaked) / max(1, n_cohort)
  struct_per_cdr3b <- parsed %>% filter(cdr3 %in% leaked) %>%
    distinct(sequence_id, cdr3) %>% count(cdr3, name = "n_struct")
  struct_lookup <- setNames(struct_per_cdr3b$n_struct, struct_per_cdr3b$cdr3)

  per_epi <- leak_df %>% group_by(epitope) %>%
    summarise(n_total = n(), n_in_stcrdab = sum(in_stcrdab),
              pct = n_in_stcrdab / n_total,
              n_struct = sum(struct_lookup[cdr3b[in_stcrdab]], na.rm = TRUE),
              .groups = "drop") %>% arrange(desc(pct))
  epi_counts <- tibble(
    set = c(rep("Jouannet_curated", nrow(per_epi)), rep("SPTCE_input", nrow(per_epi))),
    epitope = rep(per_epi$epitope, 2), n_total = rep(per_epi$n_total, 2),
    n_in_stcrdab = rep(per_epi$n_in_stcrdab, 2), pct_in_stcrdab = rep(per_epi$pct, 2))
  write_csv(epi_counts, file.path(.stcrdab_dir(), "stcrdab_overlap_per_epitope.csv"))
  write_csv(tibble(metric = "SPTCE input labels", n_total = n_cohort,
                   n_overlap = length(leaked), pct_overlap = 100 * pct_overall,
                   n_pdb_ids = n_pdb_ids, n_unique_cdr3b = n_unique_cdr3b),
            file.path(.stcrdab_dir(), "stcrdab_overlap_summary.csv"))

  n_top <- 12L
  per_epi_plot <- per_epi %>% filter(!is.na(epitope)) %>%
    slice_max(order_by = n_total, n = n_top, with_ties = FALSE) %>%
    mutate(epitope = fct_reorder(epitope, pct, .desc = TRUE),
           lab = sprintf("%d/%d\n(%.2f%%)\n[%d struct]", n_in_stcrdab, n_total,
                         100 * pct, n_struct))
  if (nrow(per_epi_plot) == 0) cat("  No labelled epitopes; skipping fig31\n")
  else {
    fills <- setNames(scales::hue_pal(l = 60, c = 70)(nlevels(per_epi_plot$epitope)),
                      levels(per_epi_plot$epitope))
    p31 <- ggplot(per_epi_plot, aes(epitope, pct, fill = epitope)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = lab), vjust = -0.25, size = 3.8, lineheight = 0.9) +
      scale_fill_manual(values = fills, guide = "none") +
      scale_y_continuous(labels = percent_format(accuracy = 0.01),
                         expand = expansion(mult = c(0, 0.22))) +
      labs(title = "TCRBuilder2 training-data overlap by epitope",
           subtitle = sprintf(paste0("Training set: %d PDB IDs total, %d unique human CDR3b ",
                                     "strings available. Cohort overlap: %d / %d = %.2f%%. Match: %s."),
                              n_pdb_ids, n_unique_cdr3b, length(leaked), n_cohort,
                              100 * pct_overall, CFG$stcrdab_match_mode %||% "cdr3_exact"),
           x = NULL, y = "Fraction of cohort TCRs found in training set",
           caption = paste0("Per-bar annotation: matched / total in cohort (% overlap) ",
                            "[N PDB structures whose CDR3b matches a cohort TCR of this epitope]")) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size=12),
            plot.caption = element_text(hjust = 0, size=14, colour = "grey45"))
    save_plot("fig31_stcrdab_overlap.png", p31, 13, 7)
  }
  cat("Stage 7 complete. Results in:", .stcrdab_dir(), "\n")
  invisible(list(leaked = leaked, per_epitope = per_epi))
}


# =================================================================
# 19. RESTORED FIGURES (stage 5 extras)
# =================================================================
.s5_figures_extra <- function(sptce_metrics, sptce_focal, sptce_all,
                              pub_clu, all_sptce, all_ld, sensitivity_all,
                              total_n, dih_map, sptce_dih,
                              sptce_dih_all_modes, members) {
  out <- CFG$jouannet_out_dir

  p4 <- sensitivity_all %>% filter(!is.na(sensitivity)) %>%
    mutate(method=order_methods(method),
           group=ifelse(str_detect(as.character(method),"SPACE-TCR"),"SPACE-TCR","Published")) %>%
    ggplot(aes(method, sensitivity, fill=group)) +
    geom_boxplot(outlier.shape=NA, alpha=0.75, linewidth=0.4) +
    geom_jitter(aes(colour=epitope), width=0.18, size=1.8, alpha=0.85) +
    scale_fill_manual(values=c("Published"="#9c9c9c","SPACE-TCR"="#e31a1c"), name=NULL) +
    scale_colour_manual(values=EPITOPE_PALETTE, name="Epitope") +
    scale_y_continuous(limits=c(0,NA), labels=percent_format()) +
    labs(title="Per-method sensitivity across 6 epitopes", x=NULL, y="Sensitivity") +
    theme(axis.text.x=element_text(angle=35,hjust=1))
  save_plot("fig04_sensitivity_box.png", p4, 11, 7)

  focal_ld <- keep_focal_rows(all_ld) %>% filter(!is.na(median_ld)) %>%
    mutate(method=order_methods(method))
  ld_cls <- focal_ld %>% mutate(class=factor(ld_class(median_ld), levels=CONV_LEVELS)) %>%
    count(method, class) %>% group_by(method) %>% mutate(pct=n/sum(n)) %>% ungroup()
  p6 <- ggplot(ld_cls, aes(fct_rev(method), pct, fill=class)) +
    geom_col(width=0.75, alpha=0.9) + coord_flip() +
    scale_fill_manual(values=PALETTES$ld_class) +
    scale_y_continuous(labels=percent_format(), expand=expansion(mult=c(0,0.01))) +
    labs(title="CDR3b LD convergence classes", x=NULL, y="Fraction of clusters", fill="CDR3b LD class") +
    theme(axis.text.y=element_text(size=13)) + guides(fill=guide_legend(nrow=2))
  save_plot("fig06_ld_convergence_stacked.png", p6, 12, 10)

  if (!is.null(sptce_dih)) {
    sptce_ld_focal <- keep_focal_rows(all_ld) %>%
      filter(!is.na(median_ld), str_detect(mode,"SPACE-TCR"))
    p7 <- ggplot(sptce_ld_focal, aes(median_ld, fill=mode)) +
      geom_histogram(binwidth=1, alpha=0.75, position="dodge", colour="white") +
      geom_vline(xintercept=c(1,3,6), linetype="dashed", colour="grey40", linewidth=0.4) +
      scale_fill_manual(values=PALETTES$mode_colour) + scale_x_continuous(breaks=0:15) +
      labs(title=sprintf("RMSD (%s A) vs RMSD+dih (%s A)", CFG$focal_threshold, .sptce_bubble_dih_thr),
           x="Median within-cluster CDR3b LD", y="Number of clusters", fill=NULL)
    save_plot("fig07_ld_rmsd_vs_dih.png", p7, 10, 5)
  }

  if (!is.null(dih_map)) {
    dist_dir   <- file.path(CFG$jouannet_sptce_out,"distance_tables")
    dist_files <- list.files(dist_dir,"\\.tsv\\.gz$",full.names=TRUE)
    labels     <- read_tsv(CFG$jouannet_labelled, show_col_types=FALSE) %>% select(id,cdr3_aa_beta,epitope)
    cat("  Dihedral sweep at", CFG$rmsd_sweep_fixed, "A...\n")
    dih_sweep <- map_dfr(CFG$dihedral_sweep, function(dc) {
      rows <- build_sptce_clusters(dist_files, labels, rmsd_cut=CFG$rmsd_sweep_fixed,
                                   dih_cut=dc, dih_map=dih_map, prefix=paste0("sw",dc))
      compute_metrics(rows, total_n) %>% mutate(dih_cutoff=dc)
    })
    rmsd_ref <- compute_metrics(filter(all_sptce, threshold==CFG$rmsd_sweep_fixed,
                                       mode=="SPACE-TCR RMSD"), total_n)
    write_csv(dih_sweep, file.path(out,"sptce_dihedral_sweep_1A.csv"))
    p8 <- ggplot(dih_sweep, aes(dih_cutoff)) +
      geom_hline(yintercept=rmsd_ref$retention, linetype="dashed", colour="#FC8D59", linewidth=0.7) +
      geom_hline(yintercept=rmsd_ref$purity,    linetype="dashed", colour="#2166AC", linewidth=0.7) +
      geom_line(aes(y=retention), colour="#FC8D59", linewidth=1.1) +
      geom_point(aes(y=retention), colour="#FC8D59", size=3) +
      geom_line(aes(y=purity), colour="#2166AC", linewidth=1.1) +
      geom_point(aes(y=purity), colour="#2166AC", size=3) +
      geom_vline(xintercept=CFG$focal_dihedral, linetype="dotted", colour="grey40", linewidth=0.5) +
      scale_x_continuous(breaks=CFG$dihedral_sweep) +
      scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
      labs(title=sprintf("Dihedral cutoff sweep (RMSD = %s A)", CFG$rmsd_sweep_fixed),
           x="Dihedral distance cutoff (radians)", y=NULL)
    save_plot("fig08_dihedral_sweep_1A.png", p8, 9, 6)
  }

  if (!is.null(sptce_dih_all_modes)) {
    all_mode_metrics <- sptce_dih_all_modes %>% group_by(threshold, mode) %>%
      group_modify(~compute_metrics(.x, total_n) %>% mutate(method=paste0(.y$mode," (",.y$threshold,"A)"))) %>%
      ungroup()
    p9 <- all_mode_metrics %>%
      pivot_longer(names(METRIC_LABELS), names_to="metric", values_to="value") %>%
      mutate(metric=factor(metric, names(METRIC_LABELS), METRIC_LABELS),
             mode_short=str_remove(mode,"^SPACE-TCR "), threshold_f=factor(threshold)) %>%
      ggplot(aes(mode_short, value, colour=threshold_f, group=threshold_f, shape=threshold_f)) +
      geom_line(linewidth=0.85, alpha=0.9) + geom_point(size=3) +
      facet_wrap(~metric, nrow=1, scales="free_y") +
      scale_colour_brewer(palette="RdYlBu", name="RMSD threshold (A)") +
      scale_shape_manual(values=c(16,17,15,18,8), name="RMSD threshold (A)") +
      scale_y_continuous(labels=percent_format()) +
      labs(title="Effect of dihedral prune across RMSD thresholds", x="SPACE-TCR mode", y=NULL) +
      theme(axis.text.x=element_text(size=12))
    save_plot("fig09_dih_modes_panel.png", p9, 15, 5)
  }

  bubble_min3 <- .s5_bubble_data(sptce_focal, pub_clu, all_sptce, total_n, 3)
  p10b <- bubble_min3 %>%
    pivot_longer(c(retention,purity,perc_clusters_gt09,perc_tcrs_gt09),
                 names_to="metric", values_to="value") %>%
    mutate(metric=factor(metric,names(METRIC_LABELS),METRIC_LABELS), method=order_methods(method)) %>%
    ggplot(aes(metric, fct_rev(method), size=value, fill=value)) +
    geom_point(shape=21, colour="grey30", stroke=0.3, alpha=0.92) +
    geom_text(aes(label=sprintf("%.2f",value)), vjust=-1.4, size=2.4, colour="grey20") +
    scale_size_area(max_size=14, guide="none") +
    scale_fill_gradient2(low="#F7F7F7",mid="#FEE08B",high="#D73027",midpoint=0.5,limits=c(0,1)) +
    labs(title="Performance metrics (cluster size >= 3)", x=NULL, y=NULL) +
    theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13),
          panel.grid.major.x=element_blank())
  save_plot("fig10b_bubble_metrics_size3plus.png", p10b, 11, 10)

  size_data_all <- bind_rows(
    filter(all_sptce, threshold==CFG$focal_threshold, mode=="SPACE-TCR RMSD") %>%
      group_by(cluster_id) %>% summarise(size=n(),.groups="drop") %>%
      mutate(method=name_sptce_rmsd(CFG$focal_threshold), colour_group="SPACE-TCR RMSD"),
    filter(all_sptce, threshold==.sptce_bubble_dih_thr, mode=="SPACE-TCR RMSD+dih") %>%
      group_by(cluster_id) %>% summarise(size=n(),.groups="drop") %>%
      mutate(method=name_sptce_focal_dih(.sptce_bubble_dih_thr), colour_group="SPACE-TCR RMSD+dih (focal)"),
    pub_clu %>% group_by(method,cluster_id) %>% summarise(size=n(),.groups="drop") %>%
      mutate(colour_group="published")
  ) %>% filter(size>=2) %>% mutate(method=order_methods(method))

  p12 <- ggplot(size_data_all, aes(log2(size), fct_rev(method), fill=colour_group)) +
    geom_violin(alpha=0.65, scale="width", linewidth=0.3, draw_quantiles=c(0.25,0.5,0.75)) +
    geom_boxplot(width=0.12, fill="white", alpha=0.85, outlier.size=0.4, linewidth=0.35) +
    geom_vline(xintercept=log2(c(2,5,10)), linetype="dashed", colour="grey50", linewidth=0.4) +
    scale_fill_manual(values=GROUP_PALETTE, name=NULL) +
    scale_x_continuous(breaks=log2(c(2,3,5,10,20,50,100,500)), labels=c(2,3,5,10,20,50,100,500)) +
    labs(title="Cluster size distribution (log2 scale)", x="Cluster size", y=NULL) +
    theme(axis.text.y=element_text(size=13))
  save_plot("fig12_cluster_size_violin.png", p12, 12, 9)

  p13 <- ggplot(size_data_all, aes(size, colour=colour_group, group=method, linetype=colour_group)) +
    stat_ecdf(linewidth=0.85, alpha=0.9) +
    geom_vline(xintercept=c(2,5,10), linetype="dotted", colour="grey40", linewidth=0.4) +
    scale_colour_manual(values=GROUP_PALETTE, name=NULL) +
    scale_linetype_manual(values=c("published"="solid","SPACE-TCR RMSD"="dashed",
                                   "SPACE-TCR RMSD+dih (focal)"="dashed"), guide="none") +
    scale_x_continuous(limits=c(2,50), breaks=c(2,5,10,20,50)) +
    scale_y_continuous(labels=percent_format()) +
    labs(title="Cumulative cluster size distribution", x="Cluster size", y="Cumulative fraction of clusters")
  save_plot("fig13_cluster_size_ecdf.png", p13, 10, 6)

  poly_stats <- bind_rows(
    filter(all_sptce, threshold==CFG$focal_threshold, mode=="SPACE-TCR RMSD") %>%
      select(cluster_id,epitope) %>% mutate(method=name_sptce_rmsd(CFG$focal_threshold), colour_group="SPACE-TCR RMSD"),
    pub_clu %>% select(method,cluster_id,epitope) %>% mutate(colour_group="published")
  ) %>% filter(!is.na(epitope)) %>%
    group_by(method,colour_group,cluster_id) %>%
    summarise(n_members=n(), n_epi=n_distinct(epitope), .groups="drop") %>%
    filter(n_members>=2) %>% mutate(method=order_methods(method))

  p15 <- ggplot(poly_stats, aes(n_epi, fct_rev(method), fill=colour_group)) +
    geom_violin(alpha=0.65, scale="width", linewidth=0.3, draw_quantiles=c(0.25,0.5,0.75)) +
    geom_boxplot(width=0.12, fill="white", alpha=0.85, outlier.size=0.4, linewidth=0.35) +
    geom_vline(xintercept=1.5, linetype="dashed", colour="grey30", linewidth=0.5) +
    scale_fill_manual(values=GROUP_PALETTE, name=NULL) +
    scale_x_log10(breaks=c(1,2,3,5,10,20,50,100,500)) +
    labs(title="Distinct epitopes per cluster (log scale)", x="Number of distinct epitopes", y=NULL) +
    theme(axis.text.y=element_text(size=13))
  save_plot("fig15_polyspecific_violin.png", p15, 12, 9)

  purity_lkp <- bind_rows(
    PUBLISHED_METRICS %>% select(method,purity) %>% mutate(colour_group="published"),
    filter(sptce_metrics, threshold==CFG$focal_threshold) %>%
      transmute(method, purity, colour_group=if_else(mode=="SPACE-TCR RMSD","SPACE-TCR RMSD","SPACE-TCR RMSD+dih (focal)"))
  )
  poly_summary <- poly_stats %>% group_by(method, colour_group) %>%
    summarise(n_clusters=n(), pct_poly=mean(n_epi>=2), .groups="drop") %>%
    inner_join(purity_lkp, by=c("method","colour_group"))
  p16 <- ggplot(poly_summary, aes(pct_poly, purity, colour=colour_group,
                                  size=colour_group=="published")) +
    geom_smooth(data=filter(poly_summary,colour_group=="published"), method="lm", se=TRUE,
                colour="grey50", linewidth=0.7, linetype="dashed", alpha=0.15) +
    geom_point(alpha=0.9) +
    geom_text(aes(label=method), size=2.4, vjust=-0.9, show.legend=FALSE) +
    scale_colour_manual(values=GROUP_PALETTE, name=NULL) +
    scale_size_manual(values=c(`TRUE`=3,`FALSE`=4), guide="none") +
    scale_x_continuous(labels=percent_format(), limits=c(0,1)) +
    scale_y_continuous(labels=percent_format(), limits=c(0.4,1.01)) +
    labs(title="% polyspecific clusters vs global purity",
         x="% polyspecific clusters", y="Global cluster purity")
  save_plot("fig16_polyspecific_vs_purity.png", p16, 10, 7)

  contam <- bind_rows(pub_clu, sptce_focal) %>%
    group_by(method, cluster_id) %>%
    summarise(size=n(), n_focal=sum(epitope==CFG$focal_epitope,na.rm=TRUE),
              n_other=sum(!is.na(epitope)&epitope!=CFG$focal_epitope),
              n_unlab=sum(is.na(epitope)), .groups="drop") %>%
    filter(n_focal>=2) %>% group_by(method) %>%
    summarise(total=sum(size), pct_focal=sum(n_focal)/total,
              pct_other=sum(n_other)/total, pct_unlab=sum(n_unlab)/total, .groups="drop") %>%
    pivot_longer(starts_with("pct_"), names_to="kind", values_to="frac") %>%
    mutate(kind=factor(recode(kind,"pct_focal"="Focal epitope",
                              "pct_other"="Other labelled","pct_unlab"="Unlabelled"),
                       levels=c("Focal epitope","Other labelled","Unlabelled")),
           method=order_methods(method))
  p19 <- ggplot(contam, aes(method, frac, fill=kind)) +
    geom_col(width=0.75, colour="white", linewidth=0.3) +
    scale_fill_manual(values=c("Focal epitope"="#1a9850","Other labelled"="#d73027",
                               "Unlabelled"="#cccccc"), name=NULL) +
    scale_y_continuous(labels=percent_format()) +
    labs(title="Composition of focal-epitope clusters", x=NULL, y=NULL) +
    theme(axis.text.x=element_text(angle=35,hjust=1,size=13))
  save_plot("fig19_focal_cluster_composition.png", p19, 11, 7)

  scatter <- bind_rows(
    pub_clu %>% group_by(method,cluster_id) %>%
      summarise(size=n(), purity={t<-table(epitope[!is.na(epitope)]);if(length(t))max(t)/sum(t) else NA_real_}, .groups="drop"),
    sptce_focal %>% group_by(method,cluster_id) %>%
      summarise(size=n(), purity={t<-table(epitope[!is.na(epitope)]);if(length(t))max(t)/sum(t) else NA_real_}, .groups="drop")
  ) %>% filter(size>=2, !is.na(purity)) %>% mutate(method=order_methods(method))
  p20 <- ggplot(scatter, aes(size, purity, colour=method)) +
    geom_jitter(alpha=0.35, size=0.9, width=0.05, height=0.005) +
    geom_smooth(method="loess", span=0.9, se=FALSE, linewidth=0.7) +
    scale_colour_manual(values=METHOD_COLOURS, drop=FALSE, name=NULL) +
    scale_x_log10(breaks=c(2,3,5,10,20,50,100,500)) +
    scale_y_continuous(labels=percent_format(), limits=c(0,1.01)) +
    facet_wrap(~method, ncol=4) +
    labs(title="Per-cluster size vs purity", x="Cluster size (log10)", y="Cluster purity") +
    theme(legend.position="none", strip.text=element_text(size=12))
  save_plot("fig20_size_vs_purity_scatter.png", p20, 14, 10)

  radar <- bind_rows(
    PUBLISHED_METRICS %>% mutate(n_clusters=9999L, colour_group="published"),
    filter_sptce_focal(all_sptce) %>% group_by(mode, threshold) %>%
      group_modify(~compute_metrics(.x,total_n) %>%
                     mutate(method=paste0(.y$mode," (",.y$threshold,"A)"))) %>% ungroup()
  ) %>% select(method, retention, purity, perc_clusters_gt09, perc_tcrs_gt09) %>%
    mutate(method=order_methods(method)) %>%
    pivot_longer(-method, names_to="metric", values_to="value") %>%
    mutate(metric=factor(recode(metric, retention="Retention", purity="Purity",
                                perc_clusters_gt09="% clusters\n>0.9", perc_tcrs_gt09="% TCRs\nin >0.9"),
                         levels=c("Retention","Purity","% clusters\n>0.9","% TCRs\nin >0.9")))
  p21 <- ggplot(radar, aes(metric, value, fill=method)) +
    geom_col(width=0.85) + coord_polar() +
    scale_fill_manual(values=METHOD_COLOURS, drop=FALSE, guide="none") +
    scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
    facet_wrap(~method, ncol=4) +
    labs(title="Per-method metric summary", x=NULL, y=NULL) +
    theme(strip.text=element_text(size=13,face="bold"), axis.text.x=element_text(size=12))
  save_plot("fig21_method_radar.png", p21, 14, 11)

  # fig32 — full cohort vs training-held-out. Skipped if Stage 7 has not
  # written its files.
  leakage_f <- file.path(.stcrdab_dir(),"leakage_per_tcr.tsv")
  if (file.exists(leakage_f)) {
    cat("  fig32: full-cohort vs training-held-out metrics\n")
    leaked_cdr3 <- read_tsv(leakage_f, show_col_types=FALSE) %>% filter(in_stcrdab) %>% pull(cdr3b)
    metric_keys <- c("retention","purity","perc_clusters_gt09","perc_tcrs_gt09")
    mode_labels <- c("SPACE-TCR RMSD"="SPTCE RMSD", "SPACE-TCR RMSD+dih"="SPTCE RMSD+dih")
    f32_df <- map_dfr(names(mode_labels), function(m) {
      thr_m <- if (m=="SPACE-TCR RMSD+dih") .sptce_bubble_dih_thr else CFG$focal_threshold
      full <- filter(all_sptce, threshold==thr_m, mode==m)
      if (!nrow(full)) return(NULL)
      held <- filter(full, !cdr3_aa_beta %in% leaked_cdr3)
      bind_rows(compute_metrics(full, total_n) %>% mutate(cohort="Full cohort"),
                compute_metrics(held, total_n) %>% mutate(cohort="Training-held-out")) %>%
        mutate(mode_lab = mode_labels[[m]])
    })
    if (nrow(f32_df)) {
      f32_long <- f32_df %>% select(mode_lab, cohort, all_of(metric_keys)) %>%
        pivot_longer(all_of(metric_keys), names_to="metric", values_to="value") %>%
        mutate(metric = factor(recode(metric, retention="Retention", purity="Purity",
                                      perc_clusters_gt09="% clusters\npurity > 0.9",
                                      perc_tcrs_gt09="% TCRs in\nhigh-purity clusters"),
                               levels=c("Retention","Purity","% clusters\npurity > 0.9",
                                        "% TCRs in\nhigh-purity clusters")),
               cohort = factor(cohort, levels=c("Full cohort","Training-held-out")),
               mode_lab = factor(mode_lab, levels=c("SPTCE RMSD","SPTCE RMSD+dih")))
      p32 <- ggplot(f32_long, aes(cohort, value, colour=mode_lab, group=mode_lab)) +
        geom_line(linewidth=0.9) + geom_point(size=2.8) +
        geom_text(aes(label=sprintf("%.2f", value)),
                  hjust=c(1.25, -0.25)[as.integer(f32_long$cohort)], size=4.0, show.legend=FALSE) +
        facet_wrap(~metric, nrow=1, scales="free_y") +
        scale_colour_manual(values=c("SPTCE RMSD"="#595959","SPTCE RMSD+dih"="#4ab04a"), name=NULL) +
        scale_y_continuous(labels=percent_format(accuracy=0.01)) +
        scale_x_discrete(labels=c("Full cohort"="Full", "Training-held-out"="Held-out"),
                         expand=expansion(mult=c(0.30,0.30))) +
        labs(title=sprintf("SPTCE metrics: full cohort vs training-held-out (%s A)", CFG$focal_threshold),
             subtitle="If held-out values drop materially, the full-cohort numbers are inflated by training leakage.",
             x=NULL, y=NULL) +
        theme(panel.grid.major.x=element_blank(), axis.text.x=element_text(size=12),
              panel.spacing=unit(1.1,"lines"), plot.margin=margin(6, 18, 6, 6))
      save_plot("fig32_metrics_full_vs_heldout.png", p32, 15, 6)
    }
  }
  invisible(NULL)
}


# ── fig25_networks_overview: dot-in-cluster panels ─────────────
.plot_networks_overview <- function(s5, ld_thr = 2, focal_epi = CFG$focal_epitope,
                                    min_cluster_size = 2L) {
  cat("  fig25_networks_overview: dot-in-cluster plot...\n")
  if (!requireNamespace("ggforce", quietly=TRUE)) {
    cat("    ggforce not available; skipping networks overview\n"); return(invisible(NULL))
  }
  sptce_all <- s5$sptce_all_thresholds %||% s5$sptce_clusters
  pub_clu   <- s5$published_clusters
  if (is.null(sptce_all) || is.null(pub_clu)) return(invisible(NULL))
  get_sptce <- function(thr) sptce_all %>%
    filter(abs(threshold - thr) < 1e-6,
           str_detect(method, regex("RMSD", ignore_case=TRUE)),
           !str_detect(method, regex("dih",  ignore_case=TRUE)))
  sptce_075 <- get_sptce(CFG$focal_threshold)
  sptce_1   <- get_sptce(1.0)
  if (!nrow(sptce_1)) {
    alt <- sort(unique(sptce_all$threshold)); alt <- alt[alt > CFG$focal_threshold]
    if (length(alt)) sptce_1 <- get_sptce(alt[1])
  }
  clustcr_clu <- filter(pub_clu, method %in% c("clusTCR","ClusTCR"))
  gliph2_clu  <- filter(pub_clu, str_detect(method, regex("GLIPH2", ignore_case=TRUE)))

  all_epi <- na.omit(unique(c(sptce_075$epitope, sptce_1$epitope,
                              clustcr_clu$epitope, gliph2_clu$epitope)))
  extra_epi  <- setdiff(all_epi, names(EPITOPE_PALETTE))
  extra_cols <- if (length(extra_epi))
    setNames(scales::hue_pal(l=55,c=80)(length(extra_epi)), extra_epi) else character(0)
  epi_pal <- c(EPITOPE_PALETTE[intersect(names(EPITOPE_PALETTE), all_epi)],
               extra_cols, "Unknown" = "#555555")
  GOLD <- unname(EPITOPE_PALETTE["GILGFVFTL"])

  sunflower_pos <- function(n) {
    if (n == 1) return(data.frame(dx=0, dy=0))
    golden <- pi * (3 - sqrt(5)); k <- seq_len(n)
    r <- sqrt((k - 0.5) / (n - 0.5)); theta <- k * golden
    data.frame(dx = r * cos(theta), dy = r * sin(theta))
  }
  layout_tcrs <- function(df, dot_r = 0.5, pad = 0.85) {
    df <- df %>% filter(!is.na(cdr3b)) %>%
      group_by(cluster_id) %>% filter(n() >= min_cluster_size) %>% ungroup()
    if (!nrow(df)) return(NULL)
    clu_info <- df %>% group_by(cluster_id) %>% summarise(n=n(), .groups="drop") %>%
      arrange(desc(n)) %>% mutate(R = sqrt(n) * dot_r + pad)
    n_clu <- nrow(clu_info)
    set.seed(CFG$seed)
    cx <- runif(n_clu, -25, 25); cy <- runif(n_clu, -25, 25); R <- clu_info$R
    for (iter in seq_len(300)) {
      dx_mat <- outer(cx, cx, "-"); dy_mat <- outer(cy, cy, "-")
      d_mat  <- sqrt(dx_mat^2 + dy_mat^2)
      md_mat <- outer(R, R, "+") + 0.6
      diag(d_mat) <- Inf
      overlap <- d_mat < md_mat
      force   <- ifelse(overlap & d_mat > 1e-6, (md_mat - d_mat) / d_mat * 0.45, 0)
      cx <- cx + rowSums(force * dx_mat); cy <- cy + rowSums(force * dy_mat)
    }
    clu_info$cx <- cx; clu_info$cy <- cy
    df %>% left_join(clu_info %>% select(cluster_id,cx,cy,R,n), by="cluster_id") %>%
      group_by(cluster_id) %>%
      group_modify(~ {
        sp <- sunflower_pos(nrow(.x))
        # pack dots into 78% of the inner radius so the cloud sits well inside
        # the boundary circle; the margin clears the rendered dot radius
        fill_r  <- (first(.x$R) - pad) * 0.78
        scale_f <- fill_r / max(max(abs(sp$dx)), max(abs(sp$dy)), 1e-6)
        .x$x <- first(.x$cx) + sp$dx * scale_f
        .x$y <- first(.x$cy) + sp$dy * scale_f
        .x
      }) %>% ungroup()
  }
  cluster_circles <- function(df_layout) df_layout %>% group_by(cluster_id) %>%
    summarise(cx=first(cx), cy=first(cy), R=first(R), .groups="drop")
  # geom_point size is a fixed physical size (mm) but the layout auto-scales to
  # the data range, so one hard-coded size overflows the denser panels.
  panel_point_size <- function(circles, df, plot_area_mm = 150) {
    span_x <- diff(range(c(circles$cx - circles$R, circles$cx + circles$R)))
    span_y <- diff(range(c(circles$cy - circles$R, circles$cy + circles$R)))
    mm_per_unit <- plot_area_mm / max(span_x, span_y, 1)
    max(2.4, min(4.5, 0.7 * mm_per_unit))
  }
  make_panel <- function(df_raw, title, subtitle, highlight=FALSE) {
    df <- layout_tcrs(df_raw)
    if (is.null(df) || !nrow(df)) return(ggplot()+theme_void()+labs(title=title))
    circles <- cluster_circles(df)
    pt_sz <- panel_point_size(circles, df)
    nc <- n_distinct(df$cluster_id); nt <- nrow(df)
    df <- df %>% mutate(
      fill_col = if (highlight) ifelse(!is.na(epitope) & epitope==focal_epi, focal_epi, "Other")
                 else ifelse(is.na(epitope), "Unknown", epitope),
      alpha_v = if (highlight) ifelse(!is.na(epitope) & epitope==focal_epi, 0.95, 0.45) else 0.88)
    pal <- if (highlight) c(setNames(GOLD, focal_epi), "Other"="#3a3a3a") else epi_pal
    ggplot() +
      ggforce::geom_circle(data=circles, aes(x0=cx, y0=cy, r=R),
                           fill=NA, colour="grey45", linewidth=0.5, inherit.aes=FALSE) +
      geom_point(data=df, aes(x, y, fill=fill_col, alpha=alpha_v),
                 shape=21, size=pt_sz, colour="white", stroke=0.15) +
      scale_fill_manual(values=pal, na.value="#555555", name="Epitope") +
      scale_alpha_identity() + coord_fixed() + theme_void() +
      theme(legend.position="right", legend.text=element_text(size=13),
            legend.key.size=unit(0.5,"cm"), legend.title=element_text(size=15, face="bold"),
            plot.title=element_text(face="bold", size=19, margin=margin(b=1)),
            plot.subtitle=element_text(size=14, colour="grey50"),
            plot.margin=margin(6,4,4,6)) +
      labs(title=sprintf("%s - %d clusters, %d TCRs", title, nc, nt), subtitle=subtitle)
  }
  sub_full <- "Each dot = one TCR, colour = its epitope. Circle = cluster boundary."
  sub_hl   <- sprintf("%s highlighted, rest faded.", focal_epi)
  lab1 <- sprintf("SPACE-TCR RMSD (%s A)",
                  if(nrow(sptce_1)) round(unique(sptce_1$threshold)[1],2) else "1.0")
  panels <- list(
    sp075_full = make_panel(sptce_075,  sprintf("SPACE-TCR RMSD (%s A)", CFG$focal_threshold), sub_full),
    sp1_full   = make_panel(sptce_1,    lab1, sub_full),
    cl_full    = make_panel(clustcr_clu,"clusTCR", sub_full),
    sp075_hl   = make_panel(sptce_075,  sprintf("SPACE-TCR RMSD (%s A)", CFG$focal_threshold), sub_hl, TRUE),
    sp1_hl     = make_panel(sptce_1,    lab1, sub_hl, TRUE),
    cl_hl      = make_panel(clustcr_clu,"clusTCR", sub_hl, TRUE))
  fig <- (panels$sp075_full | panels$sp075_hl) /
         (panels$sp1_full   | panels$sp1_hl)   /
         (panels$cl_full    | panels$cl_hl)    +
    patchwork::plot_annotation(
      title    = "Cluster network overview (focal threshold)",
      subtitle = sprintf("Min cluster size = %d.", min_cluster_size),
      theme    = theme(plot.title=element_text(face="bold", size=24, hjust=0.5),
                       plot.subtitle=element_text(size=16, colour="grey40", hjust=0.5)))
  save_plot("fig25_networks_overview.png", fig, 16, 20)
  invisible(NULL)
}


# ── Stage 5b extra figures ─────────────────────────────────────
.s5b_extra_figs <- function(s5, ld_thr=2) {
  clu <- sptce_graph_nodes(s5$sptce_clusters)
  if (is.null(clu)||!nrow(clu)) return(invisible(NULL))
  if (!all(sapply(c("igraph","ggraph","tidygraph"), requireNamespace, quietly=TRUE))) {
    message("  igraph/ggraph/tidygraph needed for network figs; skipping"); return(invisible(NULL))
  }
  make_tcr_graph <- function(df, seq_col="cdr3b") {
    rows <- filter(df, !is.na(.data[[seq_col]]))
    if (nrow(rows)<2) return(NULL)
    seqs <- rows[[seq_col]]; dm <- stringdistmatrix(seqs,seqs,method="lv")
    ei <- which(dm<=ld_thr&upper.tri(dm),arr.ind=TRUE)
    if (!nrow(ei)) g <- igraph::make_empty_graph(n=nrow(rows),directed=FALSE)
    else {
      g <- igraph::graph_from_data_frame(
        data.frame(from=ei[,1],to=ei[,2],weight=1/(1+dm[ei])),
        vertices=data.frame(name=seq_len(nrow(rows))),directed=FALSE)
      g <- igraph::simplify(g,edge.attr.comb=list(weight="max"))
      bad <- is.na(igraph::E(g)$weight)|igraph::E(g)$weight<=0
      if (any(bad)) g <- igraph::delete_edges(g,igraph::E(g)[bad])
    }
    igraph::V(g)$organism <- coalesce(EPITOPE_ORGANISM[rows$epitope],
                                      ifelse(is.na(rows$epitope),"Unlabelled","Other"))
    igraph::V(g)$cluster_id <- rows$cluster_id; g
  }
  plot_net <- function(g, title="") {
    if (is.null(g)) return(ggplot()+theme_void()+labs(title=title))
    tg <- tidygraph::as_tbl_graph(g)
    ggraph::ggraph(tg, layout="fr") +
      ggraph::geom_edge_link(alpha=0.25,edge_width=0.3,colour="grey40") +
      ggraph::geom_node_point(aes(colour=organism),size=2.6,alpha=0.95) +
      scale_colour_manual(values=ORGANISM_PALETTE,name="Organism") +
      ggraph::theme_graph(base_family="sans") +
      labs(title=title) + theme(legend.position="right",plot.title=element_text(size=14,face="bold"))
  }
  top_ids <- clu %>% count(cluster_id,sort=TRUE) %>% head(12) %>% pull(cluster_id)
  focal_clu <- filter(clu, cluster_id %in% top_ids, epitope==CFG$focal_epitope)
  if (nrow(focal_clu)>=2)
    save_plot("fig22b_focal_epitope_network.png",
              plot_net(make_tcr_graph(focal_clu), sprintf("Focal-epitope network (%s)", CFG$focal_epitope)),
              10, 8)
  biggest <- filter(clu, cluster_id==top_ids[1])
  save_plot("fig23_largest_cluster.png",
            plot_net(make_tcr_graph(biggest), sprintf("Largest cluster (n=%d)", nrow(biggest))), 10, 9)
  degs <- map_dfr(top_ids, function(cid) {
    sub <- filter(clu, cluster_id==cid)
    g <- make_tcr_graph(sub); if(is.null(g)) return(NULL)
    tibble(cluster_id=cid, degree=igraph::degree(g))
  })
  if (nrow(degs)) {
    p25 <- ggplot(degs, aes(degree)) +
      geom_histogram(binwidth=1, fill=unname(GROUP_PALETTE["SPACE-TCR RMSD"]), colour="white", linewidth=0.3) +
      facet_wrap(~cluster_id, scales="free_y", ncol=4) +
      labs(title="Per-cluster node-degree distributions",
           x="Degree (neighbours within LD threshold)", y="Count") +
      theme(strip.text=element_text(size=12))
    save_plot("fig25_degree_distributions.png", p25, 14, 9)
  }
  invisible(NULL)
}


# ── Stage 5c extra figures ─────────────────────────────────────
.s5c_extra_figs <- function(s5) {
  clu <- s5$sptce_clusters
  if (is.null(clu)||!nrow(clu)) return(invisible(NULL))
  poly <- clu %>% filter(!is.na(epitope)) %>% group_by(cluster_id) %>%
    summarise(size=n(), n_epi=n_distinct(epitope),
              dom_epi=names(sort(table(epitope),decreasing=TRUE))[1], .groups="drop") %>%
    filter(size>=5, n_epi>=2) %>% arrange(desc(n_epi),desc(size)) %>% head(12)
  if (nrow(poly) && requireNamespace("ggseqlogo",quietly=TRUE)) {
    seqs_list <- map(poly$cluster_id, function(cid) {
      s <- filter(clu,cluster_id==cid) %>% pull(cdr3b); s <- s[!is.na(s)]
      ml <- as.integer(names(sort(table(nchar(s)),decreasing=TRUE))[1]); s[nchar(s)==ml]
    })
    names(seqs_list) <- sprintf("%s\n(size=%d, epi=%d, dom=%s)",
                                substr(poly$cluster_id,1,10),poly$size,poly$n_epi,poly$dom_epi)
    seqs_list <- seqs_list[lengths(seqs_list)>=2]
    if (length(seqs_list)) {
      p36 <- ggseqlogo::ggseqlogo(seqs_list,ncol=3,method="bits") +
        labs(title="CDR3b motifs of polyspecific SPACE-TCR clusters",
             subtitle="Clusters with size >= 5 and >= 2 distinct epitopes")
      save_plot("fig36_polyspecific_motifs.png", p36, 15, 12)
    }
    comp <- filter(clu, cluster_id %in% poly$cluster_id, !is.na(epitope)) %>%
      count(cluster_id, epitope) %>% group_by(cluster_id) %>%
      mutate(frac=n/sum(n)) %>% ungroup() %>%
      mutate(cluster_id=factor(cluster_id, levels=poly$cluster_id))
    p36b <- ggplot(comp, aes(cluster_id, frac, fill=epitope)) +
      geom_col(width=0.8, colour="white", linewidth=0.2) +
      scale_fill_manual(values=EPITOPE_PALETTE, na.value="#999999", name="Epitope") +
      scale_y_continuous(labels=percent_format(), expand=expansion(mult=c(0,0.01))) +
      coord_flip() +
      labs(title="Epitope composition of polyspecific clusters", x=NULL, y="Fraction of cluster members") +
      theme(axis.text.y=element_text(size=12))
    save_plot("fig36b_polyspecific_composition.png", p36b, 11, 8)
  }
  co_pairs <- clu %>% filter(!is.na(epitope)) %>% group_by(cluster_id) %>%
    filter(n_distinct(epitope) >= 2) %>%
    summarise(epis=list(sort(unique(epitope))), .groups="drop")
  if (nrow(co_pairs)) {
    pair_tbl <- map_dfr(co_pairs$epis, function(e) {
      if (length(e) < 2) return(NULL)
      cb <- combn(e, 2); tibble(a=cb[1,], b=cb[2,])
    })
    if (nrow(pair_tbl)) {
      pair_counts <- bind_rows(pair_tbl, rename(pair_tbl, a=b, b=a)) %>% count(a, b, name="n_shared")
      p37 <- ggplot(pair_counts, aes(a, b, fill=n_shared)) +
        geom_tile(colour="white", linewidth=0.4) +
        geom_text(aes(label=n_shared), size=2.8) +
        scale_fill_gradient(low="#deebf7", high="#08306b", name="Shared\nclusters") +
        scale_x_discrete(guide=guide_axis(angle=45)) +
        labs(title="Epitope co-occurrence in polyspecific clusters", x=NULL, y=NULL) +
        theme(panel.grid=element_blank())
      save_plot("fig37_epitope_cooccurrence.png", p37, 10, 9)
    }
  }
  focal_clu <- clu %>% group_by(cluster_id) %>%
    summarise(size=n(), n_focal=sum(epitope==CFG$focal_epitope, na.rm=TRUE),
              purity={t<-table(epitope[!is.na(epitope)]); if(length(t)) max(t)/sum(t) else NA_real_},
              dom=if(any(!is.na(epitope))) names(sort(table(epitope[!is.na(epitope)]),decreasing=TRUE))[1] else NA,
              .groups="drop") %>%
    filter(size>=2, dom==CFG$focal_epitope, !is.na(purity))
  if (nrow(focal_clu)) {
    p38 <- ggplot(focal_clu, aes(size, purity, colour=purity>=0.9)) +
      geom_jitter(width=0.04, height=0.005, alpha=0.7, size=1.6) +
      scale_colour_manual(values=c(`TRUE`="#1a9850",`FALSE`="#d73027"),
                          labels=c(`TRUE`="High purity (>=0.9)",`FALSE`="Contaminated (<0.9)"), name=NULL) +
      scale_x_log10(breaks=c(2,3,5,10,20,50)) +
      scale_y_continuous(labels=percent_format(), limits=c(0,1.01)) +
      labs(title=sprintf("Purity vs size for %s-dominant SPACE-TCR clusters", CFG$focal_epitope),
           x="Cluster size (log10)", y="Cluster purity")
    save_plot("fig38_focal_purity_vs_size.png", p38, 10, 7)
  }
  invisible(NULL)
}


# ── Stage 6b — cohort difference figures (figs 26-30) ──────────
run_stage6b_cohort_diff <- function(s5 = NULL) {
  cat("\n=== Stage 6b: cohort difference figures ===\n")
  if (is.null(s5)) s5 <- run_stage5()
  sptce_cdr3 <- unique(s5$sptce_clusters$cdr3b)
  pub_cdr3   <- unique(s5$published_clusters$cdr3b)
  common     <- intersect(sptce_cdr3, pub_cdr3)
  sptce_only <- setdiff(sptce_cdr3, pub_cdr3)
  pub_only   <- setdiff(pub_cdr3, sptce_cdr3)

  overlap_df <- tibble(
    set=factor(c("SPACE-TCR only","Shared","Published only"),
               levels=c("SPACE-TCR only","Shared","Published only")),
    n=c(length(sptce_only), length(common), length(pub_only)))
  p26 <- ggplot(overlap_df, aes(set, n, fill=set)) +
    geom_col(width=0.7) + geom_text(aes(label=n), vjust=-0.4, size=3.5) +
    scale_fill_manual(values=c("SPACE-TCR only"="#e31a1c","Shared"="#9c27b0","Published only"="#9c9c9c"), guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
    labs(title="CDR3beta cohort overlap: SPACE-TCR vs published", x=NULL, y="Unique CDR3beta")
  save_plot("fig26_cohort_overlap_bar.png", p26, 9, 6)

  epi_comp <- bind_rows(
    s5$sptce_clusters %>% distinct(cdr3b, epitope) %>% mutate(src="SPACE-TCR"),
    s5$published_clusters %>% distinct(cdr3b, epitope) %>% mutate(src="Published")
  ) %>% filter(!is.na(epitope), epitope %in% CFG$jouannet_epitopes) %>%
    distinct(src, cdr3b, epitope) %>% count(src, epitope) %>%
    mutate(epitope=factor(epitope, CFG$jouannet_epitopes))
  p27 <- ggplot(epi_comp, aes(epitope, n, fill=src)) +
    geom_col(position="dodge", width=0.7) +
    scale_fill_manual(values=c("SPACE-TCR"="#e31a1c","Published"="#9c9c9c"), name=NULL) +
    labs(title="Per-epitope cohort size by source", x=NULL, y="Unique CDR3beta") +
    theme(axis.text.x=element_text(angle=20,hjust=1))
  save_plot("fig27_per_epitope_cohort.png", p27, 11, 6)

  len_df <- bind_rows(
    tibble(cdr3b=sptce_only, set="SPACE-TCR only"),
    tibble(cdr3b=common,     set="Shared"),
    tibble(cdr3b=pub_only,   set="Published only")
  ) %>% filter(!is.na(cdr3b)) %>%
    mutate(len=nchar(cdr3b), set=factor(set, levels=c("SPACE-TCR only","Shared","Published only")))
  p28 <- ggplot(len_df, aes(len, fill=set)) +
    geom_histogram(binwidth=1, position="identity", alpha=0.55, colour="white", linewidth=0.2) +
    scale_fill_manual(values=c("SPACE-TCR only"="#e31a1c","Shared"="#9c27b0","Published only"="#9c9c9c"), name=NULL) +
    scale_x_continuous(breaks=seq(6,24,2)) +
    labs(title="CDR3beta length by cohort membership", x="CDR3beta length (aa)", y="Count")
  save_plot("fig28_cdr3_length_by_cohort.png", p28, 10, 6)

  ret_df <- tibble(src=c("SPACE-TCR","Published"),
                   frac=c(mean(common %in% s5$sptce_clusters$cdr3b),
                          mean(common %in% s5$published_clusters$cdr3b)))
  p29 <- ggplot(ret_df, aes(src, frac, fill=src)) +
    geom_col(width=0.6) +
    geom_text(aes(label=percent(frac, accuracy=0.1)), vjust=-0.4, size=3.5) +
    scale_fill_manual(values=c("SPACE-TCR"="#e31a1c","Published"="#9c9c9c"), guide="none") +
    scale_y_continuous(labels=percent_format(), limits=c(0,1.05), expand=expansion(mult=c(0,0.05))) +
    labs(title="Fraction of shared CDR3beta retained in clusters", x=NULL, y="Fraction clustered")
  save_plot("fig29_shared_retention.png", p29, 8, 6)

  cohorts <- list(
    "Full SPACE-TCR" = s5$sptce_clusters,
    "Shared only"    = filter(s5$sptce_clusters, cdr3b %in% common),
    "SPACE-TCR only" = filter(s5$sptce_clusters, cdr3b %in% sptce_only))
  sens_df <- imap_dfr(cohorts, \(df, nm) tibble(
    cohort=nm, sensitivity=compute_sensitivity(df, CFG$focal_epitope))) %>%
    mutate(cohort=factor(cohort, levels=names(cohorts)))
  p30 <- ggplot(sens_df, aes(cohort, sensitivity, fill=cohort)) +
    geom_col(width=0.65) +
    geom_text(aes(label=ifelse(is.na(sensitivity),"NA",percent(sensitivity, accuracy=0.1))),
              vjust=-0.4, size=3.5) +
    scale_fill_manual(values=c("Full SPACE-TCR"="#e31a1c","Shared only"="#9c27b0","SPACE-TCR only"="#fb6a4a"), guide="none") +
    scale_y_continuous(labels=percent_format(), limits=c(0,1.05), expand=expansion(mult=c(0,0.05))) +
    labs(title=sprintf("Focal-epitope sensitivity (%s) by cohort", CFG$focal_epitope),
         x=NULL, y="Sensitivity")
  save_plot("fig30_sensitivity_by_cohort.png", p30, 9, 6)
  cat("Stage 6b complete.\n")
  invisible(sens_df)
}


# Purity of an arbitrary label (epitope or a gene call) per membership row.
compute_purity_generic <- function(df, label_col = "epitope",
                                   min_size = METRIC$min_cluster_size,
                                   hi_thresh = METRIC$high_purity_cut) {
  stopifnot(label_col %in% names(df), "cluster_id" %in% names(df))
  d <- df %>% filter(!is.na(.data[[label_col]]), .data[[label_col]] != "")
  ct <- if (nrow(d)) cluster_table(.as_members(d), unit = ".member",
                                   min_size = min_size, label = label_col) else tibble()
  if (!nrow(ct))
    return(tibble(label = label_col, n_clusters = 0L, n_clustered = 0L,
                  purity = NA_real_, perc_clusters_gt09 = NA_real_,
                  perc_tcrs_gt09 = NA_real_, mean_labels_per_cluster = NA_real_))
  tibble(label = label_col, n_clusters = nrow(ct), n_clustered = sum(ct$size),
         purity = jm_purity(ct),
         perc_clusters_gt09 = mean(ct$purity >= hi_thresh),
         perc_tcrs_gt09 = sum(ct$size[ct$purity >= hi_thresh]) / sum(ct$size),
         mean_labels_per_cluster = mean(ct$n_labels))
}


sptce_gene_purity <- function(s5, min_size = 2) {
  gl <- .load_gene_labels()
  if (is.null(gl)) return(invisible(NULL))
  combined <- bind_rows(
    s5$sptce_clusters %>% mutate(method = as.character(method)) %>%
      left_join(gl, by = "cdr3b") %>%
      select(method, cluster_id, cdr3b, epitope, v_gene_beta, j_gene_beta),
    s5$published_clusters %>% mutate(method = as.character(method)) %>%
      left_join(gl, by = "cdr3b") %>%
      select(method, cluster_id, cdr3b, epitope, v_gene_beta, j_gene_beta))

  purity_tbl <- combined %>% group_by(method) %>%
    group_modify(\(df, key) bind_rows(
      compute_purity_generic(df, "epitope",     min_size),
      compute_purity_generic(df, "v_gene_beta", min_size),
      compute_purity_generic(df, "j_gene_beta", min_size))) %>%
    ungroup() %>%
    mutate(label = recode(label, epitope="Epitope", v_gene_beta="V-gene (b)",
                          j_gene_beta="J-gene (b)"),
           label = factor(label, c("Epitope","V-gene (b)","J-gene (b)")))
  .ext_write_csv(purity_tbl, "gene_purity_by_method.csv")

  # Explicit ordering: robust to all-NA methods that break fct_reorder in
  # forcats 1.0.
  meth_order <- purity_tbl %>% filter(label == "Epitope") %>%
    arrange(desc(replace_na(perc_clusters_gt09, -Inf))) %>% pull(method)
  meth_order <- c(meth_order, setdiff(unique(purity_tbl$method), meth_order))
  purity_tbl <- purity_tbl %>% mutate(method = factor(method, levels = meth_order))

  p <- ggplot(purity_tbl, aes(x = method, y = perc_clusters_gt09, fill = label)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = scales::percent(perc_clusters_gt09, accuracy = 1)),
              position = position_dodge(width = 0.75), vjust = -0.3, size = 2.6, na.rm = TRUE) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.05),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c("Epitope"="#2C7FB8","V-gene (b)"="#7FCDBB","J-gene (b)"="#EDF8B1")) +
    labs(title = "Cluster purity by label type",
         subtitle = "% clusters at purity >= 0.9 for epitope, V-gene, and J-gene",
         x = NULL, y = "% clusters at purity >= 0.9", fill = "Label") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major.x = element_blank())
  .ext_save("fig_gene_purity_by_method.png", p, w = 11, h = 6)
  invisible(purity_tbl)
}

message("SPACE-TCR pipeline loaded. ",
        "dihedral gate = ", CFG$focal_dihedral, " rad; focal RMSD = ",
        CFG$focal_threshold, " A.")
# End of file.
