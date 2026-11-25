# =============================================================================
# run_benchmark_a.R — render Figure 9, the 10X donor-1 noise ladder
#
#   Rscript run_benchmark_a.R <jouannet_repo> <sptce_root> [out_dir]
#
#   <jouannet_repo>  the TCR_Unsupervised_Benchmark clone
#                    (third_party/fetch_jouannet.sh gets it)
#   <sptce_root>     directory holding subsetN_clusters/ and subsetN_dihedrals/
#                    for at least ONE rung -- one structural run serves all six
#   [out_dir]        where figures go (default: CFG$analysis_dir)
#
# The two 10X tables are read from the environment if set, since they are not
# redistributable and live wherever you downloaded them:
#
#   SPTCE_BINMAT=.../vdj_v1_hs_aggregated_donor1_binarized_matrix.csv
#   SPTCE_CONTIG=.../vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv
#
# SPTCE_VERSIONS=a,b,c reads all three replicate runs. The default is "a",
# matching run_benchmark_a(), and that default is deliberate: a/b/c are
# replicate RUNS, and HD/LD/iSMART/GIANA/TCRMatch are deterministic, so their
# a/b/c files are byte-identical and averaging them would fake precision. Use
# all three only when reading the stochastic tools (clusTCR, DeepTCR), and say
# in the caption which methods the spread is real for. See docs/BENCHMARK_A.md.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("usage: Rscript run_benchmark_a.R <jouannet_repo> <sptce_root> [out_dir]\n")
  quit(status = 2)
}

if (!file.exists("sptce_pipeline.R"))
  stop("run this from the repository root", call. = FALSE)

suppressMessages(source("sptce_pipeline.R"))

J   <- normalizePath(args[1], mustWork = FALSE)
SR  <- normalizePath(args[2], mustWork = FALSE)
OUT <- if (length(args) >= 3) args[3] else CFG$analysis_dir

CFG$jouannet_repo <- J
CFG$tenx_root     <- file.path(J, "Data", "10X", "Donor1")
CFG$tenx_out      <- file.path(J, "Data", "10X_Output_methods")
CFG$sptce_root    <- SR
CFG$analysis_dir  <- OUT
# .paired_betas() reads <subset_dir>/_pre_ladder_fix/donor1_subsetN_sptce.tsv,
# if present, in preference to CFG$master_file, because the structural run's
# TCR ids are keyed to that subset definition. Here both are under sptce_root.
CFG$subset_dir    <- Sys.getenv("SPTCE_SUBSET_DIR", unset = SR)

bm <- Sys.getenv("SPTCE_BINMAT", unset = NA)
ct <- Sys.getenv("SPTCE_CONTIG", unset = NA)
if (!is.na(bm) && nzchar(bm)) CFG$binmat_file <- bm
if (!is.na(ct) && nzchar(ct)) CFG$contig_file <- ct

cat("\n================================================================\n")
cat("  Benchmark A — 10X donor-1 noise ladder (Figure 9)\n")
cat("================================================================\n")
ck <- function(lbl, p, dir = FALSE) {
  ok <- !is.null(p) && length(p) == 1 && !is.na(p) && nzchar(p) &&
        (if (dir) dir.exists(p) else file.exists(p))
  cat(sprintf("  [%s] %-16s %s\n", if (ok) "ok" else "--", lbl, p))
  ok
}
o1 <- ck("tenx_root",   CFG$tenx_root,   TRUE)
o2 <- ck("tenx_out",    CFG$tenx_out,    TRUE)
o3 <- ck("sptce_root",  CFG$sptce_root,  TRUE)
ok_sub <- ck("subset_dir",  CFG$subset_dir,  TRUE)
o4 <- ck("binmat_file", CFG$binmat_file)
o5 <- ck("contig_file", CFG$contig_file)
cat("\n")

if (!(o1 && o2 && o3 && ok_sub)) {
  cat("  The benchmark repository and one structural run are required.\n")
  cat("  Get the repository with third_party/fetch_jouannet.sh.\n\n")
  quit(status = 1)
}

lab_csv <- file.path(CFG$tenx_root, "labelled_cdr3b.csv")
if (!file.exists(lab_csv) && !(o4 && o5)) {
  cat("  Epitope labels are missing and cannot be reconstructed.\n\n")
  cat("  jouannet_10x_labels() looks for\n")
  cat("    ", lab_csv, "\n")
  cat("  whose write is commented out in the benchmark's 10X_processing.R, so\n")
  cat("  it is normally absent. The fallback rebuilds the labels their way from\n")
  cat("  the raw 10X donor-1 tables, which you must supply:\n\n")
  cat("    export SPTCE_BINMAT=.../vdj_v1_hs_aggregated_donor1_binarized_matrix.csv\n")
  cat("    export SPTCE_CONTIG=.../vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv\n\n")
  cat("  Both are in the 10X Genomics donor-1 release; the benchmark's README\n")
  cat("  points at https://zenodo.org/records/6952657 . Nothing else is missing.\n")
  cat("================================================================\n\n")
  quit(status = 1)
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sptce_preflight(CFG$sptce_root)

t0 <- Sys.time()
vers <- strsplit(Sys.getenv("SPTCE_VERSIONS", unset = "a"), ",")[[1]]
vers <- trimws(vers[nzchar(trimws(vers))])
cat("  replicate runs: ", paste(vers, collapse = ", "),
    if (length(vers) > 1)
      "  (deterministic methods give identical files; spread is real only for clusTCR and DeepTCR)"
    else "", "\n\n", sep = "")
# dih_fig = NA: the field-comparison panel is drawn RMSD-only, as in Figure 9.
# Pass a gate here to draw a gated comparison instead.
res <- run_benchmark_a(versions = vers, dih_fig = NA)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

if (is.null(res)) {
  cat("\n  Benchmark A produced nothing; see the messages above.\n\n")
  quit(status = 1)
}
saveRDS(res, file.path(OUT, "benchmark_a.rds"))

figs <- list.files(OUT, "\\.(png|pdf)$", full.names = TRUE)
figs <- figs[file.info(figs)$mtime >= t0]
cat(sprintf("\n  %d file(s) written to %s in %.1f min\n", length(figs), OUT, el))
for (f in sort(basename(figs))) cat("    ", f, "\n")
cat("================================================================\n\n")
