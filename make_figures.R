# =============================================================================
# make_figures.R — render every figure for the paper in one run
#
#   Rscript make_figures.R                      # uses the paths in sptce_pipeline.R's CFG
#   Rscript make_figures.R /path/to/figures     # write somewhere else
#
# Run from the repository root. It sources sptce_pipeline.R, checks which
# inputs are present, renders what they support, and writes a manifest of what
# was produced, what was skipped, and why. Each figure is run separately, so a
# failure is recorded and the run continues.
#
# INPUTS
# ------
#   CFG$jouannet_labelled    labelled SPACE-TCR input TSV      small
#   CFG$jouannet_sptce_out   distance_tables/*.tsv.gz          large
#                            dihedrals.tsv.gz                  medium
#   CFG$jouannet_repo        the Jouannet benchmark repo       ~382 MB
#   CFG$stcrdab_parsed_csv   igblast_tcr_parsed.csv            small (Fig 8 only)
#   CFG$tenx_root/_out       10X subsets + published outputs   (Benchmark A only)
#
# Figures are written at 300 dpi by default; set SPTCE_FIG_DPI to change it.
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
dpi     <- as.numeric(Sys.getenv("SPTCE_FIG_DPI", "300"))
t_start <- Sys.time()

cat("\n================================================================\n")
cat("  SPACE-TCR figure rendering\n")
cat("================================================================\n")

if (!file.exists("sptce_pipeline.R"))
  stop("run this from the repository root (sptce_pipeline.R not found)", call. = FALSE)

suppressMessages(source("sptce_pipeline.R"))

out_dir <- if (length(args) >= 1) args[1] else CFG$analysis_dir
# out_dir is created only once the inputs have been checked, just before
# anything is written.
cat("  output:  ", normalizePath(out_dir, mustWork = FALSE), "\n")
cat("  gate:    ", CFG$focal_dihedral, "rad     focal RMSD:", CFG$focal_threshold, "A\n")
cat("  dpi:     ", dpi, "\n\n")

# ---- which inputs are present ------------------------------------------------
have <- function(p) !is.null(p) && length(p) == 1 && !is.na(p) && nzchar(p) && file.exists(p)
dist_dir  <- file.path(CFG$jouannet_sptce_out, "distance_tables")
dih_file  <- file.path(CFG$jouannet_sptce_out, "dihedrals.tsv.gz")
n_dist    <- if (dir.exists(dist_dir)) length(list.files(dist_dir, "\\.tsv\\.gz$")) else 0L

inputs <- data.frame(
  input = c("labelled input", "distance tables", "dihedrals",
            "Jouannet repo", "STCRDab parsed CSV", "10X root", "10X outputs"),
  path  = c(CFG$jouannet_labelled, dist_dir, dih_file, CFG$jouannet_repo,
            CFG$stcrdab_parsed_csv, CFG$tenx_root %||% NA, CFG$tenx_out %||% NA),
  present = c(have(CFG$jouannet_labelled), n_dist > 0, have(dih_file),
              have(CFG$jouannet_repo), have(CFG$stcrdab_parsed_csv),
              have(CFG$tenx_root %||% NA), have(CFG$tenx_out %||% NA)),
  stringsAsFactors = FALSE)
cat("  inputs:\n")
for (i in seq_len(nrow(inputs)))
  cat(sprintf("    [%s] %-20s %s\n", if (inputs$present[i]) "ok" else "--",
              inputs$input[i],
              if (inputs$input[i] == "distance tables" && n_dist > 0)
                paste0(n_dist, " files") else inputs$path[i]))
cat("\n")

if (!inputs$present[inputs$input == "labelled input"] ||
    !inputs$present[inputs$input == "distance tables"]) {
  cat("  The labelled input and the distance tables are required for every\n")
  cat("  SPACE-TCR figure. Set CFG$jouannet_labelled and CFG$jouannet_sptce_out\n")
  cat("  in sptce_pipeline.R, or run on the machine where that output lives.\n")
  cat("  Nothing was rendered.\n\n")
  quit(status = 1)
}

# Inputs are present, so an output directory is now warranted.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) {
  cat("  Could not create the output directory:", out_dir, "\n")
  cat("  Pass a writable path:  Rscript make_figures.R /path/to/figures\n\n")
  quit(status = 1)
}

# ---- bookkeeping -------------------------------------------------------------
log <- list()
step <- function(label, needs = TRUE, expr) {
  if (!isTRUE(needs)) {
    cat(sprintf("  [skip] %-42s inputs absent\n", label))
    log[[length(log) + 1]] <<- data.frame(step = label, status = "skipped",
      detail = "required inputs absent", seconds = NA_real_, stringsAsFactors = FALSE)
    return(invisible(NULL))
  }
  t0 <- Sys.time()
  # withCallingHandlers rather than tryCatch for warnings, so that a warning
  # is recorded and the stage carries on instead of being interrupted.
  res <- tryCatch(
    withCallingHandlers({ force(expr); "ok" },
                        warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (identical(res, "ok")) {
    cat(sprintf("  [ok]   %-42s %6.1fs\n", label, el))
    log[[length(log) + 1]] <<- data.frame(step = label, status = "ok", detail = "",
      seconds = el, stringsAsFactors = FALSE)
  } else {
    cat(sprintf("  [FAIL] %-42s %6.1fs  %s\n", label, el, substr(res, 1, 70)))
    log[[length(log) + 1]] <<- data.frame(step = label, status = "failed",
      detail = res, seconds = el, stringsAsFactors = FALSE)
  }
  invisible(NULL)
}

# ---- Benchmark B: the curated cohort ----------------------------------------
cat("--- Benchmark B (curated CD8 cohort) ---\n")
s5 <- NULL
step("stage 5: clusters, metrics, bubble + Tables 1-3", TRUE, { s5 <<- run_stage5() })

if (is.null(s5)) {
  cat("\n  Stage 5 produced nothing; the figures downstream of it cannot run.\n")
} else {
  step("stage 5b: networks and cluster overlays", TRUE, run_stage5b_networks(s5))
  step("stage 5c: extra figures",                 TRUE, run_stage5c_extras(s5))
  step("stage 5d: cohort accounting (per mode)",  TRUE, run_stage5d_cohort_accounting(s5))
  step("stage 6: gene composition",               TRUE, run_stage6(s5))
  step("stage 6b: cohort differences",            TRUE, run_stage6b_cohort_diff(s5))
  step("paired CDR3 sequence logos (Fig 7)",      TRUE, sptce_pairing_logos(s5))
  step("stage 11: extensions",                    TRUE, run_stage11_extensions(s5))
}

# ---- Stage 7: training-set leakage (Fig 8) ----------------------------------
cat("\n--- Stage 7 (training-set leakage) ---\n")
step("STCRDab overlap and held-out validation (Fig 8)",
     have(CFG$stcrdab_parsed_csv), run_stage7_training_leakage())

# ---- Benchmark A: the noise ladder (Fig 9) ----------------------------------
cat("\n--- Benchmark A (10X noise ladder) ---\n")
have_a <- have(CFG$tenx_root %||% NA) && have(CFG$tenx_out %||% NA)
step("noise ladder, 6 rungs (Fig 9)", have_a, run_benchmark_a())

# ---- the paper figures -----------------------------------------------------
# The stages above call the pipeline's figure functions. Figures 1 to 6 and
# Supplementary S1 to S7, as they appear in the paper, are drawn by the rebuild
# scripts, last, into a rebuild/ subfolder. They read only tables in this
# repository, so they run even when the stages above were skipped.
cat("\n--- paper figures ---\n")
rebuild_dir <- file.path(out_dir, "rebuild")
for (sc in c("scripts/fig01_rebuild.R", "scripts/fig_rebuild_all.R",
             "scripts/fig_rebuild_supp.R")) {
  step(paste("rebuild:", basename(sc)), file.exists(sc), {
    rc <- system2("Rscript", c(sc, shQuote(rebuild_dir)),
                  stdout = TRUE, stderr = TRUE)
    st <- attr(rc, "status")
    if (!is.null(st) && st != 0) stop(paste(utils::tail(rc, 12), collapse = "\n"))
    invisible(TRUE)
  })
}

# Supplementary Figures S8 and S9 are drawn in Python (matplotlib, pandas). S9
# reads the committed constant-alpha run as well as the cohort tables.
s9_in <- c("results/tables/extensions/s9_constant_alpha",
           "results/tables/extensions/s9_constant_alpha_dihedrals.tsv.gz")
for (sc in c("scripts/fig_s8_tcrvdb.py", "scripts/fig_s9_constant_alpha.py")) {
  have <- file.exists(sc) && nzchar(Sys.which("python3")) &&
          (!grepl("s9", sc) || all(file.exists(s9_in)))
  step(paste("rebuild:", basename(sc)), have, {
    rc <- system2("python3", c(sc, shQuote(rebuild_dir)),
                  stdout = TRUE, stderr = TRUE)
    st <- attr(rc, "status")
    if (!is.null(st) && st != 0) stop(paste(utils::tail(rc, 12), collapse = "\n"))
    invisible(TRUE)
  })
}

# ---- manifest ----------------------------------------------------------------
cat("\n--- manifest ---\n")
man <- if (length(log)) do.call(rbind, log) else
  data.frame(step = character(), status = character(), detail = character(),
             seconds = numeric(), stringsAsFactors = FALSE)
man$gate <- CFG$focal_dihedral
man$focal_threshold <- CFG$focal_threshold
man$rendered_at <- format(t_start, "%Y-%m-%d %H:%M:%S")

figs <- list.files(out_dir, "\\.(png|pdf|svg|tiff)$", full.names = TRUE,
                   recursive = TRUE)
# Count only files written by this run.
fresh <- figs[file.info(figs)$mtime >= t_start]
cat(sprintf("  %d figure file(s) written to %s\n", length(fresh), out_dir))
if (length(fresh)) {
  sz <- file.info(fresh)$size
  for (i in order(basename(fresh)))
    cat(sprintf("    %-56s %7.0f KB\n", basename(fresh[i]), sz[i] / 1024))
}

man_path <- file.path(out_dir, sprintf("figure_manifest_r%s_d%s.csv",
                                       CFG$focal_threshold, CFG$focal_dihedral))
write.csv(man, man_path, row.names = FALSE)
files_path <- file.path(out_dir, sprintf("figure_files_r%s_d%s.csv",
                                         CFG$focal_threshold, CFG$focal_dihedral))
if (length(fresh))
  write.csv(data.frame(file = basename(fresh), bytes = file.info(fresh)$size,
                       stringsAsFactors = FALSE), files_path, row.names = FALSE)

ok <- sum(man$status == "ok"); fl <- sum(man$status == "failed")
sk <- sum(man$status == "skipped")
cat("\n================================================================\n")
cat(sprintf("  %d ok, %d failed, %d skipped in %.1f min\n", ok, fl, sk,
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("  manifest:", basename(man_path), "\n")
if (fl > 0) {
  cat("\n  failures:\n")
  for (i in which(man$status == "failed"))
    cat(sprintf("    %-44s %s\n", man$step[i], substr(man$detail[i], 1, 90)))
}
if (sk > 0) cat("  skipped steps needed inputs that were not present (see above).\n")
cat("================================================================\n\n")
if (fl > 0) quit(status = 1)
