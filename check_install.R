# =============================================================================
# check_install.R — verify the R side of a SPTCE install, on any platform
#
#   Rscript check_install.R
#
# or from inside R / RStudio, with the repo as the working directory:
#
#   source("check_install.R")
#
# This is the portable equivalent of check_install.sh's R sections. Use it on
# Windows, where the shell script needs Git Bash or WSL. It checks the R
# toolchain, the packages, that every R file parses, and runs all three R test
# suites. It does not check the Python/SPACE2TCR side — that only matters on the
# machine that runs Stage C, which is the cluster. Use check_install.sh there.
#
# Exits non-zero if anything failed, so it can gate a script.
# =============================================================================

pass <- 0L; fail <- 0L; skip <- 0L
ok   <- function(m, x = "") { pass <<- pass + 1L; cat(sprintf("  [ok]   %-48s %s\n", m, x)) }
bad  <- function(m, x = "") { fail <<- fail + 1L; cat(sprintf("  [FAIL] %-48s %s\n", m, x)) }
skp  <- function(m, x = "") { skip <<- skip + 1L; cat(sprintf("  [skip] %-48s %s\n", m, x)) }
hdr  <- function(m) cat(sprintf("\n=== %s ===\n", m))

cat("\n================================================================\n")
cat("  SPTCE install check (R side)\n")
cat("================================================================\n")

# ---- 0. are we in the right place? -----------------------------------------
hdr("0. repository")
if (!file.exists("sptce_metrics.R")) {
  cat("\n  sptce_metrics.R not found in the working directory.\n")
  cat("  Change to the repository root first:\n")
  cat("    setwd(\"path/to/SPACE-TCR\")   # in R\n")
  cat("    cd path/to/SPACE-TCR          # in a terminal\n\n")
  quit(status = 1)
}
ok("working directory is the repo root", getwd())
for (f in c("sptce_metrics.R", "sptce_pipeline.R", "config.sh",
            "tests/test_metrics.R", "tests/test_e2e.R",
            "tests/test_pipeline.R", "make_figures.R")) {
  if (file.exists(f)) ok(f) else bad(paste(f, "is missing"))
}

# ---- 1. toolchain ----------------------------------------------------------
hdr("1. R toolchain")
ok("R version", R.version.string)
if (getRversion() < "4.1.0")
  bad("R is older than 4.1", "the lambda shorthand \\(x) will not parse")
ok("platform", paste(R.version$platform, "/", Sys.info()[["sysname"]]))

core <- c("readr", "dplyr", "tidyr", "purrr", "stringr", "tibble", "stringdist")
have_core <- vapply(core, requireNamespace, logical(1), quietly = TRUE)
if (all(have_core)) {
  ok("all 7 core packages present")
} else {
  missing <- core[!have_core]
  bad("missing core packages", paste(missing, collapse = " "))
  cat(sprintf('         install.packages(c(%s))\n',
              paste0('"', missing, '"', collapse = ", ")))
}

figs <- c("ggplot2", "patchwork", "scales", "forcats", "reshape2", "igraph",
          "ggraph", "tidygraph", "cowplot", "ggseqlogo", "ggforce", "ggrepel",
          "kableExtra")
have_figs <- vapply(figs, requireNamespace, logical(1), quietly = TRUE)
if (all(have_figs)) {
  ok("all 13 figure packages present", "make_figures.R will run")
} else {
  skp("figure packages absent", paste(figs[!have_figs], collapse = " "))
}

# ---- 2. everything parses --------------------------------------------------
hdr("2. R files parse")
rfiles <- c("sptce_metrics.R", "sptce_pipeline.R", "make_figures.R",
            Sys.glob("tests/*.R"))
for (f in rfiles) {
  e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (is.null(e)) ok(paste("parses:", f)) else bad(paste("parse error:", f), e)
}

# ---- 3. the pipeline is complete -------------------------------------------
hdr("3. sptce_pipeline.R completeness")
if (!all(have_core)) {
  skp("pipeline checks", "core packages missing")
} else {
  env <- new.env()
  e <- tryCatch({ suppressMessages(sys.source("sptce_pipeline.R", envir = env)); NULL },
                error = function(e) conditionMessage(e))
  if (!is.null(e)) {
    bad("sptce_pipeline.R sources", e)
  } else {
    fns <- ls(env, all.names = TRUE)
    fns <- fns[vapply(fns, function(n) is.function(get(n, envir = env)), logical(1))]
    ok("sources cleanly", paste(length(fns), "functions"))
    need <- c("main", "run_benchmark_a", "run_benchmark_b", "run_stage5",
              "jouannet_10x_ladder", ".size_col", ".sptce_tag",
              ".focal_universe", ".assign_overlay_factor",
              "run_stage5d_cohort_accounting")
    miss <- setdiff(need, fns)
    if (length(miss)) {
      bad("critical functions missing", paste(miss, collapse = " "))
    } else {
      ok("all critical functions defined")
    }
    # the per-mode cohort accounting must be the one in scope
    b <- paste(deparse(body(get("run_stage5d_cohort_accounting", envir = env))),
               collapse = " ")
    if (grepl("modes <- unique", b)) {
      ok("cohort accounting is the per-mode version")
    } else {
      bad("cohort accounting is the POOLED version", "re-pull the repository; per-mode accounting is required")
    }
    CFG <- get("CFG", envir = env)
    ok("dihedral gate", paste(CFG$focal_dihedral, "rad"))
    ok("focal RMSD threshold", paste(CFG$focal_threshold, "A"))
    if (!(CFG$focal_dihedral %in% CFG$dihedral_modes))
      bad("focal_dihedral is not in dihedral_modes",
          paste(CFG$dihedral_modes, collapse = ", "))
  }
}

# ---- 4. test suites --------------------------------------------------------
hdr("4. test suites")
run_suite <- function(path, pattern, label) {
  if (!file.exists(path)) { bad(paste(path, "missing")); return(invisible()) }
  out <- suppressWarnings(
    system2(file.path(R.home("bin"), "Rscript"), shQuote(path),
            stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  hit <- grep(pattern, out, value = TRUE)
  if (is.null(status) || status == 0) {
    ok(label, if (length(hit)) tail(hit, 1) else "passed")
  } else {
    bad(label, if (length(hit)) tail(hit, 1) else "see output below")
    cat(paste0("         ", grep("FAIL|Error", out, value = TRUE)[1:3], "\n"))
  }
}
if (!all(have_core)) {
  skp("test suites", "core packages missing")
} else {
  run_suite("tests/test_metrics.R",       "passed, .* failed", "metric unit tests")
  run_suite("tests/test_pipeline.R","passed, .* failed", "pipeline checks")
  run_suite("tests/test_dihedral_distance.R", "passed, .* failed", "dihedral distance (R and Python)")
  run_suite("tests/test_e2e.R",           "e2e: .* failures",  "end-to-end tests")
}

# ---- 5. data paths (informational) -----------------------------------------
hdr("5. data paths in CFG (informational)")
if (all(have_core) && exists("env") && !is.null(get0("CFG", envir = env))) {
  CFG <- get("CFG", envir = env)
  chk <- function(label, p) {
    if (is.null(p) || length(p) != 1 || is.na(p) || !nzchar(p)) {
      skp(label, "unset")
    } else if (file.exists(p)) {
      ok(label, p)
    } else {
      skp(label, paste("not on this machine:", p))
    }
  }
  chk("labelled input",   CFG$jouannet_labelled)
  chk("SPACE-TCR output", CFG$jouannet_sptce_out)
  chk("Jouannet repo",    CFG$jouannet_repo)
  chk("STCRDab CSV",      CFG$stcrdab_parsed_csv)
  cat("\n  Skips here are expected on a machine without the data. Edit the CFG\n")
  cat("  block at the top of sptce_pipeline.R.\n")
} else skp("CFG paths", "pipeline did not load")

# ---- summary ---------------------------------------------------------------
cat("\n================================================================\n")
cat(sprintf("  %d passed, %d failed, %d skipped\n", pass, fail, skip))
if (fail == 0L) {
  cat("  The R side of the install is good.\n")
  if (skip > 0L)
    cat("  Skips are optional components or absent data — see INSTALL.md.\n")
} else {
  cat("  Fix the failures above. INSTALL.md has the install steps.\n")
}
cat("================================================================\n\n")
if (!interactive() && fail > 0L) quit(status = 1)
invisible(fail == 0L)
