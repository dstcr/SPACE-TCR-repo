# =============================================================================
# install/r_packages.R — install and verify the R side, on any platform
#
#   Rscript install/r_packages.R              # core packages only (analysis)
#   Rscript install/r_packages.R --figures    # core + the plotting stack
#   Rscript install/r_packages.R --check      # verify only, install nothing
#
# or from inside R / RStudio, with the repo as the working directory:
#
#   source("install/r_packages.R")            # core
#   SPTCE_FIGURES <- TRUE; source("install/r_packages.R")
#
# Works on Windows, macOS and Linux, with no shell needed. Packages already
# installed are left alone, so it can be re-run safely.
#
#   CORE (7)      what sptce_metrics.R, validate_outputs.R and the tests need.
#                 Pure R; installs in a couple of minutes.
#   FIGURES (13)  what sptce_pipeline.R needs to draw. igraph, ggraph and
#                 tidygraph build against system libraries; see the notes
#                 printed at the end.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
want_figures <- isTRUE(get0("SPTCE_FIGURES", ifnotfound = FALSE)) ||
  "--figures" %in% args
check_only <- "--check" %in% args

CORE <- c("readr", "dplyr", "tidyr", "purrr", "stringr", "tibble", "stringdist")
FIGURES <- c("ggplot2", "patchwork", "scales", "forcats", "reshape2", "igraph",
             "ggraph", "tidygraph", "cowplot", "ggseqlogo", "ggforce",
             "ggrepel", "kableExtra")
# Only needed by ml/tests/test_metrics_match_r.py, which shuttles results between R and
# Python. Installed with --figures since it is tiny and avoids a later surprise.
EXTRA <- c("jsonlite")

# Versions this was verified against in October 2026, for reference when a
# result will not reproduce. Not enforced: pinning R packages hard tends to
# break more installs than it fixes.
VERIFIED <- c(
  readr = "2.1.5", dplyr = "1.1.4", tidyr = "1.3.1", purrr = "1.0.2",
  stringr = "1.5.1", tibble = "3.2.1", stringdist = "0.9.12",
  ggplot2 = "3.4.4", patchwork = "1.2.0", scales = "1.3.0", forcats = "1.0.0",
  reshape2 = "1.4.4", igraph = "1.6.0", ggraph = "2.1.0", tidygraph = "1.3.0",
  cowplot = "1.1.3", ggseqlogo = "0.1", ggforce = "0.4.1", ggrepel = "0.9.5",
  kableExtra = "1.4.0", jsonlite = "1.8.8"
)

ok  <- function(m, x = "") cat(sprintf("  [ok]   %-30s %s\n", m, x))
bad <- function(m, x = "") cat(sprintf("  [FAIL] %-30s %s\n", m, x))
nb  <- function(m, x = "") cat(sprintf("  [note] %-30s %s\n", m, x))

cat("\n================================================================\n")
cat("  SPTCE R packages\n")
cat("================================================================\n")
cat("  R:        ", R.version.string, "\n")
cat("  platform: ", R.version$platform, "\n")
cat("  library:  ", .libPaths()[1], "\n")
cat("  tier:     ", if (want_figures) "core + figures" else "core only", "\n\n")

if (getRversion() < "4.1.0") {
  cat("  R 4.1 or newer is required: the pipeline uses the \\(x) lambda\n")
  cat("  shorthand, which older parsers reject.\n\n")
  if (!interactive()) quit(status = 1)
}

# A repo must be set or install.packages() prompts, which breaks Rscript.
repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || !nzchar(repos["CRAN"]) ||
    identical(unname(repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  nb("CRAN mirror", "defaulted to cloud.r-project.org")
}

wanted <- c(CORE, if (want_figures) c(FIGURES, EXTRA))
have <- vapply(wanted, requireNamespace, logical(1), quietly = TRUE)

if (!check_only && any(!have)) {
  todo <- wanted[!have]
  cat("  installing", length(todo), "package(s):", paste(todo, collapse = " "), "\n\n")
  # One at a time, so a single failure does not abort the rest and the log says
  # which package died.
  for (p in todo) {
    cat("  ->", p, "\n")
    try(install.packages(p, quiet = TRUE), silent = TRUE)
  }
  cat("\n")
  have <- vapply(wanted, requireNamespace, logical(1), quietly = TRUE)
}

cat("--- core ---\n")
n_fail <- 0L
for (p in CORE) {
  if (requireNamespace(p, quietly = TRUE)) {
    v <- as.character(utils::packageVersion(p))
    tag <- if (!is.na(VERIFIED[p]) && v != VERIFIED[p])
      sprintf("(verified against %s)", VERIFIED[p]) else ""
    ok(p, paste(v, tag))
  } else {
    bad(p, "not installed"); n_fail <- n_fail + 1L
  }
}

if (want_figures) {
  cat("\n--- figures ---\n")
  for (p in c(FIGURES, EXTRA)) {
    if (requireNamespace(p, quietly = TRUE)) {
      ok(p, as.character(utils::packageVersion(p)))
    } else {
      bad(p, "not installed"); n_fail <- n_fail + 1L
    }
  }
}

# ---- does the pipeline actually load? ---------------------------------------
cat("\n--- loading the pipeline ---\n")
if (all(vapply(CORE, requireNamespace, logical(1), quietly = TRUE))) {
  for (f in c("sptce_metrics.R", "sptce_pipeline.R")) {
    if (!file.exists(f)) { nb(f, "not found (wrong working directory?)"); next }
    e <- tryCatch({
      env <- new.env(); suppressMessages(sys.source(f, envir = env))
      n <- sum(vapply(ls(env, all.names = TRUE),
                      function(n) is.function(get(n, envir = env)), logical(1)))
      ok(f, paste(n, "functions"))
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(e)) { bad(f, e); n_fail <- n_fail + 1L }
  }
} else {
  nb("skipped", "core packages missing")
}

cat("\n================================================================\n")
if (n_fail == 0L) {
  cat("  R side ready.\n")
  if (!want_figures)
    cat("  For the figures, re-run with --figures (13 more packages).\n")
  cat("  Next:  Rscript check_install.R\n")
} else {
  cat("  ", n_fail, " package(s) or file(s) failed.\n", sep = "")
  cat("\n  Common causes:\n")
  cat("   * igraph / ggraph / tidygraph compile against system libraries.\n")
  cat("     Debian/Ubuntu:  sudo apt-get install libglpk-dev libxml2-dev \\\n")
  cat("                        libudunits2-dev libgdal-dev libfontconfig1-dev\n")
  cat("     macOS:          brew install glpk libxml2 udunits gdal\n")
  cat("     Windows:        install Rtools matching your R version\n")
  cat("   * tidygraph >= 1.3.1 requires igraph >= 2.0. With an older igraph,\n")
  cat("     install the matching pair:\n")
  cat("       remotes::install_version(\"tidygraph\", \"1.3.0\")\n")
  cat("       remotes::install_version(\"ggraph\",    \"2.1.0\")\n")
  cat("   * no write access to the library: set one you own with\n")
  cat("       dir.create(\"~/R/lib\", recursive = TRUE); .libPaths(\"~/R/lib\")\n")
}
cat("================================================================\n\n")
if (!interactive() && n_fail > 0L) quit(status = 1)
invisible(n_fail == 0L)
