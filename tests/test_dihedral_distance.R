# =============================================================================
# test_dihedral_distance.R — the dihedral distance used for the gate
#
#   Rscript tests/test_dihedral_distance.R
#
# circ_dihedral() in sptce_metrics.R is the definition: mean circular distance
# between two phi/psi vectors, in radians, the shorter padded with pi.
# Stages C and D run on the cluster in Python, so scripts/dihedral_util.py
# carries the same function; this checks both against hand-computed values and
# against each other, and that the clustering step uses it unscaled.
#
# Exit status is the failure count, so it can gate a script.
# =============================================================================

ROOT <- normalizePath(file.path(dirname(sub("--file=", "",
         grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (!length(ROOT) || is.na(ROOT)) ROOT <- "."

fails <- 0L; passes <- 0L
ok <- function(label, pass, detail = "") {
  cat(sprintf("  [%s] %-56s %s\n", if (pass) "ok" else "FAIL", label, detail))
  if (pass) passes <<- passes + 1L else fails <<- fails + 1L
}

# ---- 1. hand-computed values -----------------------------------------------
source(file.path(ROOT, "sptce_metrics.R"), local = TRUE, echo = FALSE)
ok("identical vectors are 0 apart", circ_dihedral(c(1, 2), c(1, 2)) == 0)
ok("wraps: 0.1 and 2pi-0.1 are 0.2 apart",
   abs(circ_dihedral(0.1, 2 * pi - 0.1) - 0.2) < 1e-12)
ok("shorter vector padded with pi",
   abs(circ_dihedral(c(0, 0), 0) - pi / 2) < 1e-12)
ok("an absent loop is pi apart", circ_dihedral(numeric(0), 1) == pi)

# ---- 2. R against the Python used on the cluster ---------------------------

set.seed(11)
cases <- list(
  list(a = c(0.1, 0.2, -3.0),        b = c(0.1, 0.25, 3.0)),   # wraps at +-pi
  list(a = c(pi, -pi),               b = c(-pi, pi)),          # identical mod 2pi
  list(a = runif(52, -pi, pi),       b = runif(52, -pi, pi)),  # the real shape
  list(a = runif(52, -pi, pi),       b = runif(52, -pi, pi)),
  list(a = runif(10, -pi, pi),       b = runif(14, -pi, pi)),  # length mismatch
  list(a = numeric(0),               b = runif(4, -pi, pi))    # one loop absent
)
r_vals <- vapply(cases, function(c_) circ_dihedral(c_$a, c_$b), numeric(1))

payload <- file.path(tempdir(), "dihedral_cases.json")
writeLines(jsonlite::toJSON(lapply(cases, function(c_)
  list(a = as.numeric(c_$a), b = as.numeric(c_$b))), digits = 17), payload)

py <- sprintf('
import json, sys
sys.path.insert(0, %s)
from dihedral_util import dihedral_distance
for c in json.load(open(%s)):
    print(repr(dihedral_distance(c["a"], c["b"])))
', shQuote(file.path(ROOT, "scripts")), shQuote(payload))

out <- suppressWarnings(system2("python3", c("-c", shQuote(py)),
                                stdout = TRUE, stderr = TRUE))
if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
  ok("python implementation callable", FALSE, paste(out, collapse = " | "))
} else {
  p_vals <- suppressWarnings(as.numeric(out))
  if (anyNA(p_vals) || length(p_vals) != length(r_vals)) {
    ok("python returned one value per case", FALSE,
       paste(utils::head(out, 3), collapse = " | "))
  } else {
    d <- max(abs(r_vals - p_vals))
    ok("R circ_dihedral == Python dihedral_distance", d < 1e-12,
       sprintf("max |diff| = %.3g", d))
    for (k in seq_along(cases))
      ok(sprintf("  case %d (len %d vs %d)", k,
                 length(cases[[k]]$a), length(cases[[k]]$b)),
         abs(r_vals[k] - p_vals[k]) < 1e-12,
         sprintf("R %.9f  Py %.9f", r_vals[k], p_vals[k]))
  }
}

# ---- 3. the clustering step must not rescale the gate ----------------------
# A dihedral-mode cutoff is in radians, like every other statement of the gate.
src <- paste(readLines(file.path(ROOT, "scripts/cluster_one_bin.py"),
                       warn = FALSE), collapse = "\n")
ok("cluster_one_bin does not divide the dihedral distance by pi",
   !grepl("dd_beta\\s*\\)\\s*/\\s*2\\.0\\s*\\)\\s*/\\s*np\\.pi", src),
   "gate is in radians")
ok("cluster_one_bin delegates rather than duplicating",
   grepl("_shared_dihedral_distance", src), "single-sourced")

cat(sprintf("\n==== dihedral distance: %d passed, %d failed ====\n", passes, fails))
quit(status = min(fails, 100L))
