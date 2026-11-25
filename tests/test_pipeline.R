# Assertions on sptce_pipeline.R.
#
# In R a later definition of the same name silently wins, so these tests pin
# which functions are defined — above all that run_stage5d_cohort_accounting
# accounts per mode rather than pooling operating points.
#
#   Rscript tests/test_pipeline.R

pass <- 0; fail <- 0
chk <- function(label, cond, extra = "") {
  if (isTRUE(cond)) { pass <<- pass + 1; cat(sprintf("  PASS  %-56s %s\n", label, extra)) }
  else              { fail <<- fail + 1; cat(sprintf("  FAIL  %-56s %s\n", label, extra)) }
}

repo <- if (file.exists("sptce_pipeline.R")) "." else ".."
f <- file.path(repo, "sptce_pipeline.R")
chk("sptce_pipeline.R exists", file.exists(f))
if (!file.exists(f)) quit(status = 1)

chk("parses", !inherits(try(parse(f), silent = TRUE), "try-error"))

env <- new.env()
err <- tryCatch({ suppressMessages(sys.source(f, envir = env)); NULL },
                error = function(e) conditionMessage(e))
chk("sources without error", is.null(err), if (is.null(err)) "" else err)
if (!is.null(err)) quit(status = 1)

fns <- ls(env, all.names = TRUE)
fns <- fns[vapply(fns, function(n) is.function(get(n, envir = env)), logical(1))]
chk("defines a full complement of functions", length(fns) > 100, paste(length(fns), "functions"))

cat("\n== stage functions and the helpers they call ==\n")
for (n in c("run_benchmark_a", "jouannet_10x_ladder", ".size_col", ".sptce_tag",
            ".focal_universe", ".assign_overlay_factor", "run_stage5",
            "run_stage5b_networks", "run_stage5c_extras", "run_stage6",
            "run_stage6b_cohort_diff", "run_stage7_training_leakage",
            "load_jouannet_curated", "run_stage1", "sptce_pairing_logos"))
  chk(paste("defined:", n), n %in% fns)

cat("\n== consensus, leakage and layout helpers ==\n")
for (n in c("sptce_consensus_clusters", "sptce_gene_consensus",
            "compute_label_consensus", "sptce_leakage_full_cohort",
            ".spread_layout", ".relax_min_sep", "sptce_preflight"))
  chk(paste("defined:", n), n %in% fns)

# The critical one: cohort accounting must split by mode before tallying.
body_txt <- paste(deparse(body(get("run_stage5d_cohort_accounting", envir = env))),
                  collapse = " ")
chk("cohort accounting is the PER-MODE version",
    grepl("modes <- unique", body_txt),
    "(a pooled tally would count a receptor once per operating point)")
chk("cohort accounting does not pool across modes",
    !grepl("qualifying <- map_dfr\\(CFG\\$jouannet_epitopes", body_txt))

cat("\n== Benchmark C (HuARdb) is present but not wired into main() ==\n")
# These are the only implementation of Benchmark C, so they must be DEFINED.
for (n in c("run_stage10_huardb", "normalise_disease", "DISEASE_METATYPE"))
  chk(paste("defined:", n), n %in% fns || exists(n, envir = env, inherits = FALSE))
# ... but main() does not run Benchmark C, which needs a HuARdb download that
# is not part of this repository.
main_txt <- paste(deparse(body(get("main", envir = env))), collapse = " ")
chk("main() does not call run_stage10_huardb",
    !grepl("run_stage10_huardb", main_txt),
    "(it needs a HuARdb download this repo does not ship)")

cat("\n== the dihedral gate is live, not a dead knob ==\n")
CFG <- get("CFG", envir = env)
chk("CFG$focal_dihedral is set", is.numeric(CFG$focal_dihedral), CFG$focal_dihedral)
chk("CFG$focal_dihedral is in CFG$dihedral_modes",
    CFG$focal_dihedral %in% CFG$dihedral_modes,
    paste(CFG$dihedral_modes, collapse = ", "))
src <- readLines(f, warn = FALSE)
n_reads <- length(grep("CFG$focal_dihedral", src, fixed = TRUE))
chk("CFG$focal_dihedral is read more than once", n_reads > 1, paste(n_reads, "reads"))
# Benchmark A's two panels take their gates from their own settings, not from
# CFG$focal_dihedral (0.5, the value Tables 1 and 2 use for the RMSD+dih rows):
# Figure 9's thresholds panel draws BENCH_A_FIG_DIHEDRAL and its methods panel
# draws no gate. Each setting must be a named constant, not a literal in a call.
dm <- formals(get("fig_paired_methods", envir = env))$dih_keep
chk("fig_paired_methods draws no gate by default (RMSD-only, as published)",
    is.logical(dm) && is.na(dm), paste(deparse(dm), collapse = ""))
dt <- formals(get("fig_paired_thresholds", envir = env))$dih_keep
chk("fig_paired_thresholds takes its gate from a named constant",
    identical(paste(deparse(dt), collapse = ""), "BENCH_A_FIG_DIHEDRAL"),
    paste(deparse(dt), collapse = ""))
bfd <- get0("BENCH_A_FIG_DIHEDRAL", envir = env)
# Figure 9's thresholds panel draws BOTH published gates. It drew only 0.3
# while Tables 1 and 2 report 0.5 for the RMSD+dih rows, so a reader comparing
# them met two different gates with nothing saying so.
chk("BENCH_A_FIG_DIHEDRAL covers both published gates, 0.3 and 0.5",
    is.numeric(bfd) && setequal(round(bfd, 3), c(0.3, 0.5)),
    paste(bfd, collapse = ", "))
chk("the table's gate is among the gates the figure draws",
    any(abs(bfd - CFG$focal_dihedral) < 1e-9),
    paste("table", CFG$focal_dihedral, "in figure set", paste(bfd, collapse = ", ")))
# Two gates must not render as the same line type.
body_ft <- paste(deparse(body(get("fig_paired_thresholds", envir = env))), collapse = " ")
chk("the gate is a factor, not a logical, so each gate gets its own line type",
    grepl("gated = factor", body_ft) && !grepl("gated = !is.na\\(dih\\)", body_ft))
code_only <- src[!grepl("^\\s*#", src)]
chk("no anonymous numeric gate default remains in code",
    length(grep("dih_(keep|fig)\\s*=\\s*0\\.[0-9]", code_only)) == 0)

cat("\n== the six-rung ladder ==\n")
chk("CFG$subsets covers 6 rungs", length(CFG$subsets) == 6,
    paste(CFG$subsets, collapse = ","))
chk("CFG$thresholds is the published sweep",
    identical(as.numeric(CFG$thresholds), c(0.25, 0.50, 0.75, 1.00, 1.25)))
chk("CFG$focal_threshold is 0.75", CFG$focal_threshold == 0.75)

cat(sprintf("\n==== %d passed, %d failed ====\n", pass, fail))
if (fail > 0) quit(status = 1)
