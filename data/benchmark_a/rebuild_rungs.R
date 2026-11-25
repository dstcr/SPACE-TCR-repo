# Rebuild the six per-rung ladder tables from the last one plus the membership
# column, for a tool that wants them as separate files (verify_ladder.sh, the
# published methods' own input format, anything outside this repository).
#
#   Rscript data/benchmark_a/rebuild_rungs.R [out_dir]
#
# Default out_dir is a new build/ folder beside this script; nothing is written
# into data/benchmark_a/ itself.
suppressMessages({ library(readr); library(dplyr) })

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (!length(here) || is.na(here)) here <- "data/benchmark_a"
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else file.path(here, "build")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

full <- read_tsv(file.path(here, "donor1_subset6_sptce.tsv"),  show_col_types = FALSE)
lab  <- read_tsv(file.path(here, "donor1_subset6_labels.tsv"), show_col_types = FALSE)
memb <- read_tsv(file.path(here, "ladder_membership.tsv"),     show_col_types = FALSE)

stopifnot(setequal(full$id, memb$id))
expected <- c(3056L, 8089L, 21413L, 23529L, 23529L, 23529L)

for (i in 1:6) {
  keep <- memb$id[memb$first_rung <= i]
  s <- full %>% filter(id %in% keep)
  l <- lab  %>% filter(id %in% keep)
  if (nrow(s) != expected[i])
    stop("rung ", i, ": rebuilt ", nrow(s), " rows, expected ", expected[i],
         ". The membership file and the full table have drifted apart.")
  write_tsv(s, file.path(out, sprintf("donor1_subset%d_sptce.tsv",  i)))
  write_tsv(l, file.path(out, sprintf("donor1_subset%d_labels.tsv", i)))
  cat(sprintf("  rung %d: %6d receptors, %5d labelled\n", i, nrow(s),
              sum(!is.na(l$epitope) & l$epitope != "NA" & nzchar(l$epitope))))
}
cat("\n  wrote 12 files to ", out, "\n", sep = "")
