# End-to-end test of sptce_metrics.R: synthesise a tiny cohort on disk
# (distance tables + dihedrals, in the formats Stages C and D write), cluster
# it with cluster_at() exactly as the pipeline does, score it, and check the
# results are internally consistent.
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble)})
set.seed(42)
# The repo root is located from this script's path, so it runs from any clone.
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."),
         mustWork = FALSE)
if (is.na(repo_root) || !file.exists(file.path(repo_root, "sptce_metrics.R"))) repo_root <- "."
# Synthetic fixture lives in the session temp dir, so it works on any OS.
root <- file.path(tempdir(), "sptce_e2e_fixture")
unlink(root, recursive = TRUE)
dir.create(file.path(root, "sptce_out", "distance_tables"), recursive = TRUE)

# --- cohort: 60 receptors, 3 epitopes, 2 CDR3-length bins -------------------
N <- 60
epis <- rep(c("GILGFVFTL", "NLVPMVATV", "GLCTLVAML"), length.out = N)
mk_cdr3 <- function(i, e) paste0("CASS", substr(e, 1, 3), sprintf("%02d", i %% 7), "TQYF")
labels <- tibble(
  id            = sprintf("TCR%03d", 1:N),
  cdr3_aa_beta  = mapply(mk_cdr3, 1:N, epis),
  cdr3_aa_alpha = sprintf("CAVR%02dDKIIF", 1:N %% 9),
  epitope       = epis,
  v_b_gene      = rep(c("TRBV19*01","TRBV19*01","TRBV20*01"), length.out = N),
  j_b_gene      = rep(c("TRBJ2-1*01","TRBJ1-2*01"), length.out = N))
write_tsv(labels, file.path(root, "sptce_out", "..", "labelled.tsv"))

# --- distance tables: two bins. Within a bin, same-epitope pairs get small
#     RMSD, cross-epitope pairs get large RMSD, so clusters are epitope-pure at
#     tight thresholds and merge at loose ones.
for (b in 1:2) {
  ids <- labels$id[seq(b, N, by = 2)]
  g <- expand.grid(id_i = ids, id_j = ids, stringsAsFactors = FALSE) %>%
    filter(id_i < id_j)
  e <- setNames(labels$epitope, labels$id)
  g$rmsd <- ifelse(e[g$id_i] == e[g$id_j],
                   runif(nrow(g), 0.2, 0.7), runif(nrow(g), 1.1, 2.0))
  g$status <- "ok"
  write_tsv(g, gzfile(file.path(root, "sptce_out", "distance_tables",
                                sprintf("bin_%d.tsv.gz", b))))
}

# --- dihedrals: 13 CDR3 residues x 2 chains per receptor -------------------
# Angles are epitope-correlated with small jitter, so same-epitope pairs sit
# well inside a 0.3-0.5 rad gate and cross-epitope pairs do not. Purely random
# angles would average ~pi/2 apart and the gate would (correctly) reject every
# pair, leaving nothing to test.
base <- lapply(unique(labels$epitope), function(e) runif(26, -pi, pi))
names(base) <- unique(labels$epitope)
dih <- do.call(rbind, lapply(seq_len(nrow(labels)), function(i) {
  tid <- labels$id[i]; b <- base[[labels$epitope[i]]]
  a <- b + rnorm(26, 0, 0.05)
  tibble(tcr_id = tid, chain = rep(c("A","B"), each = 13),
         resi = rep(105:117, 2),
         phi_rad = a[1:26], psi_rad = a[c(14:26, 1:13)])
}))
write_tsv(dih, gzfile(file.path(root, "sptce_out", "dihedrals.tsv.gz")))

# --- score it with sptce_metrics.R, as the pipeline does --------------------
source(file.path(repo_root, "sptce_metrics.R"))
lab <- transmute(labels, receptor_id = id, cdr3b = cdr3_aa_beta, epitope,
                 v_gene = v_b_gene)
files <- list.files(file.path(root, "sptce_out", "distance_tables"), full.names = TRUE)
dmap <- load_dihedral_map(file.path(root, "sptce_out", "dihedrals.tsv.gz"))
thresholds <- c(0.25, 0.50, 0.75, 1.00, 1.25)
gates <- c(0.2, 0.3, 0.5)

asg <- bind_rows(
  lapply(thresholds, function(thr) cluster_at(files, lab, thr)),
  lapply(thresholds, function(thr) cluster_at(files, lab, thr, 0.5, dmap)))
table1 <- asg %>% group_by(method, threshold, dih) %>%
  group_modify(~ jm_metrics(.x, denom = N)) %>% ungroup()
sweep <- bind_rows(lapply(gates, function(g)
  jm_metrics(cluster_at(files, lab, 1.0, g, dmap), denom = N) %>% mutate(dih = g)))
sens <- asg %>% group_by(method, threshold) %>%
  group_modify(function(d, k) bind_rows(lapply(unique(labels$epitope), function(e) {
    n_in <- sum(labels$epitope == e)
    tibble(epitope = e, sensitivity = jm_sensitivity(d, e, n_in),
           clustered_frac = n_distinct(d$receptor_id[d$epitope == e]) / n_in)
  }))) %>% ungroup()

cat("\n######## assertions ########\n")
fail <- 0
chk <- function(lbl, cond) {
  if (isTRUE(cond)) cat("  PASS ", lbl, "\n") else { cat("  FAIL ", lbl, "\n"); fail <<- fail + 1 }
}
chk("dihedral map covers every receptor",    length(dmap) == N)
chk("table1 covers all 5 thresholds",        length(unique(table1$threshold)) == 5)
chk("table1 has a gated mode",               any(!is.na(table1$dih)))
chk("check_metrics accepts table1",          isTRUE(check_metrics(table1)))
chk("retention never above 1",               all(table1$retention <= 1))
chk("purity never above 1",                  all(table1$purity <= 1, na.rm = TRUE))
chk("clustered never exceeds denom",         all(table1$n_units_clustered <= N))
chk("retention monotone in threshold",       {
  r <- table1 %>% filter(is.na(dih)) %>% arrange(threshold) %>% pull(retention)
  all(diff(r) >= -1e-9) })
chk("purity high at tight threshold",        {
  p <- table1 %>% filter(is.na(dih), threshold == 0.50) %>% pull(purity)
  length(p) == 1 && p >= 0.99 })
chk("clusters are one bin each",             all(grepl("^bin_[12]\\.tsv\\.gz\\|", asg$cluster_id)))
chk("gate never clusters MORE than RMSD-only", {
  cmp <- table1 %>% select(threshold, dih, n = n_units_clustered) %>%
    mutate(dih = ifelse(is.na(dih), "none", "gated")) %>%
    tidyr::pivot_wider(names_from = dih, values_from = n)
  all(cmp$gated <= cmp$none, na.rm = TRUE) })
chk("tighter gate clusters no more than looser", all(diff(sweep$n_units_clustered) >= 0))
chk("sensitivity never exceeds clustered fraction",
    isTRUE(check_sensitivity(sens)))
chk("sensitivity is positive somewhere",     any(sens$sensitivity > 0))
g <- jm_gene_homogeneity(filter(asg, method == "SPACE-TCR RMSD", threshold == 0.75), "v_gene")
chk("gene homogeneity scores clusters",      nrow(g) == 1 && g$n_clusters_scored > 0)

# --- a published method's own output, scored as the reference does ----------
pub <- tibble(CDR3 = labels$cdr3_aa_beta[1:24], clusters = rep(1:6, each = 4))
ann <- distinct(transmute(labels, CDR3 = cdr3_aa_beta, Epitope = epitope))
r <- jm_reference_metrics(pub, ann, "CDR3", retention = 24 / N)
chk("published: four numbers",               length(r) == 4)
chk("published: purity in (0, 1]",           r[2] > 0 && r[2] <= 1)

cat(sprintf("\n==== e2e: %d failures ====\n", fail))
if (fail > 0) quit(status = 1)
