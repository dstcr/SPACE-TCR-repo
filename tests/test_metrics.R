# Unit tests for sptce_metrics.R against hand-computed values.
# The repo root is located from this script's path, so it runs from any clone.
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."),
         mustWork = FALSE)
if (is.na(repo_root) || !file.exists(file.path(repo_root, "sptce_metrics.R"))) repo_root <- "."
source(file.path(repo_root, "sptce_metrics.R"))

pass <- 0; fail <- 0
ok <- function(label, got, want, tol = 1e-9) {
  good <- if (is.numeric(got) && is.numeric(want)) isTRUE(abs(got - want) < tol) else isTRUE(all.equal(got, want))
  if (good) { pass <<- pass + 1; cat(sprintf("  PASS  %-52s %s\n", label, format(got))) }
  else { fail <<- fail + 1; cat(sprintf("  FAIL  %-52s got %s want %s\n", label, format(got), format(want))) }
}
expect_error <- function(label, expr) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  if (!is.null(e)) { pass <<- pass + 1; cat(sprintf("  PASS  %-52s stopped\n", label)) }
  else { fail <<- fail + 1; cat(sprintf("  FAIL  %-52s did NOT stop\n", label)) }
}

cat("\n== cluster_table / purity ==\n")
# c1: 3 members, 2x epiA 1x epiB -> size 3, n_dom 2, purity 2/3, 2 labels
# c2: 4 members, all epiA        -> size 4, n_dom 4, purity 1.0, 1 label
# c3: singleton                  -> dropped at min_size 2
asg <- tibble::tibble(
  cluster_id  = c("c1","c1","c1","c2","c2","c2","c2","c3"),
  receptor_id = c("r1","r2","r3","r4","r5","r6","r7","r8"),
  epitope     = c("A","A","B","A","A","A","A","A"))
ct <- cluster_table(asg)
ok("clusters kept (singleton dropped)", nrow(ct), 2)
ok("c1 size",   ct$size[ct$cluster_id=="c1"], 3)
ok("c1 n_dom",  ct$n_dom[ct$cluster_id=="c1"], 2)
ok("c1 purity", ct$purity[ct$cluster_id=="c1"], 2/3)
ok("c1 labels", ct$n_labels[ct$cluster_id=="c1"], 2)
ok("c2 purity", ct$purity[ct$cluster_id=="c2"], 1)
# global purity = (2+4)/(3+4) = 6/7
ok("global purity = 6/7", jm_purity(ct), 6/7)
# polyspecificity: 1 of 2 clusters multi-label
ok("polyspecificity = 0.5", jm_polyspecificity(ct), 0.5)

cat("\n== high-purity uses >= 0.9, not > 0.9 ==\n")
# A cluster of exactly 0.9 purity: 10 members, 9 dominant.
asg9 <- tibble::tibble(
  cluster_id  = rep("k", 10),
  receptor_id = paste0("x", 1:10),
  epitope     = c(rep("A", 9), "B"))
m9 <- jm_metrics(asg9, denom = 10)
ok("purity exactly 0.9 counts as high-purity", m9$pct_clusters_hi, 1)

cat("\n== retention is distinct units, not summed sizes ==\n")
# r1 appears in two clusters. Summed sizes = 5; distinct receptors = 4.
asg_dup <- tibble::tibble(
  cluster_id  = c("a","a","b","b","b"),
  receptor_id = c("r1","r2","r1","r3","r4"),
  epitope     = c("A","A","A","A","A"))
md <- jm_metrics(asg_dup, denom = 100)
ok("n_units_clustered = 4 distinct", md$n_units_clustered, 4)
ok("n_rows_clustered  = 5 rows",     md$n_rows_clustered, 5)
ok("retention = 4/100 not 5/100",    md$retention, 0.04)

cat("\n== sensitivity: reference rule, distinct numerator ==\n")
# Input has 10 epiA receptors. One qualifying cluster: size 5, 4 are A (0.8 >= 0.67),
# dominant is A, size >= 4. Recovered distinct A = 4 -> 4/10 = 0.4
asg_s <- tibble::tibble(
  cluster_id  = c(rep("q", 5), rep("p", 4)),
  receptor_id = c(paste0("a", 1:5), paste0("b", 1:4)),
  epitope     = c("A","A","A","A","B",  "B","B","B","A"))
ok("sensitivity = 4/10", jm_sensitivity(asg_s, "A", 10), 0.4)
# size filter: a 3-member all-A cluster must NOT qualify (needs size > 3)
asg_small <- tibble::tibble(cluster_id = rep("s", 3),
                            receptor_id = paste0("z", 1:3), epitope = rep("A", 3))
ok("size-3 cluster does not qualify", jm_sensitivity(asg_small, "A", 10), 0)
# 0.67 threshold, not 2/3: a 4-of-6 cluster is 0.6667 and must be EXCLUDED
asg_23 <- tibble::tibble(
  cluster_id  = rep("t", 6), receptor_id = paste0("w", 1:6),
  epitope     = c("A","A","A","A","B","B"))
ok("4-of-6 (0.6667) excluded by >= 0.67", jm_sensitivity(asg_23, "A", 10), 0)
# 5-of-7 = 0.714 qualifies
asg_57 <- tibble::tibble(
  cluster_id  = rep("u", 7), receptor_id = paste0("v", 1:7),
  epitope     = c(rep("A", 5), "B", "B"))
ok("5-of-7 (0.714) qualifies -> 5/10", jm_sensitivity(asg_57, "A", 10), 0.5)
# Duplicate-row inflation: the same receptor in two qualifying clusters must
# count ONCE. Old row-summing code returned 8/10 here.
asg_inf <- tibble::tibble(
  cluster_id  = c(rep("c1", 4), rep("c2", 4)),
  receptor_id = c("a1","a2","a3","a4", "a1","a2","a3","a4"),
  epitope     = rep("A", 8))
ok("receptor in 2 clusters counted once -> 4/10", jm_sensitivity(asg_inf, "A", 10), 0.4)

cat("\n== sensitivity cannot exceed the clustered fraction ==\n")
s <- jm_sensitivity(asg_s, "A", 10)
cf <- length(unique(asg_s$receptor_id[asg_s$epitope == "A"])) / 10
ok("sensitivity <= clustered fraction", s <= cf + 1e-9, TRUE)

cat("\n== self-checks fire ==\n")
expect_error("check_metrics rejects clustered > denom",
             check_metrics(tibble::tibble(n_units_clustered = 617, denom = 454,
                                          retention = 1.36, purity = 0.9,
                                          pct_clusters_hi = 0.9), "pooled"))
expect_error("check_metrics rejects retention > 1",
             check_metrics(tibble::tibble(n_units_clustered = 10, denom = 100,
                                          retention = 1.2, purity = 0.5,
                                          pct_clusters_hi = 0.5)))
expect_error("check_sensitivity rejects 0.83 vs ceiling 0.311",
             check_sensitivity(tibble::tibble(method = "SPACE-TCR RMSD", epitope = "GILGFVFTL",
                                              sensitivity = 0.83, clustered_frac = 0.311)))
ok("check_sensitivity passes a legal row",
   check_sensitivity(tibble::tibble(sensitivity = 0.3, clustered_frac = 0.311)), TRUE)

cat("\n== gene homogeneity ==\n")
asg_g <- tibble::tibble(
  cluster_id  = c(rep("g1", 10), rep("g2", 4)),
  receptor_id = paste0("r", 1:14),
  epitope     = rep("A", 14),
  v_gene      = c(rep("TRBV19", 9), "TRBV20", rep("TRBV7", 2), "TRBV9", "TRBV9"))
g <- jm_gene_homogeneity(asg_g, "v_gene")
# g1 = 9/10 = 0.9 (>= 0.9 -> hi); g2 = 2/4 = 0.5 (not hi) -> 1 of 2 = 0.5
ok("V-gene high-homogeneity share = 0.5", g$pct_clusters_hi, 0.5)
ok("mean distinct V genes = 2", g$mean_distinct_genes, 2)

cat("\n== within-cluster LD and classes ==\n")
ok("identical sequences: median LD 0",
   ld_stats(rep("CASSLAPGATNEKLFF", 3))$median_ld, 0)
ok("one sequence: no LD", is.na(ld_stats("CASSQDRGTGELFF")$median_ld), TRUE)
ok("two sequences: their distance",
   ld_stats(c("CASSQDRGTGELFF", "CASSQDRGTGEL"))$median_ld, 2)
ok("ld_class_of(0) is identical class", as.character(ld_class_of(0)),
   "Sequence-identical / near-identical")
ok("ld_class_of(2) is similar",  as.character(ld_class_of(2)), "Sequence-similar")
ok("ld_class_of(5) is diverse",  as.character(ld_class_of(5)), "Sequence-diverse")
ok("ld_class_of(8) is unrelated",as.character(ld_class_of(8)), "Sequence-unrelated")

cat("\n== graph helpers ==\n")
# 1-2, 2-3 connected; 4-5 separate -> two components
comp <- connected_components(5, c(1,2,4), c(2,3,5))
ok("two components", length(unique(comp)), 2)
ok("1,2,3 same component", length(unique(comp[1:3])), 1)
# circular distance wraps: 0.1 vs (2pi - 0.1) are 0.2 apart, not ~6.1
ok("circ_dihedral wraps", circ_dihedral(c(0.1), c(2*pi - 0.1)), 0.2, tol = 1e-8)
ok("circ_dihedral identical = 0", circ_dihedral(c(1,2), c(1,2)), 0)

cat("\n== a unit with two epitope labels counts towards each (reference) ==\n")
# m1: 4 CDR3b; x1 is annotated A and B. A: x1,x2,x3 = 3; B: x1,x4 = 2.
asg_ml <- tibble::tibble(
  cluster_id = rep("m1", 5),
  cdr3b      = c("x1", "x1", "x2", "x3", "x4"),
  epitope    = c("A",  "B",  "A",  "A",  "B"))
ct_ml <- cluster_table(asg_ml, unit = "cdr3b")
ok("size counts distinct units (4, not 5 rows)", ct_ml$size, 4)
ok("dominant count is 3 (x1 counts for A)",      ct_ml$n_dom, 3)
ok("two labels",                                  ct_ml$n_labels, 2)

cat("\n== published-method scoring follows purity_functions.R ==\n")
# k1: own size 4, members s1..s3 annotated (s1 A,B; s2 A; s3 A) -> nb_maj 3
# k2: own size 2, no annotated member -> dropped, as the reference's merge does
# k3: own size 1 -> below min size
out <- tibble::tibble(CDR3 = c("s1","s2","s3","s4","u1","u2","z1"),
                      clusters = c("k1","k1","k1","k1","k2","k2","k3"))
ann <- tibble::tibble(CDR3 = c("s1","s1","s2","s3","z1"),
                      Epitope = c("A","B","A","A","A"))
sc <- jm_reference_cluster_scores(out, ann, "CDR3")
ok("unannotated cluster dropped, singleton dropped", nrow(sc), 1)
ok("k1 purity = 3/4 over the method's own size", sc$purity, 0.75)
rm4 <- jm_reference_metrics(out, ann, "CDR3", retention = 0.5)
ok("returns retention, purity, % clusters hi, % TCRs hi",
   isTRUE(all.equal(rm4, c(0.5, 0.75, 0, 0))), TRUE)

cat(sprintf("\n==== %d passed, %d failed ====\n", pass, fail))
if (fail > 0) quit(status = 1)
