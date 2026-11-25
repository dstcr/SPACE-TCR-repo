# Sensitivity (Eq. 4), with each method measured in the unit it was run on:
# distinct qualifying units carrying an epitope over INPUT units carrying it.
#
# The nine published methods were NOT run on one input representation. From the
# output files in the Jouannet repository's Data/Output_methods:
#
#   chain-pooled single chain   GIANA, iSMART, TCRMatch, Hamming, Levenshtein
#     GIANA's own header names the input: .../Dataset/alpha_beta_bis/
#     liste_alpha_beta_bis.txt. Alpha and beta chains are separate records in
#     one list, so a cluster may be all-alpha, and the UNIT IS A CHAIN.
#
#   paired                      clusTCR, GLIPH2, DeepTCR, TCRdist3
#     clusTCR's junction_aa is beta+alpha concatenated; GLIPH2 has TcRb+TcRa;
#     DeepTCR has Alpha_Sequences+Beta_Sequences; TCRdist3 has CDR3a+CDR3b.
#     One beta per receptor, so beta keys the receptor. UNIT IS A RECEPTOR.
#
#   paired, structural          SPACE-TCR. UNIT IS A RECEPTOR.
#
# Each method's numerator and denominator are therefore counted in its own unit;
# mixing a chain-unit numerator with a beta-only denominator gives
# n_clustered > n_input.

source("sptce_metrics.R")

# ---- paths -------------------------------------------------------------------
#   Rscript scripts/corrected_sensitivity.R <labelled.tsv> <stage5.rds> <out_dir>
#
# stage5.rds is the object run_stage5() returns:
#     source("sptce_pipeline.R"); saveRDS(run_stage5(), "stage5.rds")
# Saving it is worth the disk: the clustering takes several minutes and every
# correction below is a recomputation over the same membership.
args <- commandArgs(trailingOnly = TRUE)
LAB <- if (length(args) >= 1) args[1] else "data/jouannet_labelled_input.tsv"
S5F <- if (length(args) >= 2) args[2] else "stage5.rds"
OUT <- if (length(args) >= 3) args[3] else "results/tables/corrected"

for (f in c(LAB, S5F)) if (!file.exists(f))
  stop("not found: ", f, "\n  usage: Rscript scripts/corrected_sensitivity.R ",
       "<labelled.tsv> <stage5.rds> <out_dir>", call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

lab <- read_tsv(LAB, show_col_types = FALSE)
s5  <- readRDS(S5F)

CHAIN_POOLED <- c("GIANA", "iSMART", "TCRMatch", "Hamming_Distance", "Levenshtein_Distance")

lab <- lab %>% filter(!is.na(epitope), nzchar(epitope))

# denominator 1: distinct CHAINS carrying the epitope (alpha and beta pooled)
den_chain <- bind_rows(
    lab %>% transmute(epitope, chain = cdr3_aa_alpha),
    lab %>% transmute(epitope, chain = cdr3_aa_beta)) %>%
  filter(!is.na(chain), nzchar(chain)) %>%
  group_by(epitope) %>% summarise(n_input = n_distinct(chain), .groups = "drop")

# denominator 2: distinct RECEPTORS carrying the epitope, keyed as the method keys them
den_beta <- lab %>% group_by(epitope) %>%
  summarise(n_input = n_distinct(cdr3_aa_beta), .groups = "drop")
den_pair <- lab %>% group_by(epitope) %>%
  summarise(n_input = n_distinct(paste(cdr3_aa_alpha, cdr3_aa_beta)), .groups = "drop")

# Eq. 4 per epitope for one method (jm_sensitivity_detail), with the clustered
# count alongside so the invariant sensitivity <= clustered fraction can be
# checked. Clusters below METRIC$min_cluster_size do not count as clustered.
sens_one <- function(memb, den, method, unit_kind) {
  memb <- memb %>% filter(!is.na(epitope), nzchar(epitope), !is.na(unit), nzchar(unit)) %>%
    distinct(cluster_id, unit, epitope)
  ct <- cluster_table(memb, unit = "unit")
  memb <- memb %>% filter(cluster_id %in% ct$cluster_id)
  clus <- memb %>% group_by(epitope) %>%
    summarise(n_clustered = n_distinct(unit), .groups = "drop")
  num <- map_dfr(den$epitope, function(e)
    tibble(epitope = e,
           n_recovered = jm_sensitivity_detail(memb, e, 1L, unit = "unit")$n_recovered))
  den %>% left_join(num, by = "epitope") %>% left_join(clus, by = "epitope") %>%
    mutate(method = method, unit = unit_kind,
           n_recovered = coalesce(as.integer(n_recovered), 0L),
           n_clustered = coalesce(n_clustered, 0L),
           sensitivity_corrected = n_recovered / n_input,
           clustered_frac        = n_clustered / n_input)
}

res <- list()
pb <- s5$published_clusters
for (m in unique(pb$method)) {
  memb <- pb %>% filter(method == m) %>% transmute(unit = cdr3b, cluster_id, epitope)
  if (m %in% CHAIN_POOLED) res[[m]] <- sens_one(memb, den_chain, m, "chain")
  else                     res[[m]] <- sens_one(memb, den_beta,  m, "receptor (beta key)")
}
sp <- s5$sptce_clusters
for (m in unique(sp$method)) {
  memb <- sp %>% filter(method == m) %>% transmute(unit = paste(cdr3a, cdr3b), cluster_id, epitope)
  res[[m]] <- sens_one(memb, den_pair, m, "receptor (paired)")
}
all <- bind_rows(res)

# ---- invariants --------------------------------------------------------------
stopifnot(all(all$sensitivity_corrected <= 1 + 1e-9))
bad <- all %>% filter(n_clustered > n_input)
cat(sprintf("rows with clustered > input: %d of %d\n", nrow(bad), nrow(all)))
if (nrow(bad)) print(as.data.frame(bad %>% mutate(x = n_clustered - n_input) %>%
                     arrange(desc(x)) %>% head(10) %>%
                     select(method, unit, epitope, n_input, n_clustered)))
stopifnot(all(all$sensitivity_corrected <= all$clustered_frac + 1e-9))
cat("invariants hold\n")

write_csv(all, file.path(OUT, "sensitivity_CORRECTED_all_methods.csv"))

foc <- unique(s5$sensitivity$epitope)
cat("\nfocal set:", paste(foc, collapse = ", "), "\n")

cat("\n==== GILGFVFTL ====\n")
print(as.data.frame(all %>% filter(epitope == "GILGFVFTL") %>%
  select(method, unit, n_input, n_clustered, n_recovered,
         sensitivity_corrected) %>%
  arrange(desc(sensitivity_corrected))), digits = 4)

rk <- all %>% filter(epitope %in% foc) %>% group_by(method, unit) %>%
  summarise(mean_corrected = mean(sensitivity_corrected, na.rm = TRUE), .groups = "drop") %>%
  mutate(rank_corrected = rank(-mean_corrected, ties.method = "min")) %>%
  arrange(rank_corrected)
cat("\n==== ranking over the focal set ====\n"); print(as.data.frame(rk), digits = 4)
write_csv(rk, file.path(OUT, "sensitivity_CORRECTED_ranking.csv"))

tab <- all %>% filter(epitope %in% foc) %>%
  select(method, unit, epitope, n_input, n_clustered, n_recovered,
         sensitivity_corrected) %>% arrange(epitope, desc(sensitivity_corrected))
write_csv(tab, file.path(OUT, "sensitivity_CORRECTED_focal.csv"))
saveRDS(list(all = all, rank = rk, focal = foc, tab = tab),
        file.path(OUT, "sens_corrected.rds"))
cat("\nwrote 3 CSVs to", OUT, "\n")
