# Sensitivity figures (Eq. 4: input-unit denominator, per-method units).
# Inputs:  the .rds written by scripts/corrected_sensitivity.R
# Outputs: four PNGs in results/figures/corrected/

suppressMessages({library(dplyr); library(readr); library(tidyr)
                  library(ggplot2); library(stringr)})

# ---- paths -------------------------------------------------------------------
#   Rscript scripts/corrected_sensitivity_figures.R <sens_corrected.rds> <out_dir>
# The .rds is written by scripts/corrected_sensitivity.R; run that first.
args <- commandArgs(trailingOnly = TRUE)
RDS <- if (length(args) >= 1) args[1] else "results/tables/corrected/sens_corrected.rds"
OUT <- if (length(args) >= 2) args[2] else "results/figures/corrected"
if (!file.exists(RDS))
  stop("not found: ", RDS, "\n  run scripts/corrected_sensitivity.R first",
       call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

X   <- readRDS(RDS)
all <- X$all; rk <- X$rank; foc <- X$focal
DPI <- 300

SP <- c("SPACE-TCR RMSD (0.75A)", "SPACE-TCR RMSD+dih (1A)")
ord <- rk$method[order(-rk$mean_corrected)]
lab_unit <- setNames(rk$unit, rk$method)

sp_bold <- function(m) ifelse(m %in% SP, paste0("**", m, "**"), m)
pal <- c("chain" = "#9E7BB5", "receptor (beta key)" = "#5B8FA8",
         "receptor (paired)" = "#C8553D")

fdat <- all %>% filter(epitope %in% foc) %>%
  mutate(method = factor(method, levels = rev(ord)),
         unit   = lab_unit[as.character(method)])

# ---- Fig 3 replacement: heatmap, 11 methods x focal epitopes ---------------
p <- ggplot(fdat, aes(epitope, method, fill = sensitivity_corrected)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", sensitivity_corrected),
                colour = sensitivity_corrected > 0.22), size = 2.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
  scale_fill_gradient(low = "#F7F7F7", high = "#1F4E5F", limits = c(0, NA),
                      name = "sensitivity") +
  labs(title = "Sensitivity by method and focal epitope",
       subtitle = "distinct qualifying units carrying the epitope, over input units carrying it",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(OUT, "fig03_sensitivity_heatmap.png"), p, width = 7.6, height = 4.4, dpi = DPI)

# ---- Fig 18: focal, ranked -------------------------------
rk2 <- rk %>% mutate(method = factor(method, levels = rev(ord)))
p <- ggplot(rk2, aes(y = method)) +
  geom_segment(aes(x = 0, xend = mean_corrected, yend = method),
               colour = "grey80", linewidth = 0.7) +
  geom_point(aes(x = mean_corrected, colour = unit), size = 3.6) +
  geom_text(aes(x = mean_corrected, label = sprintf("%.3f", mean_corrected)),
            vjust = -1.2, size = 2.8, colour = "grey20") +
  scale_colour_manual(values = pal, name = "clustering unit") +
  scale_x_continuous(expand = expansion(mult = c(0.07, 0.07))) +
  labs(title = "Mean sensitivity over the focal epitopes",
       subtitle = "mean of per-epitope sensitivity (Eq. 4) over the six focal epitopes",
       x = "mean sensitivity", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        legend.position = "bottom")
ggsave(file.path(OUT, "fig18_focal_sensitivity.png"), p, width = 7.2, height = 4.6, dpi = DPI)

# ---- all epitopes with >= 20 input units ----------------------------------
big <- all %>% filter(n_input >= 20)
p <- ggplot(big, aes(reorder(method, sensitivity_corrected, FUN = median),
                     sensitivity_corrected, fill = unit)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3, width = 0.65) +
  scale_fill_manual(values = pal, name = "clustering unit") +
  coord_flip() +
  labs(title = "Sensitivity across all epitopes with 20 or more input units",
       subtitle = sprintf("%d epitopes", n_distinct(big$epitope)),
       x = NULL, y = "sensitivity") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        legend.position = "bottom")
ggsave(file.path(OUT, "fig_sensitivity_all_epitopes.png"), p, width = 6.8, height = 4.6, dpi = DPI)

# ---- clustering unit per method ---------------------------------------------
u <- all %>% group_by(method, unit) %>%
  summarise(input_pool = max(n_input), .groups = "drop") %>%
  left_join(all %>% group_by(method) %>%
              summarise(clustered = sum(n_clustered), .groups = "drop"), by = "method")
p <- ggplot(u, aes(reorder(method, input_pool), input_pool, fill = unit)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = input_pool), hjust = -0.22, size = 2.9) +
  scale_fill_manual(values = pal, name = "clustering unit") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Clustering unit and input pool per method",
       subtitle = "input pool for the largest epitope (GILGFVFTL): 938 chains, 455 beta keys, 705 paired receptors.\nRetention and sensitivity are divided by these.",
       x = NULL, y = "input units for GILGFVFTL") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        legend.position = "bottom")
ggsave(file.path(OUT, "fig_unit_mismatch.png"), p, width = 7.0, height = 4.2, dpi = DPI)

cat("wrote:\n"); print(basename(list.files(OUT, "\\.png$")))
