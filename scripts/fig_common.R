# =============================================================================
# fig_common.R — shared style for the paper figures
#
# Sourced by scripts/fig01_rebuild.R, scripts/fig_rebuild_all.R and
# scripts/fig_rebuild_supp.R: the paper's theme and shared layout settings.
# =============================================================================

# Set a UTF-8 locale so non-ASCII characters (β, Å, ≥) render correctly when
# the session's locale is C.
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  for (lc in c("C.UTF-8", "C.utf8", "en_US.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)))) break
  }
}
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE))
  warning("no UTF-8 locale available; Greek letters and Angstrom signs will render as '.'")

suppressMessages({
  library(systemfonts); library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(ggrepel); library(forcats); library(scales); library(stringr)
})

# ---- type scale -------------------------------------------------------------
# Figures are drawn on a ~12 in canvas and placed at 6.2 in on the page, so they
# are scaled by about 0.52 in print. Type is set at 1.3x the paper's base sizes
# so it prints near body-text size; a smaller canvas would clip strips and value
# labels.
FS <- 1.3

# DejaVu Sans is used explicitly because the default device font has no Greek
# characters.
BASE_FAMILY <- if ("DejaVu Sans" %in% systemfonts::system_fonts()$family)
  "DejaVu Sans" else ""

theme_sptce <- theme_bw(base_size = 15 * FS, base_family = BASE_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    axis.title       = element_text(size = 16 * FS, face = "bold"),
    axis.text        = element_text(size = 14 * FS, colour = "grey15"),
    legend.title     = element_text(size = 14 * FS, face = "bold"),
    legend.text      = element_text(size = 13 * FS),
    legend.key.size  = unit(1.0, "lines"),
    plot.title       = element_text(size = 17 * FS, face = "bold"),
    strip.text       = element_text(size = 14 * FS, face = "bold"),
    plot.subtitle    = element_text(size = 13 * FS, colour = "grey35"),
    plot.caption     = element_text(size = 14 * FS, colour = "grey50"),
    plot.margin      = margin(8, 12, 8, 8)
  )
theme_set(theme_sptce)

# ---- palette ----------------------------------------------------------------
# The emphasis pair. #6e6e6e has at least 3:1 contrast against white; the pair
# passes colour-vision-deficiency and normal-vision checks.
#
# The house reds #e31a1c and #b30000 differ by only dE 10.6 for normal vision
# (below 15), so use them together only where shape, line type or an axis
# label also separates the two series.
SPTCE_RED  <- "#e31a1c"
SPTCE_DARK <- "#b30000"
OTHER_GREY <- "#6e6e6e"

# An ordered light-to-dark ramp for genuinely ordered categories (size bins,
# polyspecificity classes). One hue family, light = least, dark = most.
ORDINAL_4 <- c("#d9d9d9", "#fdbe85", "#e6550d", "#7f2704")
ORDINAL_6 <- c("#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#7f0000")

pretty_method <- function(x) str_replace_all(as.character(x), "_", " ")

# Whether a row is a SPACE-TCR method; drives the emphasis palette.
is_sptce <- function(x) str_detect(as.character(x), "^SPACE-TCR")

fam <- function(x) factor(ifelse(is_sptce(x), "SPACE-TCR", "Published method"),
                          levels = c("SPACE-TCR", "Published method"))

FAM_FILL <- c("SPACE-TCR" = SPTCE_RED, "Published method" = OTHER_GREY)

save_plot <- function(name, plot, w, h, dpi = 300, out = OUT) {
  ggsave(file.path(out, name), plot, width = w, height = h, dpi = dpi,
         limitsize = FALSE)
  cat("  saved:", name, sprintf("(%g x %g in)", w, h), "\n")
}

# Composite helper: stack panels with baked-in tags. The published figures used
# floating text boxes for the "A"/"B" labels, positioned absolutely, which drift
# if a panel is resized.
stack_panels <- function(..., heights = NULL, collect = FALSE) {
  stopifnot(requireNamespace("patchwork", quietly = TRUE))
  library(patchwork)
  p <- Reduce(`/`, list(...))
  if (!is.null(heights) || collect)
    p <- p + plot_layout(heights = heights,
                         guides = if (collect) "collect" else "keep")
  p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 20 * FS, face = "bold"),
          plot.tag.position = c(0, 1))
}

# Check a computed value against the number given in the paper, so a figure
# always matches its caption.
expect_close <- function(label, got, want, tol = 5e-3) {
  ok <- is.finite(got) && abs(got - want) <= tol
  cat(sprintf("    %-46s got %.4f  want %.4f  %s\n", label, got, want,
              if (ok) "OK" else "MISMATCH"))
  if (!ok) stop("figure value disagrees with the paper: ", label,
                call. = FALSE)
  invisible(TRUE)
}
