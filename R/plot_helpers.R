# Small plotting and categorical-test helpers shared by multiple analyses.
#
# Analysis-specific gene selection, statistical questions, plot construction,
# and output filenames remain in scripts/*.R. This file only centralizes visual
# conventions and low-level operations that would otherwise be copied verbatim:
# annotation colors, ComplexHeatmap annotations, and PNG device handling.
source("R/annotation_palettes.R")

make_custom_colors <- function(clin_df, cols) {
  annot_colors <- list()
  for (column in cols) {
    if (!column %in% names(clin_df)) next
    if (column %in% c("Age", "BMI", "Size", "Tumor_size", "CA19_9", "CEA",
        "TMB", "MAF_variant_count") && is.numeric(clin_df[[column]])) {
      values <- clin_df[[column]]
      if (all(is.na(values))) next
      limits <- if (column %in% c("CA19_9", "CEA")) {
        c(min(values, na.rm = TRUE), stats::quantile(values, 0.8, na.rm = TRUE))
      } else c(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
      if (length(unique(limits)) < 2L) {
        delta <- if (limits[[1]] == 0) 1 else abs(limits[[1]]) * 0.01
        limits <- limits[[1]] + c(-delta, delta)
      }
      colors <- annotation_continuous_palette(column)
      annot_colors[[column]] <- circlize::colorRamp2(limits, colors)
    } else {
      levels <- sort(unique(stats::na.omit(as.character(clin_df[[column]]))))
      if (!length(levels)) next
      annot_colors[[column]] <- annotation_discrete_colors(column,levels)
    }
  }
  annot_colors
}

make_annotation <- function(clinical, columns) {
  valid <- intersect(columns, names(clinical))
  valid <- valid[!grepl("purity|cellularity|cellarity|(^|_)loh($|_)", valid, ignore.case = TRUE)]
  if (!length(valid)) return(NULL)
  ComplexHeatmap::HeatmapAnnotation(
    df = clinical[, valid, drop = FALSE], col = make_custom_colors(clinical, valid),
    annotation_label = if (exists("figure_label")) figure_label(valid) else valid,
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(0.3, "cm"), na_col = "#BDBDBD"
  )
}

open_png <- function(filename, width, height, dpi, bg = "white") {
  grDevices::png(filename, width = width, height = height, units = "in", res = dpi, bg = bg)
  invisible(grDevices::dev.cur())
}

close_device <- function() {
  if (grDevices::dev.cur() > 1L) grDevices::dev.off()
  invisible(NULL)
}
