# Small plotting and categorical-test helpers shared by multiple analyses.
#
# Analysis-specific gene selection, statistical questions, plot construction,
# and output filenames remain in scripts/*.R. This file only centralizes visual
# conventions and low-level operations that would otherwise be copied verbatim:
# annotation colors, ComplexHeatmap annotations, and PNG device handling.

make_custom_colors <- function(clin_df, cols) {
  annot_colors <- list()
  palette_map <- list(
    Classification = "Set3", Differentiation = "Set2", NAC = "Pastel1",
    T = "Paired", N = "Dark2", N_status = "Accent", Age_group = "Set1",
    Sex = "Pastel2", Stage = "Spectral", Stage_Group = "BrBG"
  )
  defaults <- c("Set3", "Paired", "Dark2", "Accent")
  default_idx <- 1L
  for (column in cols) {
    if (!column %in% names(clin_df)) next
    if (column %in% c("BMI", "Size", "CA19_9", "CEA") && is.numeric(clin_df[[column]])) {
      values <- clin_df[[column]]
      if (all(is.na(values))) next
      limits <- if (column %in% c("CA19_9", "CEA")) {
        c(min(values, na.rm = TRUE), stats::quantile(values, 0.8, na.rm = TRUE))
      } else c(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
      if (length(unique(limits)) < 2L) {
        delta <- if (limits[[1]] == 0) 1 else abs(limits[[1]]) * 0.01
        limits <- limits[[1]] + c(-delta, delta)
      }
      colors <- switch(column,
        BMI = c("#F0F9E8", "#0868AC"), Size = c("#FEE5D9", "#A50F15"),
        CA19_9 = c("#FFFFE5", "#662506"), CEA = c("#F7FCF5", "#00441B"),
        c("white", "black")
      )
      annot_colors[[column]] <- circlize::colorRamp2(limits, colors)
    } else {
      levels <- sort(unique(stats::na.omit(as.character(clin_df[[column]]))))
      if (!length(levels)) next
      palette <- if (column %in% names(palette_map)) palette_map[[column]] else defaults[[default_idx]]
      available <- suppressWarnings(RColorBrewer::brewer.pal(max(length(levels), 3L), palette))
      assigned <- if (length(levels) <= length(available)) available[seq_along(levels)] else grDevices::colorRampPalette(available)(length(levels))
      names(assigned) <- levels
      annot_colors[[column]] <- assigned
      if (!column %in% names(palette_map)) default_idx <- if (default_idx == length(defaults)) 1L else default_idx + 1L
    }
  }
  annot_colors
}

make_annotation <- function(clinical, columns) {
  valid <- intersect(columns, names(clinical))
  if (!length(valid)) return(NULL)
  ComplexHeatmap::HeatmapAnnotation(
    df = clinical[, valid, drop = FALSE], col = make_custom_colors(clinical, valid),
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(0.3, "cm"), na_col = "grey95"
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
