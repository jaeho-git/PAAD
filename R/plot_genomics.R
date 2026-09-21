draw_custom_oncoplot <- function(maf, top_n, title, filename, anno_cols, dpi = 1000L) {
  top_genes <- maftools::getGeneSummary(maf) |> dplyr::slice_head(n = top_n) |> dplyr::pull(Hugo_Symbol)
  if (!length(top_genes)) stop("Oncoplot requires at least one mutation.", call. = FALSE)
  matrix_data <- maf@data |>
    dplyr::filter(Hugo_Symbol %in% top_genes) |>
    dplyr::transmute(
      Hugo_Symbol, Tumor_Sample_Barcode,
      Type = dplyr::case_when(
        Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "Frameshift",
        Variant_Classification == "Nonsense_Mutation" ~ "Nonsense",
        Variant_Classification == "Missense_Mutation" ~ "Missense",
        Variant_Classification == "Splice_Site" ~ "Splice_Site",
        Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins") ~ "In_Frame",
        TRUE ~ "Multi_Hit"
      )
    ) |>
    dplyr::distinct()
  matrix <- matrix_data |>
    tidyr::pivot_wider(names_from = Tumor_Sample_Barcode, values_from = Type, values_fn = function(x) paste(x, collapse = ";")) |>
    tibble::column_to_rownames("Hugo_Symbol") |>
    as.matrix()
  clinical <- as.data.frame(maf@clinical.data) |>
    dplyr::filter(Tumor_Sample_Barcode %in% colnames(matrix)) |>
    tibble::column_to_rownames("Tumor_Sample_Barcode")
  clinical <- clinical[colnames(matrix), , drop = FALSE]
  class_text <- if ("Classification" %in% names(clinical)) {
    counts <- table(clinical$Classification)
    paste(names(counts), as.integer(counts), sep = ":", collapse = ", ")
  } else ""
  final_title <- paste0(title, "\n(Total N=", ncol(matrix), if (nzchar(class_text)) paste0(" | ", class_text) else "", ")")
  annotation <- make_annotation(clinical, anno_cols)
  colors <- c(Frameshift = "#FF7F50", Nonsense = "#DC143C", Missense = "#3CB371", Splice_Site = "#9370DB", In_Frame = "#FFD700", Multi_Hit = "#000000")
  alter <- c(
    list(background = function(x, y, w, h) grid::grid.rect(x, y, w * 0.9, h * 0.9, gp = grid::gpar(fill = "#FAFAFA", col = NA))),
    lapply(names(colors), function(type) {
      force(type)
      function(x, y, w, h) grid::grid.rect(x, y, w * 0.9, h * 0.4, gp = grid::gpar(fill = colors[[type]], col = NA))
    })
  )
  names(alter)[-1] <- names(colors)
  open_png(filename, 14, 9, dpi)
  on.exit(close_device(), add = TRUE)
  heatmap <- ComplexHeatmap::oncoPrint(
    matrix, alter_fun = alter, col = colors, top_annotation = annotation,
    column_title = final_title, column_title_gp = grid::gpar(fontsize = 10),
    pct_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 9),
    show_column_names = FALSE
  )
  ComplexHeatmap::draw(heatmap, merge_legend = TRUE)
  invisible(filename)
}

draw_clustered_heatmap <- function(maf, top_n, title, filename, anno_cols, dpi = 1000L) {
  top_genes <- maftools::getGeneSummary(maf) |> dplyr::slice_head(n = top_n) |> dplyr::pull(Hugo_Symbol)
  matrix <- maf@data |>
    dplyr::filter(Hugo_Symbol %in% top_genes) |>
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) |>
    dplyr::distinct() |>
    dplyr::mutate(value = 1L) |>
    tidyr::pivot_wider(names_from = Tumor_Sample_Barcode, values_from = value, values_fill = 0L) |>
    tibble::column_to_rownames("Hugo_Symbol") |>
    as.matrix()
  if (!nrow(matrix) || !ncol(matrix)) stop("Clustered heatmap requires mutation data.", call. = FALSE)
  clinical <- as.data.frame(maf@clinical.data) |>
    dplyr::filter(Tumor_Sample_Barcode %in% colnames(matrix)) |>
    tibble::column_to_rownames("Tumor_Sample_Barcode")
  clinical <- clinical[colnames(matrix), , drop = FALSE]
  annotation <- make_annotation(clinical, anno_cols)
  open_png(filename, 14, 9, dpi)
  on.exit(close_device(), add = TRUE)
  heatmap <- ComplexHeatmap::Heatmap(
    matrix, name = "Mutation", col = c("0" = "white", "1" = "black"),
    top_annotation = annotation, column_title = paste0(title, " (Top ", top_n, ", Total N=", ncol(matrix), ")"),
    cluster_rows = nrow(matrix) > 1L, cluster_columns = ncol(matrix) > 1L,
    show_column_names = FALSE, rect_gp = grid::gpar(col = "grey90", lwd = 0.5)
  )
  ComplexHeatmap::draw(heatmap, merge_legend = TRUE)
  invisible(filename)
}

run_comparison_analysis <- function(maf, clinical, output_dir, prefix, compare_cols, top_n = 20L, dpi = 1000L) {
  top_genes <- maftools::getGeneSummary(maf) |> dplyr::slice_head(n = top_n) |> dplyr::pull(Hugo_Symbol)
  maf_long <- maf@data |> dplyr::filter(Hugo_Symbol %in% top_genes) |> dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) |> dplyr::distinct()
  outputs <- character()
  for (variable in compare_cols) {
    if (!variable %in% names(clinical)) next
    clinical_sub <- clinical |> dplyr::filter(!is.na(.data[[variable]])) |> dplyr::mutate("{variable}" := as.character(.data[[variable]]))
    levels <- get_expected_levels(variable, clinical_sub[[variable]])
    observed <- intersect(levels, unique(clinical_sub[[variable]]))
    if (!length(observed)) next
    clinical_sub <- clinical_sub |> dplyr::mutate("{variable}" := factor(.data[[variable]], levels = levels))
    sizes <- table(factor(clinical_sub[[variable]], levels = levels))
    subtitle <- paste0("Total N=", sum(sizes), " (", paste(names(sizes), as.integer(sizes), sep = "=", collapse = ", "), ")")
    joined <- dplyr::inner_join(maf_long, clinical_sub, by = "Tumor_Sample_Barcode")
    statistics <- tibble::tibble()
    if (length(observed) >= 2L) {
      statistics <- purrr::map_dfr(top_genes, function(gene) {
        mutated <- joined |> dplyr::filter(Hugo_Symbol == gene) |> dplyr::pull(Tumor_Sample_Barcode)
        mut_counts <- vapply(observed, function(group) sum(as.character(clinical_sub[[variable]]) == group & clinical_sub$Tumor_Sample_Barcode %in% mutated), numeric(1))
        wild_counts <- as.numeric(sizes[observed]) - mut_counts
        p <- tryCatch(stats::fisher.test(rbind(mut_counts, wild_counts), workspace = 2e5, simulate.p.value = TRUE)$p.value, error = function(e) NA_real_)
        tibble::tibble(Hugo_Symbol = gene, Significance = dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"))
      })
    }
    frequency <- joined |>
      dplyr::group_by(.data[[variable]], Hugo_Symbol) |>
      dplyr::summarise(Mutated = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") |>
      tidyr::complete(!!rlang::sym(variable) := factor(levels, levels = levels), Hugo_Symbol = top_genes, fill = list(Mutated = 0L)) |>
      dplyr::mutate(Total = as.numeric(sizes[as.character(.data[[variable]])]), Frequency = ifelse(Total > 0, Mutated / Total * 100, 0), Hugo_Symbol = factor(Hugo_Symbol, levels = top_genes), "{variable}" := factor(as.character(.data[[variable]]), levels = levels))
    bar <- ggplot2::ggplot(frequency, ggplot2::aes(x = Hugo_Symbol, y = Frequency, fill = .data[[variable]])) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.9), width = 0.8) +
      ggplot2::scale_fill_brewer(palette = "Pastel1", drop = FALSE) + ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::labs(title = paste0("[", prefix, "] Mutation Frequency by ", variable), subtitle = subtitle, y = "Frequency (%)", x = "Gene") +
      ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"))
    if (nrow(statistics)) {
      maxima <- frequency |> dplyr::group_by(Hugo_Symbol) |> dplyr::summarise(Max = max(Frequency), .groups = "drop")
      bar <- bar + ggplot2::geom_text(data = dplyr::inner_join(statistics, maxima, by = "Hugo_Symbol"), ggplot2::aes(x = Hugo_Symbol, y = Max + 3, label = Significance), inherit.aes = FALSE, size = 3, color = "red")
    }
    bar_file <- file.path(output_dir, paste0("4_FreqPlot_", prefix, "_", variable, ".png"))
    ggplot2::ggsave(bar_file, bar, width = 14, height = 7, dpi = dpi, bg = "white")

    counts <- maf@data |> dplyr::count(Tumor_Sample_Barcode, name = "TMB") |>
      dplyr::right_join(clinical_sub, by = "Tumor_Sample_Barcode") |>
      dplyr::mutate(TMB = tidyr::replace_na(TMB, 0L), "{variable}" := factor(as.character(.data[[variable]]), levels = levels))
    label <- "Statistical test not performed"
    if (length(observed) >= 2L) {
      formula <- stats::as.formula(paste0("TMB ~ `", variable, "`"))
      test <- tryCatch(if (length(observed) == 2L) stats::wilcox.test(formula, data = counts |> dplyr::filter(.data[[variable]] %in% observed)) else stats::kruskal.test(formula, data = counts |> dplyr::filter(.data[[variable]] %in% observed)), error = function(e) NULL)
      if (!is.null(test)) label <- paste0(if (length(observed) == 2L) "Wilcoxon" else "Kruskal", " P: ", sprintf("%.4f", test$p.value))
    }
    box <- ggplot2::ggplot(counts, ggplot2::aes(x = .data[[variable]], y = TMB, fill = .data[[variable]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.65) + ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
      ggplot2::scale_fill_brewer(palette = "Pastel1", drop = FALSE) + ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::labs(title = paste0("[", prefix, "] TMB Distribution by ", variable), subtitle = paste0(subtitle, "\n", label), y = "TMB") + ggplot2::theme_classic()
    box_file <- file.path(output_dir, paste0("4_TMB_BoxPlot_", prefix, "_", variable, ".png"))
    ggplot2::ggsave(box_file, box, width = 6, height = 6, dpi = dpi, bg = "white")
    outputs <- c(outputs, bar_file, box_file)
  }
  invisible(outputs)
}

draw_tmb_heatmap <- function(maf, title, filename, anno_cols, dpi = 1000L) {
  counts <- maf@data |> dplyr::count(Tumor_Sample_Barcode, name = "TMB")
  clinical <- as.data.frame(maf@clinical.data) |> dplyr::left_join(counts, by = "Tumor_Sample_Barcode") |> dplyr::mutate(TMB = tidyr::replace_na(TMB, 0L)) |> dplyr::arrange(dplyr::desc(TMB))
  matrix <- matrix(clinical$TMB, nrow = 1L, dimnames = list("TMB", clinical$Tumor_Sample_Barcode))
  annotation_data <- clinical |> tibble::column_to_rownames("Tumor_Sample_Barcode")
  annotation <- make_annotation(annotation_data, anno_cols)
  limits <- range(matrix, na.rm = TRUE)
  if (limits[[1]] == limits[[2]]) limits <- limits + c(-0.5, 0.5)
  open_png(filename, 14, 7, dpi)
  on.exit(close_device(), add = TRUE)
  heatmap <- ComplexHeatmap::Heatmap(
    matrix, name = "TMB", col = circlize::colorRamp2(limits, c("#F7FBFF", "#08306B")),
    top_annotation = annotation, column_title = paste0(title, " (N=", ncol(matrix), ")"),
    cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 12, fontface = "bold"), rect_gp = grid::gpar(col = "white", lwd = 0.5)
  )
  ComplexHeatmap::draw(heatmap, merge_legend = TRUE)
  invisible(filename)
}

write_clinical_summary <- function(clinical, anno_vars, filename) {
  summaries <- lapply(anno_vars, function(column) {
    if (is.numeric(clinical[[column]])) {
      clinical |> dplyr::summarise(Var = column, Mean = mean(.data[[column]], na.rm = TRUE), Median = stats::median(.data[[column]], na.rm = TRUE), SD = stats::sd(.data[[column]], na.rm = TRUE))
    } else clinical |> dplyr::count(.data[[column]]) |> dplyr::mutate(Pct = n / sum(n) * 100, Var = column)
  })
  names(summaries) <- anno_vars
  writexl::write_xlsx(summaries, filename)
  invisible(filename)
}

draw_somatic_interactions <- function(maf, figure_file, table_file, dpi = 1000L) {
  open_png(figure_file, 12, 10, dpi)
  on.exit(close_device(), add = TRUE)
  result <- maftools::somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.01), fontSize = 0.6)
  close_device()
  on.exit(NULL, add = FALSE)
  writexl::write_xlsx(as.data.frame(result), table_file)
  invisible(c(figure_file, table_file))
}

draw_stage_counts <- function(clinical, filename, dpi = 1000L) {
  levels <- get_expected_levels("Stage", clinical$Stage)
  counts <- clinical |> dplyr::filter(!is.na(Stage)) |> dplyr::mutate(Stage = factor(as.character(Stage), levels = levels)) |> dplyr::count(Stage, .drop = FALSE) |> dplyr::mutate(n = tidyr::replace_na(n, 0L))
  plot <- ggplot2::ggplot(counts, ggplot2::aes(x = Stage, y = n, fill = Stage)) +
    ggplot2::geom_col(width = 0.75) + ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.5) +
    ggplot2::scale_fill_brewer(palette = "Spectral", drop = FALSE) + ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(title = "PDAC Patient Count by Stage", y = "Count") + ggplot2::theme_classic()
  ggplot2::ggsave(filename, plot, width = 8, height = 6, dpi = dpi, bg = "white")
  invisible(filename)
}
