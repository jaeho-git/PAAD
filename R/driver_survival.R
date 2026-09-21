driver_genes <- c("KRAS", "TP53", "SMAD4", "CDKN2A")

extract_kras_subtype <- function(maf_data) {
  aa_columns <- intersect(c("HGVSp_Short", "Protein_Change", "Amino_Acid_Change", "AAChange", "HGVSp"), names(maf_data))
  kras <- maf_data |> dplyr::filter(Hugo_Symbol == "KRAS")
  if (!nrow(kras)) return(tibble::tibble(Tumor_Sample_Barcode = character(), KRAS_subtype = character(), KRAS_subtype_raw = character()))
  if (length(aa_columns)) {
    kras$aa_string <- apply(as.data.frame(kras)[, aa_columns, drop = FALSE], 1, paste, collapse = " ")
  } else kras$aa_string <- NA_character_
  kras |>
    dplyr::mutate(KRAS_subtype = dplyr::case_when(
      stringr::str_detect(aa_string, "G12D|p\\.Gly12Asp") ~ "G12D",
      stringr::str_detect(aa_string, "G12V|p\\.Gly12Val") ~ "G12V",
      stringr::str_detect(aa_string, "G12R|p\\.Gly12Arg") ~ "G12R",
      stringr::str_detect(aa_string, "G12C|p\\.Gly12Cys") ~ "G12C",
      stringr::str_detect(aa_string, "G12S|p\\.Gly12Ser") ~ "G12S",
      stringr::str_detect(aa_string, "G12A|p\\.Gly12Ala") ~ "G12A",
      stringr::str_detect(aa_string, "G13D|p\\.Gly13Asp") ~ "G13D",
      stringr::str_detect(aa_string, "Q61H|p\\.Gln61His") ~ "Q61H",
      stringr::str_detect(aa_string, "Q61R|p\\.Gln61Arg") ~ "Q61R",
      stringr::str_detect(aa_string, "Q61L|p\\.Gln61Leu") ~ "Q61L",
      as.character(Start_Position) == "25398284" & Reference_Allele == "C" & Tumor_Seq_Allele2 == "T" ~ "G12D",
      as.character(Start_Position) == "25398284" & Reference_Allele == "C" & Tumor_Seq_Allele2 == "A" ~ "G12V",
      as.character(Start_Position) == "25398285" & Reference_Allele == "C" & Tumor_Seq_Allele2 == "G" ~ "G12R",
      as.character(Start_Position) == "25398285" & Reference_Allele == "C" & Tumor_Seq_Allele2 == "A" ~ "G12C",
      as.character(Start_Position) == "25398285" & Reference_Allele == "C" & Tumor_Seq_Allele2 == "T" ~ "G12S",
      as.character(Start_Position) == "25398284" & Reference_Allele == "C" & Tumor_Seq_Allele2 == "G" ~ "G12A",
      as.character(Start_Position) == "25398281" & Reference_Allele == "C" & Tumor_Seq_Allele2 == "T" ~ "G13D",
      as.character(Start_Position) == "25380275" & Reference_Allele == "T" & Tumor_Seq_Allele2 %in% c("G", "A") ~ "Q61H/Q61L/Q61R",
      as.character(Start_Position) %in% c("25380275", "25380276", "25380277") ~ "Q61 other",
      TRUE ~ "Other KRAS mutation"
    )) |>
    dplyr::group_by(Tumor_Sample_Barcode) |>
    dplyr::summarise(KRAS_subtype_raw = paste(sort(unique(KRAS_subtype)), collapse = ";"), .groups = "drop") |>
    dplyr::mutate(KRAS_subtype = dplyr::case_when(KRAS_subtype_raw == "G12D" ~ "G12D", KRAS_subtype_raw == "G12V" ~ "G12V", KRAS_subtype_raw == "G12R" ~ "G12R", TRUE ~ "Other KRAS")) |>
    dplyr::select(Tumor_Sample_Barcode, KRAS_subtype, KRAS_subtype_raw)
}

make_driver_summary <- function(maf_data, clinical) {
  samples <- dplyr::distinct(clinical, Tumor_Sample_Barcode)$Tumor_Sample_Barcode
  statuses <- tidyr::expand_grid(Tumor_Sample_Barcode = samples, Gene = driver_genes) |>
    dplyr::left_join(maf_data |> dplyr::filter(Hugo_Symbol %in% driver_genes) |> dplyr::distinct(Tumor_Sample_Barcode, Gene = Hugo_Symbol) |> dplyr::mutate(Mutated = 1L), by = c("Tumor_Sample_Barcode", "Gene")) |>
    dplyr::mutate(Mutated = tidyr::replace_na(Mutated, 0L)) |>
    tidyr::pivot_wider(names_from = Gene, values_from = Mutated, values_fill = 0L)
  clinical |>
    dplyr::left_join(statuses, by = "Tumor_Sample_Barcode") |>
    dplyr::mutate(dplyr::across(dplyr::all_of(driver_genes), ~tidyr::replace_na(.x, 0L))) |>
    dplyr::rowwise() |>
    dplyr::mutate(driver_mutation_count = sum(dplyr::c_across(dplyr::all_of(driver_genes))), driver_mutation_status = ifelse(driver_mutation_count > 0, "Driver-mutated", "Driver wild type")) |>
    dplyr::ungroup() |>
    dplyr::left_join(extract_kras_subtype(maf_data), by = "Tumor_Sample_Barcode") |>
    dplyr::mutate(
      KRAS_subtype = ifelse(is.na(KRAS_subtype) | KRAS == 0, "WT", KRAS_subtype),
      KRAS_subtype = factor(KRAS_subtype, levels = c("G12D", "G12V", "G12R", "Other KRAS", "WT")),
      driver_count_group = dplyr::case_when(driver_mutation_count <= 1 ~ "0-1", driver_mutation_count == 2 ~ "2", driver_mutation_count >= 3 ~ "3-4"),
      driver_mutation_count_exact = factor(driver_mutation_count, levels = 0:4, labels = paste0(0:4, " Mutations")),
      driver_mutation_count_collapsed = factor(dplyr::case_when(driver_mutation_count == 0 ~ "0 Mutations", driver_mutation_count %in% 1:2 ~ "1-2 Mutations", driver_mutation_count %in% 3:4 ~ "3-4 Mutations"), levels = c("0 Mutations", "1-2 Mutations", "3-4 Mutations"))
    )
}

make_pie_chart <- function(data, group_var, title, filename, level_order = NULL, dpi = 1000L) {
  plot_data <- data |> dplyr::filter(!is.na(.data[[group_var]])) |> dplyr::count(.data[[group_var]], name = "n") |> dplyr::mutate(category = as.character(.data[[group_var]]))
  if (!is.null(level_order)) {
    plot_data <- plot_data |> tidyr::complete(category = level_order, fill = list(n = 0L)) |> dplyr::mutate(category = factor(category, levels = level_order)) |> dplyr::arrange(category)
  } else plot_data <- plot_data |> dplyr::arrange(dplyr::desc(n), category) |> dplyr::mutate(category = factor(category, levels = unique(category)))
  total <- sum(plot_data$n)
  if (!total) return(invisible(NULL))
  spread_positions <- function(y, gap = 0.22, lower = -1.10, upper = 1.10) {
    if (length(y) <= 1L) return(y)
    order <- order(y)
    values <- y[order]
    for (i in seq(2L, length(values))) if (values[i] - values[i - 1L] < gap) values[i] <- values[i - 1L] + gap
    if (max(values) > upper) values <- values - (max(values) - upper)
    if (min(values) < lower) values <- values + (lower - min(values))
    result <- numeric(length(y)); result[order] <- values; result
  }
  pie <- plot_data |> dplyr::filter(n > 0) |> dplyr::mutate(
    pct = n / total * 100, ymax = cumsum(n), ymin = dplyr::lag(ymax, default = 0),
    start = pi / 2 - 2 * pi * ymin / total, end = pi / 2 - 2 * pi * ymax / total,
    mid = (start + end) / 2,
    label = paste0(as.character(category), "\n(", n, ", ", sprintf("%.1f", pct), "%)"),
    label_inside = pct >= 12,
    label_x = ifelse(label_inside, 0.58 * cos(mid), ifelse(cos(mid) >= 0, 1.35, -1.35)),
    label_y = ifelse(label_inside, 0.58 * sin(mid), 1.12 * sin(mid)),
    anchor_x = 0.98 * cos(mid), anchor_y = 0.98 * sin(mid),
    hjust_label = ifelse(label_inside, 0.5, ifelse(label_x > 0, 0, 1))
  )
  if (any(!pie$label_inside)) {
    outside <- which(!pie$label_inside)
    adjusted <- pie[outside, ] |> dplyr::mutate(side = ifelse(label_x > 0, "right", "left")) |>
      dplyr::group_by(side) |> dplyr::mutate(label_y = spread_positions(label_y, gap = 0.24)) |> dplyr::ungroup()
    pie$label_y[outside] <- adjusted$label_y
  }
  wedges <- purrr::map_dfr(seq_len(nrow(pie)), function(index) {
    theta <- seq(pie$start[index], pie$end[index], length.out = 160)
    tibble::tibble(x = c(0, cos(theta), 0), y = c(0, sin(theta), 0), category = pie$category[index])
  })
  inside <- dplyr::filter(pie, label_inside)
  outside <- dplyr::filter(pie, !label_inside)
  plot <- ggplot2::ggplot() + ggplot2::geom_polygon(data = wedges, ggplot2::aes(x, y, group = category, fill = category), color = "white", linewidth = 0.6)
  if (nrow(outside)) plot <- plot + ggplot2::geom_segment(data = outside, ggplot2::aes(x = anchor_x, y = anchor_y, xend = label_x * 0.96, yend = label_y), color = "grey45", linewidth = 0.35)
  if (nrow(inside)) plot <- plot + ggplot2::geom_text(data = inside, ggplot2::aes(label_x, label_y, label = label), hjust = 0.5, vjust = 0.5, lineheight = 0.88, size = 4, color = "black")
  if (nrow(outside)) plot <- plot + ggplot2::geom_label(data = outside, ggplot2::aes(label_x, label_y, label = label, hjust = hjust_label), vjust = 0.5, lineheight = 0.88, size = 3.6, linewidth = 0.2, label.padding = grid::unit(0.15, "lines"), color = "black", fill = "white")
  plot <- plot + ggplot2::coord_equal(xlim = c(-1.70, 1.95), ylim = c(-1.35, 1.35), clip = "off") +
    ggplot2::scale_fill_discrete(drop = FALSE, labels = levels(plot_data$category)) +
    ggplot2::labs(title = paste0(title, " (N=", total, ")"), fill = NULL) + ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16, color = "black"), legend.position = "right", legend.text = ggplot2::element_text(size = 12), legend.key.size = grid::unit(0.55, "cm"), plot.background = ggplot2::element_rect(fill = "white", color = NA), panel.background = ggplot2::element_rect(fill = "white", color = NA), legend.background = ggplot2::element_rect(fill = "white", color = NA), plot.margin = ggplot2::margin(10, 50, 10, 10))
  ggplot2::ggsave(filename, plot, width = 10.2, height = 7.2, dpi = dpi, bg = "white")
  invisible(pie)
}

make_driver_clinicopath_table <- function(driver_data, filename) {
  variables <- list(
    "T stage" = list(column = "T", levels = c("T0", "T1", "T2", "T3", "T4")),
    "N stage" = list(column = "N", levels = c("N0", "N1", "N2")),
    "Tumor grade" = list(column = "Differentiation", levels = c("WD", "MD", "PD")),
    LVI = list(column = "LVI", levels = c("Positive", "Negative")),
    PNI = list(column = "PNI", levels = c("Positive", "Negative")),
    RM = list(column = "RM", levels = c("Positive", "Negative"))
  )
  variables <- variables[vapply(variables, function(x) x$column %in% names(driver_data), logical(1))]
  rows <- purrr::imap_dfr(variables, ~dplyr::bind_rows(
    tibble::tibble(Variable = .y, Category = "Overall", Row_type = "Overall"),
    tibble::tibble(Variable = .y, Category = .x$levels, Row_type = "Category")
  ))
  result <- rows
  for (gene in driver_genes) {
    gene_rows <- rows
    mutant_col <- paste0(gene, " Mutant, n (%)")
    wild_col <- paste0(gene, " Wild type, n (%)")
    p_col <- paste0(gene, " p-value")
    definition_col <- paste0(gene, " p-value definition")
    gene_rows[[mutant_col]] <- NA_character_
    gene_rows[[wild_col]] <- NA_character_
    gene_rows[[p_col]] <- NA_character_
    gene_rows[[definition_col]] <- NA_character_
    for (label in names(variables)) {
      spec <- variables[[label]]
      subset <- driver_data |> dplyr::filter(!is.na(.data[[spec$column]])) |>
        dplyr::mutate(gene_status = factor(ifelse(.data[[gene]] == 1, "Mutant", "Wild type"), levels = c("Mutant", "Wild type")), clinic_category = factor(as.character(.data[[spec$column]]), levels = spec$levels)) |>
        dplyr::filter(!is.na(clinic_category))
      denominator_mutant <- sum(subset$gene_status == "Mutant")
      denominator_wild <- sum(subset$gene_status == "Wild type")
      overall <- gene_rows$Variable == label & gene_rows$Row_type == "Overall"
      gene_rows[[p_col]][overall] <- format_pvalue(safe_cat_test(table(subset$clinic_category, subset$gene_status)))
      gene_rows[[definition_col]][overall] <- paste0(label, " overall distribution: ", gene, " mutant vs wild type")
      gene_rows[[mutant_col]][overall] <- gene_rows[[wild_col]][overall] <- ""
      for (category in spec$levels) {
        index <- gene_rows$Variable == label & gene_rows$Category == category & gene_rows$Row_type == "Category"
        n_mutant <- sum(subset$clinic_category == category & subset$gene_status == "Mutant")
        n_wild <- sum(subset$clinic_category == category & subset$gene_status == "Wild type")
        gene_rows[[mutant_col]][index] <- if (denominator_mutant) sprintf("%d (%.1f)", n_mutant, n_mutant / denominator_mutant * 100) else "0 (0.0)"
        gene_rows[[wild_col]][index] <- if (denominator_wild) sprintf("%d (%.1f)", n_wild, n_wild / denominator_wild * 100) else "0 (0.0)"
        binary <- factor(ifelse(subset$clinic_category == category, category, paste0("non-", category)), levels = c(category, paste0("non-", category)))
        gene_rows[[p_col]][index] <- format_pvalue(safe_cat_test(table(binary, subset$gene_status)))
        gene_rows[[definition_col]][index] <- paste0(category, " vs non-", category, ": ", gene, " mutant vs wild type")
      }
    }
    result <- dplyr::bind_cols(result, dplyr::select(gene_rows, -Variable, -Category, -Row_type))
  }
  writexl::write_xlsx(result, filename)
  display <- tibble::tibble(Characteristic = ifelse(result$Row_type == "Overall", result$Variable, paste0("  ", result$Category)))
  for (gene in driver_genes) for (suffix in c(" Mutant, n (%)", " Wild type, n (%)", " p-value")) display[[paste0(gene, suffix)]] <- result[[paste0(gene, suffix)]]
  header1 <- header2 <- as.list(rep("", ncol(display)))
  names(header1) <- names(header2) <- paste0("V", seq_len(ncol(display)))
  for (i in seq_along(driver_genes)) {
    start <- 2L + (i - 1L) * 3L
    header1[[start]] <- driver_genes[[i]]
    header2[start:(start + 2L)] <- list("Mutant, n (%)", "Wild type, n (%)", "p-value")
  }
  names(display) <- names(header1)
  display_table <- dplyr::bind_rows(tibble::as_tibble(header1), tibble::as_tibble(header2), display)
  writexl::write_xlsx(list("Driver mutation table" = display_table), stringr::str_replace(filename, "\\.xlsx$", "_display.xlsx"), col_names = FALSE)
  result
}

run_survival_plot <- function(data, group_var, time_var, event_var, title, filename, legend_title = group_var, dpi = 1000L) {
  required <- c(group_var, time_var, event_var)
  if (!all(required %in% names(data))) return(invisible(NULL))
  survival_data <- data |> dplyr::select(Tumor_Sample_Barcode, dplyr::all_of(required)) |>
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(.data[[time_var]]), !is.na(.data[[event_var]])) |>
    dplyr::transmute(Tumor_Sample_Barcode, surv_time = as.numeric(.data[[time_var]]), surv_event = as.numeric(.data[[event_var]]), group = as.factor(.data[[group_var]])) |>
    dplyr::filter(surv_time > 0, surv_event %in% c(0, 1)) |> dplyr::group_by(group) |> dplyr::filter(dplyr::n() >= 2L) |> dplyr::ungroup() |>
    dplyr::mutate(group = droplevels(group))
  if (!nrow(survival_data) || dplyr::n_distinct(survival_data$group) < 2L) return(invisible(NULL))
  fit <- survival::survfit(survival::Surv(surv_time, surv_event) ~ group, data = survival_data)
  lr <- survival::survdiff(survival::Surv(surv_time, surv_event) ~ group, data = survival_data)
  p_value <- 1 - stats::pchisq(lr$chisq, length(lr$n) - 1L)
  summary <- summary(fit)
  levels <- levels(survival_data$group)
  curve <- tibble::tibble(time = summary$time, survival = summary$surv, strata = stringr::str_replace(as.character(summary$strata), "^group=", "")) |> dplyr::mutate(strata = factor(strata, levels = levels))
  breaks <- pretty(c(0, max(survival_data$surv_time)), n = 6)
  breaks <- sort(unique(c(0, breaks[breaks >= 0 & breaks <= max(survival_data$surv_time)])))
  risk <- tidyr::expand_grid(group = levels, time = breaks) |> dplyr::rowwise() |> dplyr::mutate(n_risk = sum(survival_data$group == group & survival_data$surv_time >= time)) |> dplyr::ungroup() |> dplyr::mutate(group = factor(group, levels = rev(levels)))
  x_label <- if (stringr::str_detect(time_var, "_m$")) "Time (months)" else "Time (days)"
  legend_rows <- if (length(levels) > 3L) 2L else 1L
  p_curve <- ggplot2::ggplot(curve, ggplot2::aes(time, survival, color = strata)) + ggplot2::geom_step(linewidth = 0.9) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = paste0(title, " (N=", nrow(survival_data), ")"), subtitle = paste0("Log-rank p = ", format_pvalue(p_value)), x = x_label, y = "Survival probability", color = legend_title) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom", legend.title = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 9), legend.box = "vertical", plot.margin = ggplot2::margin(10, 10, 10, 10))
  p_risk <- ggplot2::ggplot(risk, ggplot2::aes(time, group, label = n_risk)) + ggplot2::geom_text(size = 3.3) + ggplot2::scale_x_continuous(breaks = breaks, limits = range(breaks)) + ggplot2::labs(title = "Number at risk", x = x_label, y = NULL) + ggplot2::theme_classic(base_size = 11) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11), axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  open_png(filename, 10.5, 8.8, dpi)
  on.exit(close_device(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = grid::unit(c(0.70, 0.30), "npc"))))
  print(p_curve, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_risk, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  close_device()
  on.exit(NULL, add = FALSE)
  readr::write_tsv(dplyr::mutate(risk, group = as.character(group)), stringr::str_replace(filename, "\\.png$", "_number_at_risk.tsv"))
  invisible(fit)
}

run_driver_kras_outputs <- function(maf, clinical, output_dir, dpi = 1000L) {
  driver_data <- make_driver_summary(maf@data, clinical)
  writexl::write_xlsx(driver_data, file.path(output_dir, "11_Driver_KRAS_patient_level_summary.xlsx"))
  make_pie_chart(driver_data, "KRAS_subtype", "KRAS subtype distribution", file.path(output_dir, "11_KRAS_subtype_pie.png"), c("G12D", "G12V", "G12R", "Other KRAS", "WT"), dpi)
  make_pie_chart(driver_data, "driver_mutation_count_exact", "Driver gene mutation count", file.path(output_dir, "12_Driver_mutation_status_pie.png"), paste0(0:4, " Mutations"), dpi)
  table <- make_driver_clinicopath_table(driver_data, file.path(output_dir, "13_Driver_gene_mutation_clinicopath_table.xlsx"))
  readr::write_tsv(table, file.path(output_dir, "13_Driver_gene_mutation_clinicopath_table.tsv"))
  os_time <- if ("OS_m" %in% names(driver_data) && any(!is.na(driver_data$OS_m))) "OS_m" else "OS_d"
  rfs_time <- if ("RFS_m" %in% names(driver_data) && any(!is.na(driver_data$RFS_m))) "RFS_m" else "RFS_d"
  writexl::write_xlsx(driver_data, file.path(output_dir, "14_17_PDAC_survival_input_patient_level.xlsx"))
  if (all(c(os_time, "survive") %in% names(driver_data))) {
    run_survival_plot(driver_data, "KRAS_subtype", os_time, "survive", "OS by KRAS subtype", file.path(output_dir, "14_OS_by_KRAS_subtype.png"), "KRAS subtype", dpi)
    run_survival_plot(driver_data, "driver_mutation_count_exact", os_time, "survive", "OS by driver mutation count", file.path(output_dir, "15_OS_by_driver_mutation_count_exact.png"), "Driver mutation count", dpi)
    run_survival_plot(driver_data, "driver_mutation_count_collapsed", os_time, "survive", "OS by driver mutation count", file.path(output_dir, "15_2_OS_by_driver_mutation_count_0_12_34.png"), "Driver mutation count", dpi)
  }
  if (all(c(rfs_time, "Recur") %in% names(driver_data))) {
    run_survival_plot(driver_data, "KRAS_subtype", rfs_time, "Recur", "RFS by KRAS subtype", file.path(output_dir, "16_RFS_by_KRAS_subtype.png"), "KRAS subtype", dpi)
    run_survival_plot(driver_data, "driver_mutation_count_exact", rfs_time, "Recur", "RFS by driver mutation count", file.path(output_dir, "17_RFS_by_driver_mutation_count_exact.png"), "Driver mutation count", dpi)
    run_survival_plot(driver_data, "driver_mutation_count_collapsed", rfs_time, "Recur", "RFS by driver mutation count", file.path(output_dir, "17_2_RFS_by_driver_mutation_count_0_12_34.png"), "Driver mutation count", dpi)
  }
  invisible(driver_data)
}
