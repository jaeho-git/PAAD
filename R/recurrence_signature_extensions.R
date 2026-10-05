# Prespecified recurrence-pattern and preoperative-platinum extensions.
# Patient identifiers and patient-level values are written only below private/.

run_recurrence_signature_extensions <- function(d, variants, out, savefig, tab, pal, config) {
  fmt_p <- function(x) ifelse(is.na(x), "NE", ifelse(x < .001, "<0.001", sprintf("%.3f", x)))
  dir.create(file.path(out, "private", "signature"), recursive = TRUE, showWarnings = FALSE)

  save_table_png <- function(df, title, filename, widths = NULL, font_size = 7.6,
                             width = 17, height = NULL) {
    shown <- as.data.frame(df, stringsAsFactors = FALSE)
    stripe <- rep(c("#FFFFFF", "#F4F4F4"), length.out = nrow(shown))
    theme <- gridExtra::ttheme_minimal(
      base_size = font_size, base_family = "Pretendard",
      core = list(bg_params = list(fill = stripe, col = NA), fg_params = list(hjust = 0, x = .02)),
      colhead = list(bg_params = list(fill = "#E6E6E6", col = NA),
        fg_params = list(fontface = "bold", hjust = 0, x = .02)))
    grob <- gridExtra::tableGrob(shown, rows = NULL, theme = theme)
    if (!is.null(widths)) grob$widths <- grid::unit(widths, "null")
    if (is.null(height)) height <- max(3.5,
      grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE) + .8)
    canvas <- gridExtra::arrangeGrob(grob, top = grid::textGrob(title,
      gp = grid::gpar(fontfamily = "Pretendard", fontsize = 16, fontface = "bold"),
      x = .01, hjust = 0))
    ragg::agg_png(file.path(out, "figures", paste0(filename, ".png")),
      width = width, height = height, units = "in", res = 300, background = "white")
    grid::grid.draw(canvas); grDevices::dev.off()
  }

  make_km <- function(z, group_var, family_name, title) {
    z <- z[complete.cases(z[, c(group_var, "post_recurrence_months", "survive")]), , drop = FALSE]
    z <- z[is.finite(z$post_recurrence_months) & z$post_recurrence_months > 0 & z$survive %in% 0:1, , drop = FALSE]
    z$group <- droplevels(factor(z[[group_var]]))
    if (nlevels(z$group) < 2) stop("At least two recurrence groups are required for ", family_name)
    fit <- survival::survfit(survival::Surv(post_recurrence_months, survive) ~ group, data = z,
      conf.type = "log-log")
    lr <- survival::survdiff(survival::Surv(post_recurrence_months, survive) ~ group, data = z)
    global_p <- pchisq(lr$chisq, nlevels(z$group) - 1, lower.tail = FALSE)
    pw <- pairwise_logrank(z$post_recurrence_months, z$survive, z$group)
    st <- summary(fit, censored = TRUE)
    curve <- data.frame(time = st$time, prob = st$surv, censor = st$n.censor,
      group = sub("^group=", "", as.character(st$strata))) |>
      dplyr::bind_rows(data.frame(time = 0, prob = 1, censor = 0, group = levels(z$group))) |>
      dplyr::arrange(group, time)
    max_time <- max(z$post_recurrence_months)
    breaks <- seq(0, floor(max_time / 12) * 12, 12)
    if (length(breaks) > 9) breaks <- pretty(c(0, max_time), n = 7)
    risk <- tidyr::expand_grid(group = levels(z$group), time = breaks) |>
      dplyr::rowwise() |>
      dplyr::mutate(n_risk = sum(z$group == group & z$post_recurrence_months >= time)) |>
      dplyr::ungroup()
    cols <- setNames(rep(pal, length.out = nlevels(z$group)), levels(z$group))
    p1 <- ggplot2::ggplot(curve, ggplot2::aes(time, prob, color = group)) +
      ggplot2::geom_step(linewidth = .9) +
      ggplot2::geom_point(data = curve[curve$censor > 0, ], shape = 3, size = .7) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::scale_x_continuous(breaks = breaks) +
      ggplot2::labs(title = title,
        subtitle = sprintf("N=%d, deaths=%d; global log-rank p %s", nrow(z), sum(z$survive), fmt_p(global_p)),
        x = NULL, y = "Post-recurrence overall survival", color = NULL,
        caption = pairwise_caption(pw, 110)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    risk$group <- factor(risk$group, levels = rev(levels(z$group)))
    p2 <- ggplot2::ggplot(risk, ggplot2::aes(time, group, label = n_risk, color = group)) +
      ggplot2::geom_text(size = 3) + ggplot2::scale_color_manual(values = cols) +
      ggplot2::scale_x_continuous(breaks = breaks) + ggplot2::labs(x = "Months", y = NULL) +
      ggplot2::theme_classic(base_family = "Pretendard", base_size = 10) +
      ggplot2::theme(legend.position = "none", axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank())
    global <- data.frame(family = family_name, variable = group_var,
      endpoint = "death after recorded recurrence", n = nrow(z), deaths = sum(z$survive),
      chisq = unname(lr$chisq), df = nlevels(z$group) - 1, p = global_p, method = "Global log-rank")
    tab(global, paste0("Input_", family_name, "_logrank.tsv"))
    tab(pw, paste0("Input_", family_name, "_pairwise_logrank.tsv"))
    list(plot = cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(.77, .23), align = "v", axis = "lr"),
      data = z, global = global, pairwise = pw)
  }

  make_adjusted_cox <- function(z, group_var, family_name) {
    covariates <- c("splines::ns(Age, df = 3)", "Sex", "T_stage", "N_stage", "M_stage",
      "strata(Differentiation)", "LVI", "PNI", "R_status")
    needed <- unique(c(group_var, "post_recurrence_months", "survive", unlist(lapply(covariates,
      function(x) all.vars(stats::reformulate(x))))))
    z <- droplevels(z[complete.cases(z[, needed, drop = FALSE]), , drop = FALSE])
    z$group <- droplevels(factor(z[[group_var]]))
    z$time <- z$post_recurrence_months; z$event <- z$survive
    covariates <- covariates[vapply(covariates, function(a)
      all(vapply(all.vars(stats::reformulate(a)), function(v) length(unique(z[[v]])) > 1, logical(1))), logical(1))]
    full <- survival::coxph(stats::reformulate(c("group", covariates),
      response = "survival::Surv(time,event)"), data = z, x = TRUE, model = TRUE, ties = "efron")
    reduced <- survival::coxph(stats::reformulate(covariates,
      response = "survival::Surv(time,event)"), data = z, x = TRUE, model = TRUE, ties = "efron")
    lrt <- stats::anova(reduced, full, test = "LRT")
    pairs <- pairwise_cox(full, z, "group")
    ph <- tryCatch(survival::cox.zph(full)$table, error = function(e) NULL)
    pairs$family <- family_name
    pairs$n <- nrow(z); pairs$deaths <- sum(z$event)
    pairs$global_LRT_p <- lrt$`Pr(>|Chi|)`[2]
    pairs$PH_global_p <- if (is.null(ph)) NA_real_ else ph["GLOBAL", "p"]
    pairs$comparison_group <- pairs$group1; pairs$reference_group <- pairs$group2
    pairs$ratio_definition <- "Death hazard after recurrence in comparison_group / reference_group"
    pairs$direction <- ratio_direction(pairs$HR, pairs$group1, pairs$group2,
      "death after recorded recurrence", pairs$p_holm < .05)
    pairs$formula <- paste(deparse(formula(full)), collapse = " ")
    tab(pairs, paste0("Input_", family_name, "_Cox_models.tsv"))
    pairs
  }

  recurrence_extension <- function(z, group_var, family_name, title, file_stub) {
    dist <- z |> dplyr::filter(!is.na(.data[[group_var]])) |>
      dplyr::count(.data[[group_var]], name = "n") |>
      dplyr::mutate(percent = 100 * n / sum(n))
    names(dist)[1] <- "group"
    tab(dist, paste0("Input_", family_name, "_distribution.tsv"))
    km <- make_km(z, group_var, family_name, title)
    cx <- make_adjusted_cox(km$data, group_var, family_name)
    forest <- forest_result_plot(cx, paste0(title, ": adjusted Cox comparisons"),
      paste0(cx$comparison_group, "\nReference: ", cx$reference_group),
      p_column = "p_holm", p_heading = "Holm p") +
      ggplot2::labs(caption = paste0(
        "HR = comparison / reference death hazard after recorded recurrence. HR > 1 indicates higher hazard in the comparison group.\n",
        "Adjusted for age (spline), sex, T/N/M category, differentiation strata, LVI, PNI, and resection margin; pointwise 95% CI."))
    display <- cx |> dplyr::transmute(
      Comparison = comparison_group, Reference = reference_group,
      `n / deaths` = paste0(n, " / ", deaths),
      `HR [95% CI]` = sprintf("%.2f [%.2f, %.2f]", HR, lower95, upper95),
      `Raw p` = fmt_p(p), `Holm p` = fmt_p(p_holm), Direction = direction)
    savefig(km$plot, paste0(file_stub, "_KM"), 11, 8)
    savefig(forest, paste0(file_stub, "_forest"), 13, max(5, 2.4 + .55 * nrow(cx)))
    save_table_png(display, paste0(title, ". Adjusted Cox model"), paste0(file_stub, "_Cox_table"),
      widths = c(1.5, 1.5, 1, 1.5, .8, .8, 4), width = 18)
    list(distribution = dist, km = km, cox = cx, forest = forest, display = display)
  }

  recurrent <- d[d$Recur == 1 & !is.na(d$Recurrence_pattern), , drop = FALSE]
  recurrent$post_recurrence_months <- recurrent$OS_m - recurrent$RFS_m
  recurrent <- recurrent[is.na(recurrent$post_recurrence_months) | recurrent$post_recurrence_months > 0, , drop = FALSE]

  distant <- recurrent[recurrent$Recurrence_pattern == "Distant only" & !is.na(recurrent$Distant_pattern), , drop = FALSE]
  distant$Distant_pattern_group <- factor(ifelse(distant$Distant_pattern == "Liver only", "Liver only",
    ifelse(distant$Distant_pattern == "Lung only", "Lung only", "Other distant")),
    levels = c("Liver only", "Lung only", "Other distant"))
  distant_map <- data.frame(source_value = c("Liver only", "Lung only", "Peritoneum only", "Multiple distant", "Other distant"),
    analysis_group = c("Liver only", "Lung only", "Other distant", "Other distant", "Other distant"))
  tab(distant_map, "Input_Distant_pattern_grouping_rule.tsv")
  distant_result <- recurrence_extension(distant, "Distant_pattern_group", "Distant_only_pattern_post_OS",
    "Post-recurrence survival within distant-only recurrence", "Figure4_Distant_only_pattern")

  recurrent$Expanded_recurrence_group <- dplyr::case_when(
    recurrent$Recurrence_pattern == "Local only" ~ "Local only",
    recurrent$Recurrence_pattern == "Local + distant" ~ "Local + distant",
    recurrent$Recurrence_pattern == "Distant only" & recurrent$Distant_pattern == "Liver only" ~ "Liver only",
    recurrent$Recurrence_pattern == "Distant only" & recurrent$Distant_pattern == "Lung only" ~ "Lung only",
    recurrent$Recurrence_pattern == "Distant only" & !is.na(recurrent$Distant_pattern) ~ "Other distant",
    TRUE ~ NA_character_)
  recurrent$Expanded_recurrence_group <- factor(recurrent$Expanded_recurrence_group,
    levels = c("Local only", "Local + distant", "Liver only", "Lung only", "Other distant"))
  expanded <- recurrent[!is.na(recurrent$Expanded_recurrence_group), , drop = FALSE]
  expanded_result <- recurrence_extension(expanded, "Expanded_recurrence_group", "Expanded_recurrence_post_OS",
    "Post-recurrence survival by expanded recurrence phenotype", "Figure4_Expanded_recurrence")
  savefig(cowplot::plot_grid(distant_result$km$plot, distant_result$forest,
    expanded_result$km$plot, expanded_result$forest, ncol = 1,
    rel_heights = c(.27, .20, .29, .24)), "Figure4_Recurrence_extensions", 15, 28)

  # HRD/MMR-list retained-variant comparisons by preoperative platinum exposure.
  platinum <- d[d$Preop_platinum_exposure %in% c("No", "Yes"), , drop = FALSE]
  platinum$Preop_platinum <- factor(platinum$Preop_platinum_exposure, levels = c("No", "Yes"))
  sets <- repair_gene_sets()
  vv <- variants
  vv$paper_gene <- repair_paper_symbol(vv$Hugo_Symbol)
  gene_rows <- list()
  for (set_name in names(sets)) for (gene in sets[[set_name]]) {
    flag <- as.integer(platinum$Tumor_Sample_Barcode %in%
      vv$Tumor_Sample_Barcode[vv$paper_gene == gene])
    tbl <- table(factor(flag, levels = 0:1), platinum$Preop_platinum)
    ft <- tryCatch(stats::fisher.test(tbl), error = function(e) NULL)
    gene_rows[[paste(set_name, gene)]] <- data.frame(gene_set = set_name, gene = gene,
      no_n = sum(platinum$Preop_platinum == "No"), yes_n = sum(platinum$Preop_platinum == "Yes"),
      no_variant_n = sum(flag == 1 & platinum$Preop_platinum == "No"),
      yes_variant_n = sum(flag == 1 & platinum$Preop_platinum == "Yes"),
      no_variant_percent = 100 * mean(flag[platinum$Preop_platinum == "No"]),
      yes_variant_percent = 100 * mean(flag[platinum$Preop_platinum == "Yes"]),
      odds_ratio_yes_vs_no = if (is.null(ft)) NA_real_ else unname(ft$estimate),
      lower95 = if (is.null(ft)) NA_real_ else ft$conf.int[1],
      upper95 = if (is.null(ft)) NA_real_ else ft$conf.int[2],
      p = if (is.null(ft)) NA_real_ else ft$p.value)
  }
  gene_results <- dplyr::bind_rows(gene_rows) |>
    dplyr::group_by(gene_set) |> dplyr::mutate(q_BH_within_gene_set = p.adjust(p, "BH")) |>
    dplyr::ungroup()
  tab(gene_results, "Input_Preop_platinum_HRD_MMR_gene_tests.tsv")

  status_rows <- list()
  for (set_name in names(sets)) {
    status_var <- paste0(set_name, "_variant_status")
    count_var <- paste0(set_name, "_gene_count")
    status <- factor(platinum[[status_var]], levels = c("No retained variant", "Variant detected"))
    ft <- fisher.test(table(status, platinum$Preop_platinum))
    w <- suppressWarnings(wilcox.test(platinum[[count_var]] ~ platinum$Preop_platinum, exact = FALSE))
    tab2 <- table(status, platinum$Preop_platinum)
    status_rows[[paste0(set_name, "_status")]] <- data.frame(gene_set = set_name,
      outcome = "Any retained variant", method = "Fisher exact", p = ft$p.value,
      effect = unname(ft$estimate), effect_definition = "Odds ratio: Yes / No",
      no_n = colSums(tab2)["No"], yes_n = colSums(tab2)["Yes"],
      no_value = 100 * tab2["Variant detected", "No"] / colSums(tab2)["No"],
      yes_value = 100 * tab2["Variant detected", "Yes"] / colSums(tab2)["Yes"])
    status_rows[[paste0(set_name, "_count")]] <- data.frame(gene_set = set_name,
      outcome = "Retained variant gene count", method = "Wilcoxon rank-sum", p = w$p.value,
      effect = unname(w$statistic), effect_definition = "Wilcoxon W",
      no_n = sum(platinum$Preop_platinum == "No"), yes_n = sum(platinum$Preop_platinum == "Yes"),
      no_value = median(platinum[[count_var]][platinum$Preop_platinum == "No"], na.rm = TRUE),
      yes_value = median(platinum[[count_var]][platinum$Preop_platinum == "Yes"], na.rm = TRUE))
  }
  status_results <- dplyr::bind_rows(status_rows)
  status_results$p_holm <- p.adjust(status_results$p, "holm")
  tab(status_results, "Input_Preop_platinum_HRD_MMR_set_tests.tsv")

  status_plot_data <- dplyr::bind_rows(lapply(names(sets), function(set_name) {
    status_var <- paste0(set_name, "_variant_status")
    platinum |> dplyr::count(Preop_platinum, status = .data[[status_var]], name = "n") |>
      dplyr::group_by(Preop_platinum) |> dplyr::mutate(percent = 100 * n / sum(n), gene_set = set_name) |>
      dplyr::ungroup()
  }))
  ps <- ggplot2::ggplot(status_plot_data,
    ggplot2::aes(Preop_platinum, percent, fill = status)) +
    ggplot2::geom_col(width = .65) +
    ggplot2::geom_text(ggplot2::aes(label = n), position = ggplot2::position_stack(vjust = .5), size = 3.4) +
    ggplot2::facet_wrap(~gene_set) +
    ggplot2::scale_fill_manual(values = c("No retained variant" = "#D8D8D8", "Variant detected" = pal[2])) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
    ggplot2::labs(x = "Preoperative platinum exposure", y = "Patients (%)", fill = NULL,
      title = "HRD- and MMR-list retained variants by preoperative platinum exposure",
      subtitle = "Variant-list membership is not a functional HRD or dMMR diagnosis",
      caption = paste(vapply(names(sets), function(s) {
        x <- status_results[status_results$gene_set == s & status_results$outcome == "Any retained variant", ]
        paste0(s, ": Fisher p ", fmt_p(x$p), "; Holm p ", fmt_p(x$p_holm),
          "; OR (Yes/No) ", sprintf("%.2f", x$effect))
      }, character(1)), collapse = " | ")) + ggplot2::theme(legend.position = "top")

  gene_plot <- gene_results |> dplyr::filter(no_variant_n + yes_variant_n > 0) |>
    dplyr::mutate(label = paste0(gene, " (", gene_set, ")"),
      label = factor(label, levels = rev(label[order(odds_ratio_yes_vs_no, na.last = TRUE)])))
  pg <- ggplot2::ggplot(gene_plot, ggplot2::aes(odds_ratio_yes_vs_no, label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "#777777") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower95, xmax = upper95), orientation = "y", width = .18) +
    ggplot2::geom_point(ggplot2::aes(color = gene_set), size = 2.6) +
    ggplot2::scale_x_log10() + ggplot2::scale_color_manual(values = c(HRD = pal[2], MMR = pal[3])) +
    ggplot2::labs(x = "Odds ratio for variant detection (Yes / No, log scale)", y = NULL, color = NULL,
      title = "Gene-level retained-variant comparisons",
      caption = "Pointwise 95% CI; q values use BH correction within each prespecified gene list. Sparse estimates are imprecise.") +
    ggplot2::theme(legend.position = "top")
  savefig(ps, "Figure3_Platinum_HRD_MMR_status", 10, 6)
  savefig(pg, "Figure3_Platinum_HRD_MMR_genes", 10, max(6, .28 * nrow(gene_plot) + 2.5))
  savefig(cowplot::plot_grid(ps, pg, ncol = 1, rel_heights = c(.38, .62)),
    "Figure3_Platinum_HRD_MMR", 12, 15)

  # Signature analysis is delegated to a local, explicitly configured Python
  # environment. No packages or reference data are installed by this script.
  signature_status <- "Disabled by configuration"
  if (isTRUE(config$signature_analysis)) {
    if (is.null(config$signature_python) || !file.exists(config$signature_python)) {
      stop("signature_analysis=TRUE but signature_python does not exist.")
    }
    if (is.null(config$signature_reference_volume) || !dir.exists(config$signature_reference_volume)) {
      stop("signature_analysis=TRUE but signature_reference_volume does not exist.")
    }
    cohort_file <- file.path(out, "private", "signature", "signature_cohort.tsv")
    readr::write_tsv(dplyr::transmute(d, sample_id = Tumor_Sample_Barcode,
      patient_id, Preop_platinum_exposure), cohort_file, na = "NA")
    args <- c("scripts/mutational_signature_analysis.py",
      "--maf", config$maf_file, "--cohort", cohort_file,
      "--clinical", config$clinical_file, "--clinical-sheet", config$clinical_sheet,
      "--output", file.path(out, "private", "signature"),
      "--public-output", out, "--reference-volume", config$signature_reference_volume,
      "--min-snv", config$signature_min_snv, "--bootstraps", config$signature_bootstraps,
      "--min-count", config$signature_positive_min_count,
      "--min-proportion", config$signature_positive_min_proportion,
      "--min-cosine", config$signature_min_cosine,
      "--min-stability", config$signature_min_stability,
      "--seed", ifelse(is.null(config$random_seed), 20260928L, config$random_seed))
    log_file <- file.path(out, "private", "signature", "signature_analysis.log")
    exit_status <- system2(config$signature_python, args = args, stdout = log_file, stderr = log_file)
    if (!identical(exit_status, 0L)) stop("Mutational signature analysis failed. See ", log_file)
    signature_status <- "Completed"
  }

  invisible(list(distant_only = distant_result, expanded = expanded_result,
    platinum_repair_sets = status_results, platinum_repair_genes = gene_results,
    signature_status = signature_status))
}
