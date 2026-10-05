#!/usr/bin/env Rscript

# 목적: 임상 그룹별 유전자 변이 빈도와 retained MAF row count("MAF_variant_count")를 비교합니다.
# 입력: config 파일에 지정된 임상 XLSX, 대상 환자 TXT, MAF 파일
# 출력: 각 임상 변수별 4_FreqPlot_PDAC_<변수>.png 및 4_TMB_BoxPlot_PDAC_<변수>.png
# 실행(저장소 루트에서):
#   Rscript --vanilla scripts/clinical_group_comparisons.R --config=config/local.R
# 합성 예제 실행:
#   Rscript --vanilla scripts/clinical_group_comparisons.R --config=config/config.synthetic.R
# 주의: 이 프로젝트의 기존 "MAF_variant_count" 값은 mutations/Mb가 아니라 샘플별 retained MAF 행 수입니다.

# 이 스크립트는 경로 해석을 단순하게 유지하기 위해 반드시 저장소 루트에서 실행합니다.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c("R/config.R", "R/data.R")
if (!all(file.exists(file.path(PROJECT_ROOT, required_project_files)))) {
  stop(
    "Run this script from the PAAD repository root.\n",
    "Example: Rscript --vanilla scripts/clinical_group_comparisons.R --config=config/local.R",
    call. = FALSE
  )
}

# 공통 설정 읽기와 입력 전처리 helper만 재사용합니다.
source(file.path(PROJECT_ROOT, "R", "config.R"))
source(file.path(PROJECT_ROOT, "R", "data.R"))
source(file.path(PROJECT_ROOT, "R", "pairwise_tests.R"))

comparison_level_order <- function(variable, observed) {
  # Group order is part of this analysis definition and is therefore kept next
  # to the comparison code instead of being hidden in a shared helper file.
  preset_levels <- list(
    Differentiation = c("WD", "MD", "PD"),
    NAC = c("No", "Yes"),
    T = c("T0", "T1", "T2", "T3", "T4"),
    N = c("N0", "N1", "N2"),
    N_status = c("N0", "N1"),
    Stage_Group = c("I", "II", "III")
  )
  observed <- sort(unique(stats::na.omit(as.character(observed))))

  levels <- if (variable %in% names(preset_levels)) {
    unique(c(preset_levels[[variable]], observed))
  } else if (variable == "Stage") {
    stage_order <- c("0", "IA", "IB", "IIA", "IIB", "III", "IV", "I", "II")
    unique(c(stage_order[stage_order %in% observed], observed))
  } else {
    observed
  }

  levels[!is.na(levels) & nzchar(levels)]
}

compare_one_clinical_variable <- function(
    maf,
    clinical,
    variable,
    top_genes,
    output_dir,
    output_prefix,
    dpi) {
  clinical_grouped <- clinical |>
    dplyr::filter(!is.na(.data[[variable]])) |>
    dplyr::mutate("{variable}" := as.character(.data[[variable]]))

  all_levels <- comparison_level_order(variable, clinical_grouped[[variable]])
  observed_levels <- intersect(all_levels, unique(clinical_grouped[[variable]]))
  if (!length(observed_levels)) return(character())

  clinical_grouped <- clinical_grouped |>
    dplyr::mutate("{variable}" := factor(.data[[variable]], levels = all_levels))
  group_sizes <- table(factor(clinical_grouped[[variable]], levels = all_levels))
  cohort_subtitle <- paste0(
    "Total N=", sum(group_sizes), " (",
    paste(names(group_sizes), as.integer(group_sizes), sep = "=", collapse = ", "),
    ")"
  )

  top_gene_mutations <- maf@data |>
    dplyr::filter(Hugo_Symbol %in% top_genes) |>
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) |>
    dplyr::distinct()
  mutations_with_group <- dplyr::inner_join(
    top_gene_mutations,
    clinical_grouped,
    by = "Tumor_Sample_Barcode"
  )

  # 각 유전자에서 임상 그룹 간 mutated/wild-type 비율을 Fisher 검정으로 비교합니다.
  significance <- tibble::tibble()
  if (length(observed_levels) >= 2L) {
    significance <- purrr::map_dfr(top_genes, function(gene) {
      mutated_samples <- mutations_with_group |>
        dplyr::filter(Hugo_Symbol == gene) |>
        dplyr::pull(Tumor_Sample_Barcode)
      mutated_counts <- vapply(observed_levels, function(group) {
        sum(
          as.character(clinical_grouped[[variable]]) == group &
            clinical_grouped$Tumor_Sample_Barcode %in% mutated_samples
        )
      }, numeric(1))
      wild_type_counts <- as.numeric(group_sizes[observed_levels]) - mutated_counts
      p_value <- tryCatch(
        stats::fisher.test(
          rbind(mutated_counts, wild_type_counts),
          workspace = 2e5,
          simulate.p.value = TRUE
        )$p.value,
        error = function(error) NA_real_
      )
      tibble::tibble(
        Hugo_Symbol = gene,
        Significance = dplyr::case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
    })
  }

  mutation_frequency <- mutations_with_group |>
    dplyr::group_by(.data[[variable]], Hugo_Symbol) |>
    dplyr::summarise(Mutated = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") |>
    tidyr::complete(
      !!rlang::sym(variable) := factor(all_levels, levels = all_levels),
      Hugo_Symbol = top_genes,
      fill = list(Mutated = 0L)
    ) |>
    dplyr::mutate(
      Total = as.numeric(group_sizes[as.character(.data[[variable]])]),
      Frequency = ifelse(Total > 0, Mutated / Total * 100, 0),
      Hugo_Symbol = factor(Hugo_Symbol, levels = top_genes),
      "{variable}" := factor(as.character(.data[[variable]]), levels = all_levels)
    )

  frequency_plot <- ggplot2::ggplot(
    mutation_frequency,
    ggplot2::aes(x = Hugo_Symbol, y = Frequency, fill = .data[[variable]])
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.9), width = 0.8) +
    ggplot2::scale_fill_brewer(palette = "Pastel1", drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      title = paste0("[", output_prefix, "] Variant frequency by ", figure_label(variable)),
      subtitle = cohort_subtitle,
      y = "Frequency (%)",
      x = "Gene", fill = figure_label(variable)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold")
    )

  if (nrow(significance)) {
    maximum_frequency <- mutation_frequency |>
      dplyr::group_by(Hugo_Symbol) |>
      dplyr::summarise(Max = max(Frequency), .groups = "drop")
    frequency_plot <- frequency_plot +
      ggplot2::geom_text(
        data = dplyr::inner_join(significance, maximum_frequency, by = "Hugo_Symbol"),
        ggplot2::aes(x = Hugo_Symbol, y = Max + 3, label = Significance),
        inherit.aes = FALSE,
        size = 3,
        color = "red"
      )
  }

  if (identical(config$clinical_schema, "curated_v3")) {
    pw <- dplyr::bind_rows(lapply(top_genes, function(gene) {
      calls <- clinical_grouped$Tumor_Sample_Barcode %in% top_gene_mutations$Tumor_Sample_Barcode[top_gene_mutations$Hugo_Symbol == gene]
      p <- pairwise_categorical(calls, clinical_grouped[[variable]])
      if (nrow(p)) cbind(data.frame(gene, variable), p) else NULL
    }))
    readr::write_tsv(pw, file.path(output_dir, paste0("4_Frequency_pairwise_", variable, ".tsv")))
    frequency_plot <- frequency_plot + ggplot2::labs(caption = paste0("Stars: nominal omnibus Fisher p. All pairs with Holm p: 4_Frequency_pairwise_", variable, ".tsv"))
  }
  frequency_file <- file.path(
    output_dir,
    paste0("4_FreqPlot_", output_prefix, "_", variable, ".png")
  )
  ggplot2::ggsave(
    frequency_file,
    frequency_plot,
    width = 14,
    height = 7,
    dpi = dpi,
    bg = "white"
  )

  # 기존 v19 분석과 같은 정의: MAF에서 유지된 행을 샘플별로 셉니다.
  mutation_counts <- maf@data |>
    dplyr::count(Tumor_Sample_Barcode, name = "MAF_variant_count") |>
    dplyr::right_join(clinical_grouped, by = "Tumor_Sample_Barcode") |>
    dplyr::mutate(
      MAF_variant_count = tidyr::replace_na(MAF_variant_count, 0L),
      "{variable}" := factor(as.character(.data[[variable]]), levels = all_levels)
    )

  test_label <- "Statistical test not performed"
  if (length(observed_levels) >= 2L) {
    formula <- stats::as.formula(paste0("MAF_variant_count ~ `", variable, "`"))
    test <- tryCatch(
      if (length(observed_levels) == 2L) {
        stats::wilcox.test(
          formula,
          data = mutation_counts |> dplyr::filter(.data[[variable]] %in% observed_levels)
        )
      } else {
        stats::kruskal.test(
          formula,
          data = mutation_counts |> dplyr::filter(.data[[variable]] %in% observed_levels)
        )
      },
      error = function(error) NULL
    )
    if (!is.null(test)) {
      test_name <- if (length(observed_levels) == 2L) "Wilcoxon" else "Kruskal"
      test_label <- paste0(test_name, " P: ", sprintf("%.4f", test$p.value))
    }
  }

  tmb_plot <- ggplot2::ggplot(
    mutation_counts,
    ggplot2::aes(x = .data[[variable]], y = MAF_variant_count, fill = .data[[variable]])
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.65) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    ggplot2::scale_fill_brewer(palette = "Pastel1", drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      title = paste0("[", output_prefix, "] Variant count by ", figure_label(variable)),
      subtitle = paste0(cohort_subtitle, "\n", test_label),
      y = "Variant count", x = figure_label(variable), fill = figure_label(variable)
    ) +
    ggplot2::theme_classic()

  if (identical(config$clinical_schema, "curated_v3")) {
    pw <- pairwise_numeric(mutation_counts$MAF_variant_count, mutation_counts[[variable]])
    readr::write_tsv(pw, file.path(output_dir, paste0("4_Variant_count_pairwise_", variable, ".tsv")))
    tmb_plot <- tmb_plot + ggplot2::labs(caption = pairwise_caption(pw, 75)) +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 8, hjust = 0))
  }
  tmb_file <- file.path(
    output_dir,
    paste0("4_TMB_BoxPlot_", output_prefix, "_", variable, ".png")
  )
  ggplot2::ggsave(
    tmb_file,
    tmb_plot,
    width = 9,
    height = 8,
    dpi = dpi,
    bg = "white"
  )

  c(frequency_file, tmb_file)
}

run_clinical_group_comparisons <- function(data, config) {
  # random_seed가 있으면 개별 실행과 전체 실행의 simulated Fisher 결과를 맞춥니다.
  if (is.null(config$random_seed)) {
    message(
      "Compatibility note: random_seed is NULL, so simulated Fisher tests ",
      "remain unseeded as in v19."
    )
  } else {
    set.seed(config$random_seed)
  }

  comparison_variables <- intersect(
    c("Differentiation", "NAC", "T", "N", "N_status", "Stage", "Stage_Group"),
    names(data$clinical)
  )
  if (identical(config$clinical_schema, "curated_v3")) comparison_variables <- c(
    "Differentiation", "Neoadjuvant", "T_stage", "N_stage", "LN_positive", "AJCC_stage", "Stage_Group")
  top_genes <- maftools::getGeneSummary(data$maf) |>
    dplyr::slice_head(n = config$top_n) |>
    dplyr::pull(Hugo_Symbol)

  outputs <- unlist(lapply(comparison_variables, function(variable) {
    message("Comparing mutation patterns by ", variable)
    compare_one_clinical_variable(
      maf = data$maf,
      clinical = data$clinical,
      variable = variable,
      top_genes = top_genes,
      output_dir = config$output_dir,
      output_prefix = "PDAC",
      dpi = config$plot_dpi
    )
  }), use.names = FALSE)

  invisible(outputs)
}

# ---- 명시적 실행부 ---------------------------------------------------------
config <- load_config_from_command_line(PROJECT_ROOT)
analysis_data <- load_analysis_data(config)
run_clinical_group_comparisons(analysis_data, config)
