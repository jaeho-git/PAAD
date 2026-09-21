#!/usr/bin/env Rscript

# 목적: PDAC 변이 유무 heatmap(top 20/top 5)과 샘플별 retained MAF row count heatmap을 만듭니다.
# 입력: config 파일에 지정된 임상 XLSX, 대상 환자 TXT, MAF 파일
# 출력: 5_Heatmap_Top20_PDAC.png, 6_Heatmap_Top5_PDAC.png, 7_TMB_Heatmap_PDAC.png
# 실행(저장소 루트에서):
#   Rscript --vanilla scripts/mutation_heatmaps.R --config=config/local.R
# 합성 예제 실행:
#   Rscript --vanilla scripts/mutation_heatmaps.R --config=config/config.synthetic.R
# 주의: 이 프로젝트의 기존 "TMB" 값은 mutations/Mb가 아니라 샘플별 retained MAF 행 수입니다.

# 이 스크립트는 경로 해석을 단순하게 유지하기 위해 반드시 저장소 루트에서 실행합니다.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c("R/config.R", "R/data.R", "R/plot_helpers.R")
if (!all(file.exists(file.path(PROJECT_ROOT, required_project_files)))) {
  stop(
    "Run this script from the PAAD repository root.\n",
    "Example: Rscript --vanilla scripts/mutation_heatmaps.R --config=config/local.R",
    call. = FALSE
  )
}

# 공통 설정 읽기, 입력 전처리, annotation/PNG helper만 재사용합니다.
source(file.path(PROJECT_ROOT, "R", "config.R"))
source(file.path(PROJECT_ROOT, "R", "data.R"))
source(file.path(PROJECT_ROOT, "R", "plot_helpers.R"))

draw_mutation_heatmap <- function(maf, top_n, title, filename, annotation_columns, dpi) {
  top_genes <- maftools::getGeneSummary(maf) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Hugo_Symbol)

  mutation_matrix <- maf@data |>
    dplyr::filter(Hugo_Symbol %in% top_genes) |>
    dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) |>
    dplyr::distinct() |>
    dplyr::mutate(value = 1L) |>
    tidyr::pivot_wider(
      names_from = Tumor_Sample_Barcode,
      values_from = value,
      values_fill = 0L
    ) |>
    tibble::column_to_rownames("Hugo_Symbol") |>
    as.matrix()
  if (!nrow(mutation_matrix) || !ncol(mutation_matrix)) {
    stop("Clustered heatmap requires mutation data.", call. = FALSE)
  }

  clinical <- as.data.frame(maf@clinical.data) |>
    dplyr::filter(Tumor_Sample_Barcode %in% colnames(mutation_matrix)) |>
    tibble::column_to_rownames("Tumor_Sample_Barcode")
  clinical <- clinical[colnames(mutation_matrix), , drop = FALSE]
  top_annotation <- make_annotation(clinical, annotation_columns)

  open_png(filename, width = 14, height = 9, dpi = dpi)
  on.exit(close_device(), add = TRUE)
  heatmap <- ComplexHeatmap::Heatmap(
    mutation_matrix,
    name = "Mutation",
    col = c("0" = "white", "1" = "black"),
    top_annotation = top_annotation,
    column_title = paste0(title, " (Top ", top_n, ", Total N=", ncol(mutation_matrix), ")"),
    cluster_rows = nrow(mutation_matrix) > 1L,
    cluster_columns = ncol(mutation_matrix) > 1L,
    show_column_names = FALSE,
    rect_gp = grid::gpar(col = "grey90", lwd = 0.5)
  )
  ComplexHeatmap::draw(heatmap, merge_legend = TRUE)

  invisible(filename)
}

draw_retained_maf_count_heatmap <- function(maf, title, filename, annotation_columns, dpi) {
  mutation_counts <- maf@data |>
    dplyr::count(Tumor_Sample_Barcode, name = "TMB")
  clinical <- as.data.frame(maf@clinical.data) |>
    dplyr::left_join(mutation_counts, by = "Tumor_Sample_Barcode") |>
    dplyr::mutate(TMB = tidyr::replace_na(TMB, 0L)) |>
    dplyr::arrange(dplyr::desc(TMB))

  count_matrix <- matrix(
    clinical$TMB,
    nrow = 1L,
    dimnames = list("TMB", clinical$Tumor_Sample_Barcode)
  )
  annotation_data <- clinical |>
    tibble::column_to_rownames("Tumor_Sample_Barcode")
  top_annotation <- make_annotation(annotation_data, annotation_columns)

  color_limits <- range(count_matrix, na.rm = TRUE)
  if (color_limits[[1]] == color_limits[[2]]) {
    color_limits <- color_limits + c(-0.5, 0.5)
  }

  open_png(filename, width = 14, height = 7, dpi = dpi)
  on.exit(close_device(), add = TRUE)
  heatmap <- ComplexHeatmap::Heatmap(
    count_matrix,
    name = "TMB",
    col = circlize::colorRamp2(color_limits, c("#F7FBFF", "#08306B")),
    top_annotation = top_annotation,
    column_title = paste0(title, " (N=", ncol(count_matrix), ")"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    rect_gp = grid::gpar(col = "white", lwd = 0.5)
  )
  ComplexHeatmap::draw(heatmap, merge_legend = TRUE)

  invisible(filename)
}

run_mutation_heatmaps <- function(data, config) {
  # Clinical tracks shown above each heatmap, in their intended display order.
  annotation_columns <- intersect(
    c(
      "Sex", "Age_group", "Differentiation", "NAC", "T", "N", "N_status",
      "BMI", "CA19_9", "CEA", "Stage", "Stage_Group", "Size"
    ),
    names(data$maf@clinical.data)
  )

  message("[1/3] Top-", config$top_n, " mutation heatmap")
  draw_mutation_heatmap(
    maf = data$maf,
    top_n = config$top_n,
    # Keep the v19 publication title and filename; production top_n is 20.
    title = "PDAC Clustered Heatmap (Top 20)",
    filename = output_path(config, "5_Heatmap_Top20_PDAC.png"),
    annotation_columns = annotation_columns,
    dpi = config$plot_dpi
  )

  message("[2/3] Top-", min(5L, config$top_n), " mutation heatmap")
  draw_mutation_heatmap(
    maf = data$maf,
    top_n = min(5L, config$top_n),
    title = "PDAC Clustered Heatmap (Top 5)",
    filename = output_path(config, "6_Heatmap_Top5_PDAC.png"),
    annotation_columns = annotation_columns,
    dpi = config$plot_dpi
  )

  message("[3/3] Retained MAF row count heatmap")
  draw_retained_maf_count_heatmap(
    maf = data$maf,
    title = "TMB Heatmap (PDAC)",
    filename = output_path(config, "7_TMB_Heatmap_PDAC.png"),
    annotation_columns = annotation_columns,
    dpi = config$plot_dpi
  )

  invisible(TRUE)
}

# ---- 명시적 실행부 ---------------------------------------------------------
config <- load_config_from_command_line(PROJECT_ROOT)
analysis_data <- load_analysis_data(config)
run_mutation_heatmaps(analysis_data, config)
