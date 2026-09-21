#!/usr/bin/env Rscript

# 목적: PDAC 코호트의 주요 변이 유전자를 두 개의 oncoplot으로 시각화합니다.
# 입력: config 파일에 지정된 임상 XLSX, 대상 환자 TXT, MAF 파일
# 출력: 1_Oncoplot_Main_PDAC.png, 1_2_Oncoplot_TN_Filtered_PDAC.png
# 실행(저장소 루트에서):
#   Rscript --vanilla scripts/oncoplots.R --config=config/local.R
# 합성 예제 실행:
#   Rscript --vanilla scripts/oncoplots.R --config=config/config.synthetic.R

# 이 스크립트는 경로 해석을 단순하게 유지하기 위해 반드시 저장소 루트에서 실행합니다.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c("R/config.R", "R/data.R", "R/plot_helpers.R")
if (!all(file.exists(file.path(PROJECT_ROOT, required_project_files)))) {
  stop(
    "Run this script from the PAAD repository root.\n",
    "Example: Rscript --vanilla scripts/oncoplots.R --config=config/local.R",
    call. = FALSE
  )
}

# 공통 설정 읽기, 입력 전처리, 색상/PNG helper만 재사용합니다.
source(file.path(PROJECT_ROOT, "R", "config.R"))
source(file.path(PROJECT_ROOT, "R", "data.R"))
source(file.path(PROJECT_ROOT, "R", "plot_helpers.R"))

draw_oncoplot <- function(maf, top_n, title, filename, annotation_columns, dpi) {
  top_genes <- maftools::getGeneSummary(maf) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Hugo_Symbol)
  if (!length(top_genes)) stop("Oncoplot requires at least one mutation.", call. = FALSE)

  mutation_long <- maf@data |>
    dplyr::filter(Hugo_Symbol %in% top_genes) |>
    dplyr::transmute(
      Hugo_Symbol,
      Tumor_Sample_Barcode,
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

  mutation_matrix <- mutation_long |>
    tidyr::pivot_wider(
      names_from = Tumor_Sample_Barcode,
      values_from = Type,
      values_fn = function(x) paste(x, collapse = ";")
    ) |>
    tibble::column_to_rownames("Hugo_Symbol") |>
    as.matrix()

  clinical <- as.data.frame(maf@clinical.data) |>
    dplyr::filter(Tumor_Sample_Barcode %in% colnames(mutation_matrix)) |>
    tibble::column_to_rownames("Tumor_Sample_Barcode")
  clinical <- clinical[colnames(mutation_matrix), , drop = FALSE]

  cohort_text <- ""
  if ("Classification" %in% names(clinical)) {
    cohort_counts <- table(clinical$Classification)
    cohort_text <- paste(names(cohort_counts), as.integer(cohort_counts), sep = ":", collapse = ", ")
  }
  plot_title <- paste0(
    title,
    "\n(Total N=", ncol(mutation_matrix),
    if (nzchar(cohort_text)) paste0(" | ", cohort_text) else "",
    ")"
  )

  top_annotation <- make_annotation(clinical, annotation_columns)
  mutation_colors <- c(
    Frameshift = "#FF7F50", Nonsense = "#DC143C", Missense = "#3CB371",
    Splice_Site = "#9370DB", In_Frame = "#FFD700", Multi_Hit = "#000000"
  )
  alter_functions <- c(
    list(background = function(x, y, w, h) {
      grid::grid.rect(x, y, w * 0.9, h * 0.9, gp = grid::gpar(fill = "#FAFAFA", col = NA))
    }),
    lapply(names(mutation_colors), function(type) {
      force(type)
      function(x, y, w, h) {
        grid::grid.rect(
          x, y, w * 0.9, h * 0.4,
          gp = grid::gpar(fill = mutation_colors[[type]], col = NA)
        )
      }
    })
  )
  names(alter_functions)[-1] <- names(mutation_colors)

  open_png(filename, width = 14, height = 9, dpi = dpi)
  on.exit(close_device(), add = TRUE)
  oncoplot <- ComplexHeatmap::oncoPrint(
    mutation_matrix,
    alter_fun = alter_functions,
    col = mutation_colors,
    top_annotation = top_annotation,
    column_title = plot_title,
    column_title_gp = grid::gpar(fontsize = 10),
    pct_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 9),
    show_column_names = FALSE
  )
  ComplexHeatmap::draw(oncoplot, merge_legend = TRUE)
  invisible(filename)
}

run_oncoplots <- function(data, config) {
  # Clinical tracks shown above the oncoplot, in their intended display order.
  annotation_columns <- intersect(
    c(
      "Sex", "Age_group", "Differentiation", "NAC", "T", "N", "N_status",
      "BMI", "CA19_9", "CEA", "Stage", "Stage_Group", "Size"
    ),
    names(data$maf@clinical.data)
  )

  message("[1/2] Main PDAC oncoplot")
  draw_oncoplot(
    maf = data$maf,
    top_n = config$top_n,
    # Keep the publication title used by v19; production top_n is 20.
    title = "Top 20 Mutated Genes in PDAC",
    filename = output_path(config, "1_Oncoplot_Main_PDAC.png"),
    annotation_columns = annotation_columns,
    dpi = config$plot_dpi
  )

  complete_tn_samples <- data$clinical |>
    dplyr::filter(!is.na(T), !is.na(N)) |>
    dplyr::pull(Tumor_Sample_Barcode)

  if (length(complete_tn_samples)) {
    message("[2/2] T/N-filtered PDAC oncoplot")
    tn_maf <- maftools::subsetMaf(data$maf, tsb = complete_tn_samples)
    draw_oncoplot(
      maf = tn_maf,
      top_n = config$top_n,
      title = "PDAC Oncoplot (Filtered by T & N)",
      filename = output_path(config, "1_2_Oncoplot_TN_Filtered_PDAC.png"),
      annotation_columns = annotation_columns,
      dpi = config$plot_dpi
    )
  } else {
    message("Skipping T/N-filtered oncoplot: no samples have both T and N.")
  }

  invisible(TRUE)
}

# ---- 명시적 실행부 ---------------------------------------------------------
config <- load_config_from_command_line(PROJECT_ROOT)
analysis_data <- load_analysis_data(config)
run_oncoplots(analysis_data, config)
