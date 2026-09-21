#!/usr/bin/env Rscript

# 목적: PDAC 분석 코호트의 병기(Stage)별 환자 수를 막대그래프로 표시합니다.
# 입력: config 파일에 지정된 임상 XLSX, 대상 환자 TXT, MAF 파일
# 출력: 10_Stage_Counts_PDAC.png
# 실행(저장소 루트에서):
#   Rscript --vanilla scripts/stage_distribution.R --config=config/local.R
# 합성 예제 실행:
#   Rscript --vanilla scripts/stage_distribution.R --config=config/config.synthetic.R

# 이 스크립트는 경로 해석을 단순하게 유지하기 위해 반드시 저장소 루트에서 실행합니다.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c("R/config.R", "R/data.R")
if (!all(file.exists(file.path(PROJECT_ROOT, required_project_files)))) {
  stop(
    "Run this script from the PAAD repository root.\n",
    "Example: Rscript --vanilla scripts/stage_distribution.R --config=config/local.R",
    call. = FALSE
  )
}

# 공통 설정 읽기와 입력 전처리 helper만 재사용합니다.
source(file.path(PROJECT_ROOT, "R", "config.R"))
source(file.path(PROJECT_ROOT, "R", "data.R"))

run_stage_distribution <- function(data, config) {
  output_file <- output_path(config, "10_Stage_Counts_PDAC.png")
  observed_stages <- sort(unique(stats::na.omit(as.character(data$clinical$Stage))))
  expected_stage_order <- c("0", "IA", "IB", "IIA", "IIB", "III", "IV", "I", "II")
  stage_levels <- unique(c(
    expected_stage_order[expected_stage_order %in% observed_stages],
    observed_stages
  ))
  stage_counts <- data$clinical |>
    dplyr::filter(!is.na(Stage)) |>
    dplyr::mutate(Stage = factor(as.character(Stage), levels = stage_levels)) |>
    dplyr::count(Stage, .drop = FALSE) |>
    dplyr::mutate(n = tidyr::replace_na(n, 0L))

  stage_plot <- ggplot2::ggplot(
    stage_counts,
    ggplot2::aes(x = Stage, y = n, fill = Stage)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.5) +
    ggplot2::scale_fill_brewer(palette = "Spectral", drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(title = "PDAC Patient Count by Stage", y = "Count") +
    ggplot2::theme_classic()

  ggplot2::ggsave(
    output_file,
    stage_plot,
    width = 8,
    height = 6,
    dpi = config$plot_dpi,
    bg = "white"
  )
  message("Saved: ", output_file)

  invisible(output_file)
}

# ---- 명시적 실행부 ---------------------------------------------------------
config <- load_config_from_command_line(PROJECT_ROOT)
analysis_data <- load_analysis_data(config)
run_stage_distribution(analysis_data, config)
