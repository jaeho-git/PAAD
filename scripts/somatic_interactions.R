#!/usr/bin/env Rscript

# 목적: 상위 변이 유전자 사이의 공발생(co-occurrence) 및 상호배타 관계를 분석합니다.
# 입력: config 파일에 지정된 임상 XLSX, 대상 환자 TXT, MAF 파일
# 출력: 9_Somatic_Interactions_PDAC.png, 9_Somatic_Interactions_Results_PDAC.xlsx
# 실행(저장소 루트에서):
#   Rscript --vanilla scripts/somatic_interactions.R --config=config/local.R
# 합성 예제 실행:
#   Rscript --vanilla scripts/somatic_interactions.R --config=config/config.synthetic.R

# 이 스크립트는 경로 해석을 단순하게 유지하기 위해 반드시 저장소 루트에서 실행합니다.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c("R/config.R", "R/data.R", "R/plot_helpers.R")
if (!all(file.exists(file.path(PROJECT_ROOT, required_project_files)))) {
  stop(
    "Run this script from the PAAD repository root.\n",
    "Example: Rscript --vanilla scripts/somatic_interactions.R --config=config/local.R",
    call. = FALSE
  )
}

# 공통 설정 읽기, 입력 전처리, PNG helper만 재사용합니다.
source(file.path(PROJECT_ROOT, "R", "config.R"))
source(file.path(PROJECT_ROOT, "R", "data.R"))
source(file.path(PROJECT_ROOT, "R", "plot_helpers.R"))

run_somatic_interactions <- function(data, config) {
  figure_file <- output_path(config, "9_Somatic_Interactions_PDAC.png")
  table_file <- output_path(config, "9_Somatic_Interactions_Results_PDAC.xlsx")

  message("Analyzing somatic interactions among the top 25 mutated genes.")
  open_png(figure_file, width = 12, height = 10, dpi = config$plot_dpi)
  on.exit(close_device(), add = TRUE)
  interaction_results <- maftools::somaticInteractions(
    maf = data$maf,
    top = 25,
    pvalue = c(0.05, 0.01),
    fontSize = 0.6
  )
  close_device()
  on.exit(NULL, add = FALSE)

  writexl::write_xlsx(as.data.frame(interaction_results), table_file)
  message("Saved: ", figure_file)
  message("Saved: ", table_file)

  invisible(c(figure_file, table_file))
}

# ---- 명시적 실행부 ---------------------------------------------------------
config <- load_config_from_command_line(PROJECT_ROOT)
analysis_data <- load_analysis_data(config)
run_somatic_interactions(analysis_data, config)
