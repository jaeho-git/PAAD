#!/usr/bin/env Rscript

# 목적: PDAC 분석 코호트의 임상 변수별 기술통계와 빈도표를 Excel로 저장합니다.
# 입력: config 파일에 지정된 임상 XLSX, 대상 환자 TXT, MAF 파일
# 출력: 2_Clinical_Table_PDAC.xlsx
# 실행(저장소 루트에서):
#   Rscript --vanilla scripts/clinical_summary.R --config=config/local.R
# 합성 예제 실행:
#   Rscript --vanilla scripts/clinical_summary.R --config=config/config.synthetic.R

# 이 스크립트는 경로 해석을 단순하게 유지하기 위해 반드시 저장소 루트에서 실행합니다.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c("R/config.R", "R/data.R")
if (!all(file.exists(file.path(PROJECT_ROOT, required_project_files)))) {
  stop(
    "Run this script from the PAAD repository root.\n",
    "Example: Rscript --vanilla scripts/clinical_summary.R --config=config/local.R",
    call. = FALSE
  )
}

# 공통 설정 읽기와 입력 전처리 helper만 재사용합니다.
source(file.path(PROJECT_ROOT, "R", "config.R"))
source(file.path(PROJECT_ROOT, "R", "data.R"))

write_clinical_summary <- function(clinical, variables, filename) {
  summary_sheets <- lapply(variables, function(variable) {
    if (is.numeric(clinical[[variable]])) {
      clinical |>
        dplyr::summarise(
          Var = variable,
          Mean = mean(.data[[variable]], na.rm = TRUE),
          Median = stats::median(.data[[variable]], na.rm = TRUE),
          SD = stats::sd(.data[[variable]], na.rm = TRUE)
        )
    } else {
      clinical |>
        dplyr::count(.data[[variable]]) |>
        dplyr::mutate(
          Pct = n / sum(n) * 100,
          Var = variable
        )
    }
  })
  names(summary_sheets) <- clinical_label(variables)

  writexl::write_xlsx(summary_sheets, filename)
  invisible(filename)
}

run_clinical_summary <- function(data, config) {
  # These are the variables included in the clinical workbook. Keeping the
  # list here makes the table definition visible in the analysis file itself.
  variables <- intersect(
    c(
      "Sex", "Age_group", "Differentiation", "NAC", "T", "N", "N_status",
      "BMI", "CA19_9", "CEA", "Stage", "Stage_Group", "Size"
    ),
    names(data$clinical)
  )
  if (identical(config$clinical_schema, "curated_v3")) variables <- c(curated_continuous(), curated_categorical(), "Age_group", "Stage_Group")
  variables <- intersect(variables, names(data$clinical))
  output_file <- output_path(config, "2_Clinical_Table_PDAC.xlsx")

  message("Writing clinical summary for ", nrow(data$clinical), " PDAC samples.")
  write_clinical_summary(
    clinical = data$clinical,
    variables = variables,
    filename = output_file
  )
  message("Saved: ", output_file)

  invisible(output_file)
}

# ---- 명시적 실행부 ---------------------------------------------------------
config <- load_config_from_command_line(PROJECT_ROOT)
analysis_data <- load_analysis_data(config)
run_clinical_summary(analysis_data, config)
