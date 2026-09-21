# Shared input preparation only.
#
# This file does not define any figure, statistical comparison, driver rule, or
# survival analysis. It performs the steps every analysis needs in exactly the
# same way: package checks, reading the three inputs, v19-compatible cohort
# filtering/column normalization, and MAF construction. Analysis-specific
# variables and rules are kept in the descriptively named files under scripts/.

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         "\nInstall packages before running the analysis; scripts never install packages automatically.",
         call. = FALSE)
  }
  invisible(TRUE)
}

analysis_packages <- function() {
  c("maftools", "readxl", "writexl", "tidyverse", "ComplexHeatmap",
    "circlize", "RColorBrewer", "grid", "survival")
}

load_analysis_packages <- function() {
  packages <- analysis_packages()
  assert_packages(packages)
  suppressPackageStartupMessages({
    for (package in packages) library(package, character.only = TRUE)
  })
  invisible(packages)
}

rename_if_present <- function(data, old, new) {
  if (old %in% names(data) && !new %in% names(data)) names(data)[names(data) == old] <- new
  data
}

prepare_clinical_data <- function(clinical_data, target_patients, sex_column_policy = "legacy") {
  # Reproduce the v19 cohort construction in one place so every standalone
  # analysis receives the same patients, column names, and derived categories.
  clinical_data <- rename_if_present(clinical_data, "CA 19-9", "CA19_9")
  clinical_data <- rename_if_present(clinical_data, "N status", "N_status")

  if (!"Tumor_Sample_Barcode" %in% names(clinical_data) ||
      !"Tumor_Sample_Barcode" %in% names(target_patients)) {
    stop("Both clinical and target-patient inputs require Tumor_Sample_Barcode.", call. = FALSE)
  }
  if (sex_column_policy == "clinical_korean") {
    if (!"성별코드" %in% names(clinical_data)) stop("sex_column_policy='clinical_korean' requires 성별코드.", call. = FALSE)
    clinical_data$Sex <- as.character(clinical_data[["성별코드"]])
  }

  duplicate_target_columns <- c(
    "ID_NGS", "OS_m", "RFS_m", "OS_d", "RFS_d", "survive", "Recur",
    "LVI", "PNI", "RM", "age", "Age"
  )
  clinical_base <- clinical_data |>
    dplyr::filter(Tumor_Sample_Barcode %in% target_patients$Tumor_Sample_Barcode) |>
    dplyr::select(-dplyr::any_of(duplicate_target_columns))

  join_columns <- c(
    "Tumor_Sample_Barcode", "ID_NGS", "age", "OS_m", "RFS_m", "OS_d",
    "RFS_d", "survive", "Recur", "LVI", "PNI", "RM"
  )
  if (sex_column_policy == "target_patients") {
    if (!"sex" %in% names(target_patients)) stop("sex_column_policy='target_patients' requires sex.", call. = FALSE)
    join_columns <- c(join_columns, "sex")
  }
  target_join <- target_patients |>
    dplyr::mutate(Tumor_Sample_Barcode = stringr::str_trim(as.character(Tumor_Sample_Barcode))) |>
    dplyr::select(dplyr::any_of(join_columns)) |>
    dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
  if (sex_column_policy == "target_patients") names(target_join)[names(target_join) == "sex"] <- "Sex"

  clinical_joined <- dplyr::left_join(clinical_base, target_join, by = "Tumor_Sample_Barcode")
  if (!"age" %in% names(clinical_joined)) clinical_joined$age <- NA_real_
  for (column in c("LVI", "PNI", "RM")) {
    if (!column %in% names(clinical_joined)) clinical_joined[[column]] <- NA_character_
  }

  clinical_processed <- clinical_joined |>
    dplyr::mutate(
      Age = suppressWarnings(as.numeric(stringr::str_replace_all(as.character(age), "[^0-9.]", ""))),
      Age_group = dplyr::case_when(
        !is.na(Age) & Age <= 55 ~ "Age <=55",
        !is.na(Age) & Age >= 70 ~ "Age >=70",
        !is.na(Age) & Age > 55 & Age < 70 ~ "Age 56-69",
        TRUE ~ NA_character_
      ),
      LVI = as.character(LVI),
      PNI = as.character(PNI),
      RM = as.character(RM)
    )

  required_clinical <- c("Classification", "Differentiation", "NAC", "T", "N", "Stage", "Size", "BMI", "CA19_9", "CEA")
  missing <- setdiff(required_clinical, names(clinical_processed))
  if (length(missing)) stop("Clinical input is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!"N_status" %in% names(clinical_processed)) clinical_processed$N_status <- NA_character_

  normalize_binary_pathology <- function(x) {
    raw <- stringr::str_trim(as.character(x))
    lower <- stringr::str_to_lower(raw)
    dplyr::case_when(
      raw %in% c("1", "1.0") | lower %in% c("positive", "pos", "yes", "y") ~ "Positive",
      raw %in% c("0", "0.0") | lower %in% c("negative", "neg", "no", "n") ~ "Negative",
      TRUE ~ raw
    )
  }

  clinical_processed <- clinical_processed |>
    dplyr::filter(Classification != "Exception") |>
    dplyr::mutate(Tumor_Sample_Barcode = stringr::str_trim(as.character(Tumor_Sample_Barcode))) |>
    dplyr::mutate(
      Differentiation = dplyr::case_when(
        grepl("wd", Differentiation, ignore.case = TRUE) ~ "WD",
        grepl("md", Differentiation, ignore.case = TRUE) ~ "MD",
        grepl("pd", Differentiation, ignore.case = TRUE) ~ "PD",
        TRUE ~ as.character(Differentiation)
      ),
      NAC = dplyr::case_when(NAC == "y" ~ "Yes", NAC == "n" ~ "No", TRUE ~ as.character(NAC)),
      T = dplyr::case_when(as.character(T) %in% c("1", "1.0") ~ "T1", as.character(T) %in% c("2", "2.0") ~ "T2", as.character(T) %in% c("3", "3.0") ~ "T3", TRUE ~ as.character(T)),
      N = dplyr::case_when(as.character(N) %in% c("0", "0.0") ~ "N0", as.character(N) %in% c("1", "1.0") ~ "N1", as.character(N) %in% c("2", "2.0") ~ "N2", TRUE ~ as.character(N)),
      N_status = dplyr::case_when(as.character(N_status) %in% c("0", "0.0") ~ "N0", as.character(N_status) %in% c("1", "1.0") ~ "N1", TRUE ~ as.character(N_status)),
      Stage_Group = dplyr::case_when(Stage %in% c("IA", "IB") ~ "I", Stage %in% c("IIA", "IIB") ~ "II", Stage == "III" ~ "III", TRUE ~ NA_character_),
      Size = suppressWarnings(as.numeric(Size)),
      BMI = suppressWarnings(as.numeric(BMI)),
      CA19_9 = suppressWarnings(as.numeric(stringr::str_replace_all(as.character(CA19_9), "[<>=\\s]", ""))),
      CEA = suppressWarnings(as.numeric(stringr::str_replace_all(as.character(CEA), "[<>=\\s]", ""))),
      dplyr::across(dplyr::any_of(c("OS_m", "RFS_m", "OS_d", "RFS_d", "survive", "Recur")), ~suppressWarnings(as.numeric(.x))),
      LVI = normalize_binary_pathology(LVI),
      PNI = normalize_binary_pathology(PNI),
      RM = normalize_binary_pathology(RM)
    )

  target_cols <- c(
    "Tumor_Sample_Barcode", "Sex", "Classification", "Differentiation", "NAC",
    "Age", "Age_group", "Size", "T", "N", "N_status", "BMI", "CA19_9", "CEA",
    "Stage", "Stage_Group", "OS_m", "RFS_m", "OS_d", "RFS_d", "survive", "Recur",
    "LVI", "PNI", "RM"
  )
  clinical_processed |>
    dplyr::select(dplyr::all_of(intersect(target_cols, names(clinical_processed)))) |>
    dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE) |>
    dplyr::filter(grepl("PDAC", Classification, ignore.case = TRUE))
}

load_analysis_data <- function(config) {
  # Read and connect the clinical workbook, target-patient table, and MAF.
  # The returned list is the complete common input used by scripts/*.R.
  load_analysis_packages()
  clinical_raw <- readxl::read_excel(config$clinical_file, sheet = config$clinical_sheet)
  target_patients <- readr::read_tsv(
    config$target_patients_file,
    locale = readr::locale(encoding = config$target_encoding),
    show_col_types = FALSE
  )
  clinical <- prepare_clinical_data(clinical_raw, target_patients, config$sex_column_policy)
  if (!nrow(clinical)) stop("No PDAC patients remain after v19 cohort filtering.", call. = FALSE)
  message("PDAC analysis cohort: ", nrow(clinical), " distinct patient/sample rows.")
  if (identical(config$sex_column_policy, "legacy") && !"Sex" %in% names(clinical)) {
    message("Compatibility note: Sex annotation is absent under v19 legacy policy. ",
            "Choose an explicit sex_column_policy to map a source column.")
  }
  message("Terminology note: v19's 'TMB' outputs count retained MAF rows per sample; ",
          "they are not normalized mutations/Mb.")

  maf_raw <- maftools::read.maf(maf = config$maf_file, clinicalData = clinical, verbose = FALSE)
  maf <- maftools::subsetMaf(maf = maf_raw, tsb = clinical$Tumor_Sample_Barcode)
  list(clinical = clinical, maf = maf, config = config)
}
