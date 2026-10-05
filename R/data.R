# Shared input preparation only.

# Display-only labels: raw column names and original data are preserved.
# Approved workbook mapping: original column names in data, short display labels.
clinical_label <- function(x) {
  labels <- c(Preop_platinum_exposure = "Preop_platinum", Operation = "OP",
    Adjuvant_treatment = "Adjuvant_Tx", Adjuvant_chemotherapy = "Adjuvant_CT",
    Adjuvant_radiotherapy = "Adjuvant_RT")
  out <- x; i <- x %in% names(labels); out[i] <- unname(labels[x[i]]); out
}
# Clear figure labels are separate from approved analysis-column aliases.
figure_label <- function(x) {
  labels <- c(Age = "Age (years)", Age10 = "Age (per 10 years)", Age_group = "Age group",
    Tumor_size = "Tumor size", TMB = "Reported TMB", MAF_variant_count = "Variant count",
    Preop_platinum_exposure = "Preoperative platinum", Operation = "Surgery type",
    Adjuvant_treatment = "Adjuvant treatment", Adjuvant_chemotherapy = "Adjuvant chemotherapy",
    Adjuvant_radiotherapy = "Adjuvant radiotherapy", Neoadjuvant = "Neoadjuvant treatment",
    T_stage = "T category", N_stage = "N category", LN_positive = "Lymph-node involvement",
    M_stage = "M category", AJCC_stage = "AJCC stage", Stage_Group = "Stage group",
    R_status = "Resection margin", LVI = "Lymphovascular invasion", PNI = "Perineural invasion",
    NGS_group = "NGS group (source category)", MSI = "MSI (clinical report)",
    OS_months = "Overall survival (months)", DFS_months = "Recorded recurrence follow-up (months)",
    KRAS_subtype = "Clinical KRAS subtype", KRAS_subtype_MAF = "Workbook MAF-based KRAS subtype",
    KRAS_subtype_raw = "Re-derived MAF KRAS subtype",
    KRAS_clinical_group = "Clinical KRAS subtype", KRAS_MAF_group = "Re-derived MAF KRAS subtype",
    KRAS_workbook_MAF_group = "Workbook MAF-based KRAS subtype")
  out <- gsub("_", " ", x); i <- x %in% names(labels); out[i] <- unname(labels[x[i]]); out
}
excluded_analysis_fields <- function(x) {
  x[grepl("purity|cellularity|cellarity|(^|_)loh($|_)", x, ignore.case = TRUE)]
}
curated_continuous <- function() c("Age", "BMI", "CA19_9", "CEA", "Tumor_size", "TMB", "OS_months", "DFS_months")
curated_categorical <- function() c("Sex", "Diabetes", "Other_cancer", "NGS_group", "Neoadjuvant",
  "Preop_platinum_exposure", "Operation", "ASA", "Differentiation", "LVI", "PNI",
  "T_stage", "N_stage", "LN_positive", "M_stage", "AJCC_stage", "R_status",
  "HG_PanIN", "IPMN_HGD", "Adjuvant_treatment", "Adjuvant_chemotherapy",
  "Adjuvant_radiotherapy", "KRAS_subtype", "MSI", "Death_event", "Recurrence_event",
  "Recurrence_pattern", "Distant_pattern")
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

prepare_curated_clinical <- function(raw, dictionary) {
  # The updated workbook is authoritative. Never replace its outcomes or
  # pathology with the old target-patient TXT. Names/meaning were checked against
  # Data_dictionary; all source fields remain available for the new analyses.
  required <- c("Tumor_Sample_Barcode", "patient_id", "Age", "Sex", "BMI",
    "CA19_9", "CEA", "Neoadjuvant", "Differentiation", "Tumor_size",
    "T_stage", "N_stage", "LN_positive", "M_stage", "AJCC_stage", "R_status",
    "LVI", "PNI", "KRAS_subtype", "KRAS_subtype_MAF", "TMB", "Recurrence_event", "DFS_months",
    "Death_event", "OS_months")
  missing <- setdiff(required, names(raw))
  if (length(missing)) stop("Curated clinical columns missing: ", paste(missing, collapse = ", "))
  if (!all(c("Analysis_variable", "Definition/Coding") %in% names(dictionary)) ||
      !all(required %in% dictionary$Analysis_variable)) {
    stop("Data_dictionary does not document all required analysis columns.")
  }
  clean <- function(x) {
    y <- trimws(as.character(x)); y[y %in% c("", "NA", "N/A")] <- NA_character_; y
  }
  binary <- function(x, variable) {
    y <- clean(x)
    if (any(!is.na(y) & !y %in% c("y", "n"))) stop("Unexpected y/n code in ", variable)
    ifelse(is.na(y), NA_integer_, as.integer(y == "y"))
  }
  raw <- as.data.frame(raw)
  raw[] <- lapply(raw, function(x) if (is.character(x)) clean(x) else x)
  raw$Tumor_Sample_Barcode <- clean(raw$Tumor_Sample_Barcode)
  if (anyNA(raw$Tumor_Sample_Barcode) || anyDuplicated(raw$Tumor_Sample_Barcode)) {
    stop("Curated barcode must be nonmissing and unique; resolve duplicates before analysis.")
  }
  raw$patient_id <- clean(raw$patient_id)
  if (anyNA(raw$patient_id)) stop("Curated patient_id is missing.")
  raw$patient_duplicate <- duplicated(raw$patient_id) | duplicated(raw$patient_id, fromLast = TRUE)
  for (v in c("Age", "BMI", "CA19_9", "CEA", "Tumor_size", "TMB", "DFS_months", "OS_months")) {
    parsed <- suppressWarnings(as.numeric(raw[[v]]))
    if (any(!is.na(raw[[v]]) & (is.na(parsed) | !is.finite(parsed)))) stop("Non-numeric value in ", v)
    raw[[v]] <- parsed
  }
  for (v in c("T_stage", "N_stage", "LN_positive", "M_stage")) raw[[v]] <- as.character(raw[[v]])
  allowed <- list(Sex = c("F", "M"), T_stage = c("1", "2", "3", "4", "0"),
    N_stage = c("0", "1", "2"), LN_positive = c("0", "1"), M_stage = c("0", "1"),
    R_status = c("R0", "R1", "R2"), Differentiation = c("wd", "md", "pd"))
  for (v in names(allowed)) if (any(!is.na(raw[[v]]) & !raw[[v]] %in% allowed[[v]])) stop("Unrecognized curated category: ", v)
  raw$Classification <- "PDAC" # Curated PDAC cohort; old workbook cross-check is exported by the audit.
  raw$NAC <- ifelse(binary(raw$Neoadjuvant, "Neoadjuvant") == 1L, "Yes", "No")
  raw$Size <- raw$Tumor_size # Source units are unconfirmed: no unit conversion.
  raw$T <- ifelse(is.na(raw$T_stage), NA_character_, paste0("T", raw$T_stage))
  raw$N <- ifelse(is.na(raw$N_stage), NA_character_, paste0("N", raw$N_stage))
  raw$N_status <- ifelse(is.na(raw$LN_positive), NA_character_, paste0("N", raw$LN_positive))
  raw$Stage <- raw$AJCC_stage
  raw$Stage_Group <- sub("[AB]$", "", raw$Stage)
  raw$Differentiation <- toupper(raw$Differentiation)
  for (v in c("LVI", "PNI")) raw[[v]] <- ifelse(binary(raw[[v]], v) == 1L, "Positive", "Negative")
  raw$RM <- ifelse(raw$R_status == "R0", "Negative", ifelse(raw$R_status %in% c("R1", "R2"), "Positive", NA_character_))
  raw$OS_m <- raw$OS_months
  raw$survive <- binary(raw$Death_event, "Death_event")
  raw$RFS_m <- raw$DFS_months # Legacy alias only. Plot labels explicitly say recurrence-only.
  raw$Recur <- binary(raw$Recurrence_event, "Recurrence_event")
  if (any(raw$OS_m <= 0, na.rm = TRUE) || any(raw$RFS_m <= 0, na.rm = TRUE)) {
    stop("Nonpositive survival/recurrence times require review before analysis.")
  }
  if (any(raw$Recur == 1 & raw$RFS_m > raw$OS_m + 0.05, na.rm = TRUE)) {
    stop("A recorded recurrence occurs after the OS follow-up/death time; review source times.")
  }
  raw$RFS_m[raw$M_stage == "1"] <- NA_real_
  raw$Recur[raw$M_stage == "1"] <- NA_integer_
  raw$KRAS_report <- raw$KRAS_subtype
  # Both workbook KRAS fields are preserved. KRAS_subtype_MAF is the selected
  # manuscript source when kras_source = "workbook_maf"; KRAS_subtype remains
  # available for private audit only and is never mixed into that analysis.
  raw$TMB_report <- raw$TMB
  # MAF_variant_count is computed separately in the genomic analyses.
  raw$Age_group <- ifelse(raw$Age <= 55, "Age <=55", ifelse(raw$Age >= 70, "Age >=70", "Age 56-69"))
  for (v in intersect(c("ASA", "T_stage", "N_stage", "LN_positive", "M_stage"), names(raw))) raw[[v]] <- factor(raw[[v]])
  raw$primary_eligible <- !raw$patient_duplicate
  raw
}

load_analysis_data <- function(config) {
  # Read and connect the clinical workbook, target-patient table, and MAF.
  # The returned list is the complete common input used by scripts/*.R.
  load_analysis_packages()
  clinical_raw <- readxl::read_excel(config$clinical_file, sheet = config$clinical_sheet)
  if (identical(config$clinical_schema, "curated_v3")) {
    dictionary <- readxl::read_excel(config$clinical_file, sheet = "Data_dictionary")
    clinical_all <- prepare_curated_clinical(clinical_raw, dictionary)
    keep <- if (config$updated_cohort == "localized_unique") (clinical_all$primary_eligible & !is.na(clinical_all$M_stage) & clinical_all$M_stage == "0") else clinical_all$primary_eligible
    clinical <- clinical_all[keep, , drop = FALSE]
    message("Curated workbook: ", nrow(clinical_all), " samples; analysis cohort: ", nrow(clinical),
            "; policy: ", config$updated_cohort, ". Old target list is NOT used to select or overwrite patients.")
    # Establish MAF coverage before any absent mutation is coded as zero.
    maf_table <- readr::read_tsv(config$maf_file, comment = "#", show_col_types = FALSE)
    missing_maf <- setdiff(clinical_all$Tumor_Sample_Barcode, maf_table$Tumor_Sample_Barcode)
    if (length(missing_maf)) stop(length(missing_maf), " curated barcodes have no MAF records; do not infer wild type.")
    maf_raw <- maftools::read.maf(maf = config$maf_file, clinicalData = clinical, verbose = FALSE)
    maf <- maftools::subsetMaf(maf = maf_raw, tsb = clinical$Tumor_Sample_Barcode)
    message("MAF calls are retained variants, not confirmed pathogenic alterations; no CNA/LOH or germline data.")
    return(list(clinical = clinical, clinical_all = clinical_all, dictionary = dictionary, maf = maf, config = config))
  }
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
