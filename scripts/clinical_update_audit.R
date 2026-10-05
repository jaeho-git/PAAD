#!/usr/bin/env Rscript
# Compare the curated workbook with the previous inputs, without editing either.
# All row-level reports are private outputs. Run from the repository root.
source("R/config.R")
source("R/data.R")
config <- load_config_from_command_line()
if (!identical(config$clinical_schema, "curated_v3")) {
  message("Clinical update audit skipped: legacy schema.")
  quit(status = 0)
}
load_analysis_packages()
new_raw <- readxl::read_excel(config$clinical_file, sheet = config$clinical_sheet)
dictionary <- readxl::read_excel(config$clinical_file, sheet = "Data_dictionary")
new <- prepare_curated_clinical(new_raw, dictionary)
old_file <- file.path(config$project_root, "data/raw/PDAC_clinical_info.xlsx")
if (!is.null(config$previous_clinical_file)) old_file <- normalize_config_path(config$previous_clinical_file, config$project_root, TRUE)
if (!file.exists(old_file)) stop("Previous clinical workbook required for the comparison audit: ", old_file)
old <- readxl::read_excel(old_file)
target <- readr::read_tsv(config$target_patients_file, locale = readr::locale(encoding = config$target_encoding), show_col_types = FALSE)
maf <- readr::read_tsv(config$maf_file, comment = "#", show_col_types = FALSE)
key <- "Tumor_Sample_Barcode"
write_report <- function(x, name) readr::write_tsv(as.data.frame(x), output_path(config, name), na = "NA")

linkage <- new[, c(key, "patient_id", "NGS_report_id", "patient_duplicate", "M_stage", "NGS_group", "primary_eligible")]
linkage$old_clinical_match <- new[[key]] %in% old[[key]]
linkage$old_target_match <- new[[key]] %in% target[[key]]
linkage$maf_match <- new[[key]] %in% maf[[key]]
linkage$old_classification <- old$Classification[match(new[[key]], old[[key]])]
write_report(linkage, "00_Barcode_linkage_private.tsv")
write_report(new_raw[new$patient_duplicate, ], "00_Duplicate_patient_records_private.tsv")
old_pdac <- old[old[[key]] %in% target[[key]] & !is.na(old$Classification) & grepl("PDAC", old$Classification, ignore.case = TRUE), ]
write_report(old_pdac[!old_pdac[[key]] %in% new[[key]], c(key, "Classification")], "00_Previous_cohort_not_in_update_private.tsv")

mapping <- data.frame(
  new = c("Age", "Sex", "BMI", "CA19_9", "CEA", "Neoadjuvant", "Differentiation", "Tumor_size", "T_stage", "N_stage", "LN_positive", "M_stage", "AJCC_stage", "R_status", "LVI", "PNI", "OS_months", "DFS_months", "Death_event", "Recurrence_event"),
  old = c("age", "성별코드", "BMI", "CA 19-9", "CEA", "NAC", "Differentiation", "Size", "T", "N", "N status", "M", "Stage", "RM", "LVI", "PNI", "OS_m", "RFS_m", "survive", "Recur"),
  source = c("target", rep("clinical", 12), rep("target", 7)),
  meaning = c("Pathology age vs previous target age; definitions can differ", "F/M", "kg/m^2", "Diagnosis value; assay unit unconfirmed", "Diagnosis value; assay unit unconfirmed", "y/n", "wd/md/pd", "Source units unconfirmed", "Curated T", "Curated N", "Node-positive indicator", "M category at diagnosis", "Curated stage; AJCC edition unconfirmed", "R0/R1 vs negative/positive", "y/n vs 0/1", "y/n vs 0/1", "New surgery-origin OS; previous time definition unconfirmed", "Recurrence-only follow-up; NOT death-inclusive DFS", "y=death vs legacy event=1", "y=recurrence vs legacy event=1")
)
normalize <- function(x, field) {
  x <- tolower(trimws(as.character(x))); x[x %in% c("", "na", "n/a", "#null!")] <- NA_character_
  if (field %in% c("LVI", "PNI", "Death_event", "Recurrence_event", "Neoadjuvant")) {
    x[x %in% c("y", "yes", "positive", "1", "1.0")] <- "1"
    x[x %in% c("n", "no", "negative", "0", "0.0")] <- "0"
  }
  if (field == "R_status") {
    x[x %in% c("r0", "negative", "0", "0.0")] <- "0"
    x[x %in% c("r1", "positive", "1", "1.0")] <- "1"
  }
  if (field %in% c("Age", "BMI", "CA19_9", "CEA", "Tumor_size", "T_stage", "N_stage", "LN_positive", "M_stage", "OS_months", "DFS_months")) {
    x <- suppressWarnings(as.numeric(gsub("[<>=[:space:]]", "", x)))
  }
  x
}
summaries <- details <- list()
for (i in seq_len(nrow(mapping))) {
  map <- mapping[i, ]; previous <- if (map$source == "target") target else old
  pos <- match(new_raw[[key]], previous[[key]])
  a <- normalize(new_raw[[map$new]], map$new); b <- normalize(previous[[map$old]][pos], map$new)
  matched <- !is.na(pos); comparable <- matched & !is.na(a) & !is.na(b)
  equal <- if (is.numeric(a)) abs(a - b) < 1e-6 else a == b
  changed <- comparable & !is.na(equal) & !equal
  summaries[[i]] <- cbind(map, matched = sum(matched), both_present = sum(comparable),
    changed_values = sum(changed), new_available_old_missing = sum(matched & !is.na(a) & is.na(b)),
    new_missing_old_available = sum(matched & is.na(a) & !is.na(b)))
  indices <- which(changed | (matched & xor(is.na(a), is.na(b))))
  details[[i]] <- data.frame(Tumor_Sample_Barcode = new_raw[[key]][indices], variable = rep(map$new, length(indices)),
    old_value = as.character(previous[[map$old]][pos[indices]]), new_value = as.character(new_raw[[map$new]][indices]))
}
write_report(dplyr::bind_rows(summaries), "00_Column_mapping_and_changes.tsv")
write_report(dplyr::bind_rows(details), "00_Changed_values_private.tsv")
old_pos <- match(new_raw[[key]], old[[key]])
target_pos <- match(new_raw[[key]], target[[key]])
equal_text <- function(x, y) sum(!is.na(x) & !is.na(y) & trimws(as.character(x)) == trimws(as.character(y)))
identity_checks <- data.frame(check = c("Barcode to old NGS report ID", "Barcode to old clinical age (검사나이)", "Barcode to target NGS report ID", "Barcode to target patient ID"),
  comparable_n = c(nrow(new_raw), nrow(new_raw), sum(!is.na(target_pos)), sum(!is.na(target_pos))),
  equal_n = c(equal_text(new_raw$NGS_report_id, old$임상리포트ID_NGS[old_pos]),
    equal_text(new_raw$Age, old$검사나이[old_pos]), equal_text(new_raw$NGS_report_id, target$ID_NGS[target_pos]),
    equal_text(new_raw$patient_id, target$ID[target_pos])))
write_report(identity_checks, "00_Identifier_and_age_crosschecks.tsv")
missing <- data.frame(variable = names(new_raw), n = nrow(new_raw), missing_n = vapply(new_raw, function(x) sum(is.na(x) | trimws(as.character(x)) %in% c("", "NA", "N/A")), integer(1)))
missing$missing_percent <- missing$missing_n / missing$n * 100
write_report(missing, "00_Missing_data.tsv")

counts <- data.frame(metric = c("Updated sample rows", "Unique patient IDs", "Duplicate-patient rows excluded", "M1 sample rows", "Primary M0+M1 unique-patient cohort", "Old clinical matches", "MAF matches", "Old target-list matches", "Previous PDAC cohort", "Previous PDAC absent from update", "New samples not in old target list", "Deaths without recorded recurrence", "Death later than last recurrence follow-up"),
 n = c(nrow(new), dplyr::n_distinct(new$patient_id), sum(new$patient_duplicate), sum(new$M_stage == "1"), sum(new$primary_eligible), sum(linkage$old_clinical_match), sum(linkage$maf_match), sum(linkage$old_target_match), nrow(old_pdac), sum(!old_pdac[[key]] %in% new[[key]]), sum(!linkage$old_target_match), sum(new$Recur == 0 & new$survive == 1, na.rm = TRUE), sum(new$Recur == 0 & new$survive == 1 & new$OS_m > new$RFS_m + 0.05, na.rm = TRUE)))
write_report(counts, "00_Cohort_audit_counts.tsv")
issues <- data.frame(item = c("Repeated patient identifiers", "Recurrence endpoint", "Competing risks", "Tumor size and assay units", "Stage/margin definitions", "NGS timing", "Variant interpretation"),
  handling = c("Retain all source rows for descriptive audit; exclude all rows from duplicated patient IDs in primary inference until representative specimens are selected.",
  "Use recorded recurrence event and follow-up; label recurrence-only. Do not call this death-inclusive DFS/RFS.",
  "Not estimated: death can occur after last recurrence assessment; recurrence-free status at death and first-event ascertainment require confirmation.",
  "No conversion and no cm/U-per-mL label until source units are confirmed.",
  "Use curated T/N/stage and R0/R1; do not infer AJCC edition or margin-distance cutoff.",
  "NGS_group definitions/timing unconfirmed. No interpretation as anatomic resectability or delayed-entry times.",
  "MAF has no pathogenicity, CNA, LOH or germline annotations. Absent retained calls are not proof of biological wild type."))
write_report(issues, "00_Definition_issues.tsv")
files <- c(config$clinical_file, old_file, config$target_patients_file, config$maf_file)
write_report(data.frame(file = basename(files), md5 = unname(tools::md5sum(files))), "00_Input_checksums.tsv")
message("Clinical update audit complete. Reports: ", config$output_dir)
