#!/usr/bin/env Rscript
# Run after the curated analysis. Cross-check outputs against untouched source
# cells and independently integrate a Kaplan-Meier curve for the RMST result.
source("R/config.R")
config <- load_config_from_command_line()
stopifnot(identical(config$clinical_schema, "curated_v3"))
read_table <- function(name) readr::read_tsv(output_path(config, name), show_col_types = FALSE)
raw <- as.data.frame(readxl::read_excel(config$clinical_file, sheet = "Main", na = c("", "NA")))
dup <- duplicated(raw$patient_id) | duplicated(raw$patient_id, fromLast = TRUE)
keep <- !dup
data <- read_table("18_Primary_analysis_dataset_private.tsv")
idx <- match(data$Tumor_Sample_Barcode, raw$Tumor_Sample_Barcode)
stopifnot(nrow(data) == sum(keep), setequal(data$Tumor_Sample_Barcode, raw$Tumor_Sample_Barcode[keep]),
  !anyDuplicated(data$patient_id), setequal(as.character(unique(data$M_stage)), c("0", "1")),
  identical(as.numeric(data$OS_m), as.numeric(raw$OS_months[idx])),
  identical(as.numeric(data$Age), as.numeric(raw$Age[idx])),
  all(data$survive == as.integer(raw$Death_event[idx] == "y")),
  isTRUE(all.equal(data$RFS_m, ifelse(raw$M_stage[idx] == 1, NA_real_, as.numeric(raw$DFS_months[idx])))),
  isTRUE(all.equal(data$Recur, as.numeric(raw$Recurrence_event[idx] == "y"))),
  all(data$driver_count == rowSums(data[, c("KRAS", "TP53", "SMAD4", "CDKN2A")])) )
legacy <- readxl::read_excel(output_path(config, "11_Driver_KRAS_patient_level_summary.xlsx"))
pos <- match(data$Tumor_Sample_Barcode, legacy$Tumor_Sample_Barcode)
stopifnot(nrow(legacy) == nrow(data), all(!is.na(pos)),
  all(as.character(data$KRAS_MAF) == as.character(legacy$KRAS_MAF_group[pos])),
  all(data$driver_count == legacy$driver_mutation_count[pos]))
risk <- read_table("21_OS_number_at_risk.tsv")
for (i in which(risk$variable == "KRAS")) {
  group <- if (risk$group[i] == "Retained variant") 1 else 0
  stopifnot(risk$n_risk[i] == sum(data$KRAS == group & data$OS_m >= risk$time[i]))
}
manual_rmst <- function(t, event, tau) {
  deaths <- sort(unique(t[event == 1 & t <= tau])); area <- 0; s <- 1; previous <- 0
  for (time in deaths) {
    area <- area + (time - previous) * s
    s <- s * (1 - sum(t == time & event == 1) / sum(t >= time))
    previous <- time
  }
  area + (tau - previous) * s
}
rmst <- read_table("21_Unadjusted_RMST_36months.tsv")
for (g in 0:1) {
  z <- data[data$KRAS == g, ]
  expected <- manual_rmst(z$OS_m, z$survive, 36)
  got <- rmst$rmst_months[rmst$variable == "KRAS" & rmst$group == as.character(g)]
  stopifnot(length(got) == 1, abs(got - expected) < 1e-8)
}
rec <- read_table("22_Observed_recurrence_patterns.tsv")
stopifnot(sum(rec$n) == sum(data$Recur == 1, na.rm = TRUE), all(rec$denominator_recurrent == sum(rec$n)))
co <- read_table("23_KRAS_report_MAF_concordance_private.tsv")
ct <- read_table("23_KRAS_report_MAF_concordance_counts.tsv")
stopifnot(nrow(co) == sum(keep), sum(ct$n) == sum(keep),
  all(co$comparable_group_agreement == (co$KRAS_report_group == co$KRAS_MAF)))
cox <- read_table("20_OS_Cox_coefficients.tsv")
stopifnot(all(cox$lower95 <= cox$HR), all(cox$HR <= cox$upper95), all(cox$p >= 0 & cox$p <= 1))
source_hash <- read_table("00_Input_checksums.tsv")
stopifnot(unname(tools::md5sum(config$clinical_file)) == source_hash$md5[source_hash$file == basename(config$clinical_file)])
integrated <- read_table("29_Integrated_analysis_dataset_private.tsv")
ii <- match(integrated$Tumor_Sample_Barcode, raw$Tumor_Sample_Barcode)
stopifnot(nrow(integrated) == sum(keep), !anyDuplicated(integrated$patient_id),
  !any(integrated$patient_duplicate), all(integrated$KRAS_subtype == raw$KRAS_subtype[ii]),
  isTRUE(all.equal(integrated$TMB, as.numeric(raw$TMB[ii]))),
  all(integrated$Stage_Group[integrated$M_stage == 1] == "IV"),
  !anyNA(integrated$MAF_variant_count), !anyNA(legacy$KRAS_MAF_group),
  all(as.character(legacy$KRAS_MAF_group[pos]) == integrated$KRAS_MAF_group[match(data$Tumor_Sample_Barcode, integrated$Tumor_Sample_Barcode)]))
maf <- maftools::read.maf(config$maf_file, verbose = FALSE)
cnt <- table(as.character(maf@data$Tumor_Sample_Barcode))
stopifnot(all(integrated$MAF_variant_count == as.numeric(cnt[integrated$Tumor_Sample_Barcode])))
pairs <- read_table("27_TMB_MAF_pairs_and_missingness_private.tsv")
corr <- read_table("27_TMB_MAF_rank_correlation.tsv")
ok <- complete.cases(pairs[, c("TMB", "MAF_variant_count")])
stopifnot(corr$n == sum(ok),
  abs(corr$spearman_rho - cor(pairs$TMB[ok], pairs$MAF_variant_count[ok], method = "spearman")) < 1e-12,
  corr$lower95 <= corr$spearman_rho, corr$upper95 >= corr$spearman_rho)
ka <- read_table("25_KRAS_source_comparison_private.tsv")
stopifnot(nrow(ka) == nrow(integrated),
  all(ka$KRAS_group_agreement == (ka$KRAS_clinical_group == ka$KRAS_MAF_group)),
  all(is.na(ka$KRAS_exact_agreement[!ka$KRAS_exact_comparable])))
mapping <- read_table("28_Applied_variable_mapping.tsv")
stopifnot(all(mapping$type[mapping$variable %in% c("Age", "Tumor_size")] == "continuous"),
  all(mapping$type[mapping$variable %in% c("ASA", "T_stage", "N_stage", "LN_positive", "M_stage")] == "categorical"),
  setequal(mapping$display_label[mapping$variable %in% c("Preop_platinum_exposure", "Operation",
    "Adjuvant_treatment", "Adjuvant_chemotherapy", "Adjuvant_radiotherapy")],
    c("Preop_platinum", "OP", "Adjuvant_Tx", "Adjuvant_CT", "Adjuvant_RT")))
sr <- read_table("26_KRAS_KM_number_at_risk.tsv")
for (i in seq_len(nrow(sr))) {
  v <- sr$variable[i]; time <- if (sr$endpoint[i] == "OS") integrated$OS_months else integrated$RFS_m
  event <- if (sr$endpoint[i] == "OS") integrated$survive else integrated$Recur
  stopifnot(sr$n_risk[i] == sum(integrated[[v]] == sr$group[i] & time >= sr$time[i] & !is.na(event), na.rm = TRUE))
}
source_cox <- read_table("29_Source_specific_Cox_coefficients.tsv")
common <- subset(source_cox, cohort == "Same paired patients")
for (e in unique(common$endpoint)) stopifnot(length(unique(common$n[common$endpoint == e])) == 1)
for (v in c("KRAS_clinical_group", "KRAS_MAF_group")) for (e in c("OS", "Recorded_recurrence")) {
  stopifnot(file.exists(output_path(config, paste0("26_", v, "_", e, ".png"))))
}
cat("Curated output validation passed: source values/cohort, legacy-vs-new genomics, risk tables, independent KM integration, recurrence totals, KRAS concordance, Cox bounds, unchanged source hash.\n")
