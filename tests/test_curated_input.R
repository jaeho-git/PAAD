#!/usr/bin/env Rscript
# Isolated schema/endpoint regression tests using wholly synthetic in-memory data.
source("R/data.R")
n <- 6L
raw <- data.frame(Tumor_Sample_Barcode = paste0("SYN-V3-", seq_len(n)),
  patient_id = c("SYN-P1", "SYN-P1", "SYN-P3", "SYN-P4", "SYN-P5", "SYN-P6"),
  Age = 50:55, Sex = rep(c("F", "M"), 3), BMI = 22, CA19_9 = NA_real_, CEA = 1,
  Neoadjuvant = rep(c("n", "y"), 3), Differentiation = c("wd", "md", "pd", "NA", "wd", "md"),
  Tumor_size = 2.5, T_stage = 2, N_stage = c(0, 0, 1, 2, 0, 1), LN_positive = c(0, 0, 1, 1, 0, 1),
  M_stage = c(0, 0, 1, 0, 0, 0), AJCC_stage = c("IB", "IB", "IV", "III", "IB", "IIB"),
  R_status = c("R0", "R1", "R0", "R0", "R0", "R1"), LVI = "y", PNI = c("y", "n", "NA", "n", "y", "y"),
  KRAS_subtype = c("G12D", "G12D", "WT/Not detected", "G12V", "G12R", "Q61H"),
  KRAS_subtype_MAF = c("G12D", "G12D", "Not detected", "G12V", "G12R", "Other"),
  TMB = c(2, NA, 5, 8, 7, 4),
  Recurrence_event = c("n", "n", "NA", "y", "n", "y"), DFS_months = c(8, 8, 9, 7, 20, 6),
  Death_event = c("n", "n", "y", "y", "n", "y"), OS_months = c(12, 12, 10, 13, 22, 10))
dict <- data.frame(Analysis_variable = names(raw), `Definition/Coding` = "Synthetic test definition", check.names = FALSE)
d <- prepare_curated_clinical(raw, dict)
stopifnot(sum(d$patient_duplicate) == 2L, sum(d$primary_eligible) == 4L,
  identical(d$Age, as.numeric(raw$Age)), identical(d$OS_m, as.numeric(raw$OS_months)),
  d$survive[4] == 1L, d$Recur[4] == 1L, is.na(d$Recur[3]), is.na(d$RFS_m[3]),
  is.na(d$PNI[3]), is.na(d$Differentiation[4]), d$RM[1] == "Negative", d$RM[2] == "Positive",
  d$TMB_report[4] == 8, "TMB" %in% names(d), "KRAS_subtype" %in% names(d),
  "KRAS_subtype_MAF" %in% names(d), d$KRAS_subtype_MAF[6] == "Other",
  is.factor(d$T_stage), is.factor(d$M_stage), d$KRAS_report[3] == "WT/Not detected")
must_fail <- function(expr) stopifnot(inherits(tryCatch({ force(expr); NULL }, error = identity), "error"))
bad <- raw; bad$Tumor_Sample_Barcode[2] <- bad$Tumor_Sample_Barcode[1]
must_fail(prepare_curated_clinical(bad, dict))
bad <- raw; bad$Death_event[2] <- "unknown-code"
must_fail(prepare_curated_clinical(bad, dict))
bad <- raw; bad$OS_months[2] <- "not numeric"
must_fail(prepare_curated_clinical(bad, dict))
must_fail(prepare_curated_clinical(raw, dict[-1, ]))
example <- read.delim("data/example/curated_v3_main.example.tsv", check.names = FALSE, na.strings = c("", "NA"))
example_dict <- read.delim("data/example/curated_v3_dictionary.example.tsv", check.names = FALSE)
stopifnot(ncol(example) == 41L, nrow(example_dict) == 41L,
  identical(names(example), example_dict$Analysis_variable), all(grepl("^SYN-", example$patient_id)),
  sum(prepare_curated_clinical(example, example_dict)$primary_eligible) == sum(!(duplicated(example$patient_id) | duplicated(example$patient_id, fromLast = TRUE))))
cat("Curated schema/endpoint regression tests passed.\n")
