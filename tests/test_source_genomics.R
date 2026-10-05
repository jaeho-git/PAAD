#!/usr/bin/env Rscript
# Synthetic-only tests for exact-vs-grouped provenance and multiple calls.
source("R/genomics.R")
source("R/plot_helpers.R")
annotation <- data.frame(Tumor_size = c(1.5, 3), Age = c(55, 70), ASA = factor(c(1, 2)))
colors <- make_custom_colors(annotation, names(annotation))
stopifnot(is.function(colors$Tumor_size), is.function(colors$Age), is.character(colors$ASA))
stopifnot(identical(normalize_kras(c("G12V, G12D", "G12D;G12V", "WT/Not detected", "Not detected", NA)),
  c("G12D;G12V", "G12D;G12V", "No KRAS call", "No KRAS call", NA_character_)))
stopifnot(as.character(kras_analysis_group(c("Other", "Not detected"))) == c("Other KRAS", "No KRAS call"))
v <- data.frame(Tumor_Sample_Barcode = c("S1", "S2", "S3", "S4", "S4", "S5"),
  Hugo_Symbol = "KRAS", NCBI_Build = c(rep("GRCh37", 5), "GRCh38"), Chromosome = "12",
  Start_Position = c(25398284, 25398285, 25380275, 25398284, 25398285, 25398284),
  Reference_Allele = c("C", "C", "T", "C", "C", "C"), Tumor_Seq_Allele2 = c("T", "G", "G", "T", "G", "T"))
k <- derive_kras_maf(v)
stopifnot(k$KRAS_subtype_raw[k$Tumor_Sample_Barcode == "S1"] == "G12D",
  k$KRAS_subtype_raw[k$Tumor_Sample_Barcode == "S2"] == "G12R",
  k$KRAS_subtype_raw[k$Tumor_Sample_Barcode == "S3"] == "Q61 (unresolved)",
  !k$KRAS_MAF_exact_resolved[k$Tumor_Sample_Barcode == "S3"],
  k$KRAS_subtype_raw[k$Tumor_Sample_Barcode == "S4"] == "G12D;G12R",
  k$KRAS_MAF_group[k$Tumor_Sample_Barcode == "S4"] == "Other KRAS",
  !k$KRAS_MAF_exact_resolved[k$Tumor_Sample_Barcode == "S5"],
  is.character(k$KRAS_MAF_group),
  ifelse(TRUE, k$KRAS_MAF_group[1], "No KRAS call") == "G12D")
v$HGVSp_Short <- c(NA, NA, "p.Q61H", NA, NA, "p.G12V")
k <- derive_kras_maf(v)
stopifnot(k$KRAS_subtype_raw[k$Tumor_Sample_Barcode == "S3"] == "Q61H",
  k$KRAS_MAF_exact_resolved[k$Tumor_Sample_Barcode == "S3"],
  k$KRAS_subtype_raw[k$Tumor_Sample_Barcode == "S5"] == "G12V")
clinical <- data.frame(Tumor_Sample_Barcode = c("S1", "S6"),
  KRAS_subtype = c("G12V", "WT/Not detected"),
  KRAS_subtype_MAF = c("G12D", "Not detected"), TMB = c(4, NA))
a <- add_source_genomics(clinical, v)
stopifnot(identical(a$KRAS_subtype, clinical$KRAS_subtype), identical(a$TMB, clinical$TMB),
  as.character(a$KRAS_workbook_MAF_group) == c("G12D", "No KRAS call"),
  a$KRAS_subtype_raw[2] == "No KRAS call", !a$KRAS_MAF_exact_resolved[2],
  is.na(a$MAF_variant_count[2]))
cat("Synthetic KRAS provenance/normalization tests passed.\n")
