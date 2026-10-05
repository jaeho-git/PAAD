#!/usr/bin/env Rscript
# Independent patient-count, source-list, denominator and exclusion checks.
source("R/config.R")
source("R/data.R")
source("R/repair_gene_sets.R")
config <- load_config_from_command_line()
read <- function(f) readr::read_tsv(output_path(config, f), show_col_types = FALSE)
sets <- repair_gene_sets()
stopifnot(length(sets$HRD) == 18L, length(sets$MMR) == 4L,
  !length(intersect(sets$HRD, sets$MMR)),
  identical(repair_paper_symbol(c("ABRAXAS1", "RAD51")), c("FAM175A", "RAD51")),
  !"Tumor_cellularity" %in% curated_continuous())
# A repeated variant must count as one patient and one altered gene.
toy <- data.frame(Tumor_Sample_Barcode = c("S1", "S2", "S3"))
calls <- data.frame(Tumor_Sample_Barcode = c("S1", "S1", "S2", "S2"),
  Hugo_Symbol = c("ABRAXAS1", "ABRAXAS1", "ATM", "MLH1"))
check <- repair_group_status(toy, calls)
stopifnot(identical(check$HRD_gene_count, c(1, 1, 0)),
  identical(check$MMR_gene_count, c(0, 1, 0)),
  identical(as.character(check$Repair_pattern), c("HRD-list only", "Both gene lists", "Neither gene list")))
input <- load_analysis_data(config)
maf <- as.data.frame(input$maf@data)
maf$Hugo_Symbol <- repair_paper_symbol(maf$Hugo_Symbol)
actual <- read("31_Repair_patient_status_private.tsv")
stopifnot(nrow(actual) == nrow(input$clinical), !anyDuplicated(actual$patient_id),
  setequal(actual$Tumor_Sample_Barcode, input$clinical$Tumor_Sample_Barcode))
freq <- read("31_Repair_gene_patient_frequencies.tsv")
for (i in seq_len(nrow(freq))) {
  v <- maf[maf$Hugo_Symbol == freq$gene[i], ]
  stopifnot(freq$patients[i] == length(unique(v$Tumor_Sample_Barcode)),
    freq$variant_rows[i] == nrow(v), freq$total_patients[i] == nrow(actual),
    abs(freq$percent[i] - 100 * freq$patients[i] / nrow(actual)) < 1e-10)
}
summ <- read("31_Repair_gene_list_summary.tsv")
for (s in names(sets)) {
  ids <- unique(as.character(maf$Tumor_Sample_Barcode[maf$Hugo_Symbol %in% sets[[s]]]))
  stopifnot(all((actual[[paste0(s, "_gene_count")]] > 0) == (actual$Tumor_Sample_Barcode %in% ids)),
    summ$variant_positive_n[summ$gene_list == s] == length(ids))
  for (suffix in c("all_patients", "variant_positive_patients"))
    stopifnot(file.exists(output_path(config, paste0("31_", s, "_oncoplot_", suffix, ".png"))))
}
overlap <- read("31_Repair_gene_list_overlap.tsv")
stopifnot(sum(overlap$n) == nrow(actual), abs(sum(overlap$percent) - 100) < 1e-10)
desc <- read("32_Repair_group_clinical_distributions.tsv")
stopifnot(!length(excluded_analysis_fields(desc$variable)))
for (z in split(desc, paste(desc$grouping, desc$variable, desc$group))) {
  stopifnot(all(z$total_n == z$observed_n + z$missing_n))
  if (z$category[1] != "Continuous") stopifnot(sum(z$n) == z$total_n[1], abs(sum(z$percent) - 100) < 1e-9)
}
tt <- read("32_Repair_group_clinical_tests.tsv")
for (z in split(tt, tt$grouping)) {
  expected <- p.adjust(z$p, "BH")
  stopifnot(isTRUE(all.equal(z$q_BH, expected)))
}
mapping <- read("28_Applied_variable_mapping.tsv")
stopifnot(!length(excluded_analysis_fields(mapping$variable)), "figure_label" %in% names(mapping))
cat("Repair outputs validated: lists, alias, unique patient counts, all denominators, BH, exclusion rules and oncoplots.\n")
