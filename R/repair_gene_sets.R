# McIntyre et al., Cancer 2020, Methods p3940, DOI 10.1002/cncr.33038.
# A variant in a listed gene does NOT establish functional HRD or dMMR.
repair_gene_sets <- function() list(
  HRD = c("ARID1A", "ATM", "BAP1", "BARD1", "BLM", "BRCA1", "BRCA2",
    "BRIP1", "CHEK2", "FAM175A", "FANCA", "FANCC", "NBN", "PALB2",
    "RAD50", "RAD51", "RAD51C", "RTEL1"),
  MMR = c("MLH1", "MSH2", "MSH6", "PMS2"))
# The paper uses the former symbol FAM175A; ABRAXAS1 is its current synonym.
# RAD51B/D are NOT substituted for RAD51 or RAD51C.
repair_paper_symbol <- function(x) {
  x <- as.character(x)
  x[x == "ABRAXAS1"] <- "FAM175A"
  x
}
repair_group_status <- function(clinical, variants) {
  sets <- repair_gene_sets()
  variants$paper_gene <- repair_paper_symbol(variants$Hugo_Symbol)
  d <- clinical
  for (set in names(sets)) {
    m <- matrix(FALSE, nrow(d), length(sets[[set]]), dimnames = list(d$Tumor_Sample_Barcode, sets[[set]]))
    for (g in sets[[set]]) m[, g] <- rownames(m) %in%
      variants$Tumor_Sample_Barcode[variants$paper_gene == g]
    d[[paste0(set, "_gene_count")]] <- rowSums(m)
    d[[paste0(set, "_variant_status")]] <- factor(ifelse(rowSums(m) > 0,
      "Variant detected", "No retained variant"), levels = c("No retained variant", "Variant detected"))
  }
  d$Repair_pattern <- factor(ifelse(d$HRD_gene_count > 0 & d$MMR_gene_count > 0, "Both gene lists",
    ifelse(d$HRD_gene_count > 0, "HRD-list only", ifelse(d$MMR_gene_count > 0, "MMR-list only", "Neither gene list"))),
    levels = c("Neither gene list", "HRD-list only", "MMR-list only", "Both gene lists"))
  d
}
