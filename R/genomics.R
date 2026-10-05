# Explicit common KRAS derivation, shared to prevent different subtype rules
# across figures. No analysis is run by sourcing this file.
normalize_kras <- function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(NA_character_)
    if (s %in% c("WT/Not detected", "No KRAS call", "Not detected", "WT")) return("No KRAS call")
    a <- trimws(unlist(strsplit(s, "[,;|]")))
    a <- sub("^p[.]", "", a)
    paste(sort(unique(a[nzchar(a)])), collapse = ";")
  }, character(1), USE.NAMES = FALSE)
}
kras_analysis_group <- function(x) {
  x <- normalize_kras(x)
  factor(ifelse(is.na(x), NA, ifelse(x %in% c("G12D", "G12V", "G12R", "No KRAS call"), x, "Other KRAS")),
    levels = c("No KRAS call", "G12D", "G12V", "G12R", "Other KRAS"))
}
derive_kras_maf <- function(variants) {
  k <- as.data.frame(variants[variants$Hugo_Symbol == "KRAS", ])
  if (!nrow(k)) return(data.frame(Tumor_Sample_Barcode = character(), KRAS_subtype_raw = character(),
    KRAS_MAF_group = character(), KRAS_MAF_exact_resolved = logical()))
  aa_cols <- intersect(c("HGVSp_Short", "Protein_Change", "Amino_Acid_Change", "AAChange", "HGVSp"), names(k))
  k$call <- "Other KRAS (unresolved)"
  k$resolved <- FALSE
  # Preserve the original validated G12/G13 coordinate rules. Do NOT infer an
  # exact Q61 subtype from the old ambiguous Q61H/Q61L/Q61R coordinate rule.
  rules <- c("25398284:C:T" = "G12D", "25398284:C:A" = "G12V", "25398285:C:G" = "G12R",
    "25398285:C:A" = "G12C", "25398285:C:T" = "G12S", "25398284:C:G" = "G12A", "25398281:C:T" = "G13D")
  valid_build <- as.character(k$NCBI_Build) %in% c("GRCh37", "37", "hg19")
  valid_chr <- as.character(k$Chromosome) %in% c("12", "chr12")
  key <- paste(k$Start_Position, k$Reference_Allele, k$Tumor_Seq_Allele2, sep = ":")
  ok <- valid_build & valid_chr & key %in% names(rules)
  k$call[ok] <- unname(rules[key[ok]]); k$resolved[ok] <- TRUE
  q61 <- valid_build & valid_chr & k$Start_Position %in% c(25380275, 25380276, 25380277) & !ok
  k$call[q61] <- "Q61 (unresolved)"
  # An explicit short protein annotation is preferred when available.
  if (length(aa_cols)) for (i in seq_len(nrow(k))) {
    a <- unique(unlist(lapply(aa_cols, function(v) {
      text <- as.character(k[[v]][i])
      if (is.na(text)) return(character())
      m <- regmatches(text, gregexpr("(?:p[.])?[A-Z][0-9]+[A-Z*]", text, perl = TRUE))[[1]]
      sub("^p[.]", "", m)
    })))
    if (length(a) == 1L) { k$call[i] <- a; k$resolved[i] <- TRUE }
  }
  k |>
    dplyr::group_by(Tumor_Sample_Barcode) |>
    dplyr::summarise(KRAS_subtype_raw = paste(sort(unique(call)), collapse = ";"),
      KRAS_MAF_exact_resolved = all(resolved), .groups = "drop") |>
    dplyr::mutate(KRAS_MAF_group = as.character(kras_analysis_group(KRAS_subtype_raw)))
}
add_source_genomics <- function(clinical, variants) {
  out <- dplyr::left_join(clinical, derive_kras_maf(variants), by = "Tumor_Sample_Barcode")
  no_call <- is.na(out$KRAS_subtype_raw)
  out$KRAS_subtype_raw[no_call] <- "No KRAS call"
  out$KRAS_MAF_exact_resolved[no_call] <- FALSE
  out$KRAS_MAF_group <- kras_analysis_group(out$KRAS_subtype_raw)
  out$KRAS_clinical_group <- kras_analysis_group(out$KRAS_subtype)
  out$KRAS_workbook_MAF_group <- kras_analysis_group(out$KRAS_subtype_MAF)
  counts <- dplyr::count(variants, Tumor_Sample_Barcode, name = "MAF_variant_count")
  dplyr::left_join(out, counts, by = "Tumor_Sample_Barcode")
}
