# Count variant occurrences, not occupied oncoplot cells or alteration classes.
# One event in two patients counts twice. Repeated annotation of the same event
# in the same patient counts once. A multi-gene annotation can contribute once
# to each gene, but only once to that patient's total.
oncoplot_variant_counts <- function(variants, sample_ids, genes) {
  v <- as.data.frame(variants)
  keys <- c("Tumor_Sample_Barcode", "Chromosome", "Start_Position",
            "Reference_Allele", "Tumor_Seq_Allele2")
  required <- c(keys, "Hugo_Symbol")
  if (!all(required %in% names(v))) stop("Variant counting requires: ", paste(required, collapse = ", "))
  if (anyDuplicated(sample_ids)) stop("Oncoplot requires unique sample IDs.")
  v <- v[v$Tumor_Sample_Barcode %in% sample_ids, , drop = FALSE]
  if (anyNA(v[, required])) stop("Missing variant identity fields; count cannot be verified.")
  if ("NCBI_Build" %in% names(v)) keys <- c(keys, "NCBI_Build")
  if ("End_Position" %in% names(v)) keys <- c(keys, "End_Position")
  v$Chromosome <- sub("^chr", "", as.character(v$Chromosome), ignore.case = TRUE)
  events <- v[!duplicated(v[, keys, drop = FALSE]), , drop = FALSE]
  gene_events <- v[!duplicated(v[, c(keys, "Hugo_Symbol"), drop = FALSE]), , drop = FALSE]
  sample_counts <- as.integer(table(factor(events$Tumor_Sample_Barcode, levels = sample_ids)))
  gene_counts <- as.integer(table(factor(gene_events$Hugo_Symbol, levels = genes)))
  carriers <- unique(gene_events[, c("Tumor_Sample_Barcode", "Hugo_Symbol")])
  patient_counts <- as.integer(table(factor(carriers$Hugo_Symbol, levels = genes)))
  list(
    patients = data.frame(sample_id = sample_ids, variant_count = sample_counts),
    genes = data.frame(gene = genes, variant_count = gene_counts,
      patients = patient_counts, displayed_n = length(sample_ids),
      percent = 100 * patient_counts / length(sample_ids)),
    events = events, gene_events = gene_events,
    audit = data.frame(retained_rows = nrow(v), unique_patient_events = nrow(events),
      unique_patient_gene_events = nrow(gene_events),
      removed_duplicate_patient_gene_rows = nrow(v) - nrow(gene_events))
  )
}

oncoplot_count_annotations <- function(counts, font_size = 8) {
  list(
    top = ComplexHeatmap::HeatmapAnnotation(
      `Variant count` = ComplexHeatmap::anno_barplot(counts$patients$variant_count,
        border = FALSE, gp = grid::gpar(fill = "#487990", col = NA),
        height = grid::unit(13, "mm"), axis_param = list(gp = grid::gpar(fontsize = font_size))),
      annotation_name_gp = grid::gpar(fontsize = font_size)),
    right = ComplexHeatmap::rowAnnotation(
      `Variant count` = ComplexHeatmap::anno_barplot(counts$genes$variant_count,
        border = FALSE, gp = grid::gpar(fill = "#487990", col = NA),
        width = grid::unit(22, "mm"), axis_param = list(gp = grid::gpar(fontsize = font_size))),
      annotation_name_gp = grid::gpar(fontsize = font_size)),
    left = ComplexHeatmap::rowAnnotation(
      `Patients n (%)` = ComplexHeatmap::anno_text(
        sprintf("%d (%.1f%%)", counts$genes$patients, counts$genes$percent),
        gp = grid::gpar(fontsize = font_size), just = "right",
        location = grid::unit(1, "npc"), width = grid::unit(26, "mm")),
      annotation_name_gp = grid::gpar(fontsize = font_size))
  )
}
