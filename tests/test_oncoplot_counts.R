source("R/oncoplot_counts.R")
v <- data.frame(Tumor_Sample_Barcode = c("A", "A", "A", "B", "A"),
  Hugo_Symbol = c("TP53", "TP53", "TP53", "TP53", "ATM"),
  Chromosome = "17", Start_Position = c(10, 10, 20, 10, 30),
  Reference_Allele = "C", Tumor_Seq_Allele2 = "T")
x <- oncoplot_variant_counts(v, c("A", "B", "C"), c("TP53", "ATM", "BRCA2"))
stopifnot(identical(x$patients$variant_count, c(3L, 1L, 0L)),
          identical(x$genes$variant_count, c(3L, 1L, 0L)),
          identical(x$genes$patients, c(2L, 1L, 0L)),
          x$audit$removed_duplicate_patient_gene_rows == 1)
# Same genomic event with two gene annotations is one patient event.
v <- rbind(v, transform(v[1, ], Hugo_Symbol = "OTHER"))
y <- oncoplot_variant_counts(v, c("A", "B"), c("TP53", "ATM", "OTHER"))
stopifnot(sum(y$patients$variant_count) == 4,
          sum(y$genes$variant_count) == 5)
cat("Oncoplot event / gene / patient counting tests passed.\n")
