#!/usr/bin/env Rscript
# HRD-related (18) and MMR (4) gene lists from McIntyre 2020.
# Run from PAAD root: Rscript --vanilla scripts/repair_gene_profiles.R --config=config/local.R
# No CNV/LOH/purity/cellularity analyses; no functional HRD or dMMR inference.
source("R/config.R")
source("R/data.R")
source("R/genomics.R")
source("R/plot_helpers.R")
source("R/pairwise_tests.R")
source("R/repair_gene_sets.R")
config <- load_config_from_command_line()
if (!identical(config$clinical_schema, "curated_v3")) {
  message("Repair-gene profiles skipped: curated_v3 required."); quit(status = 0)
}
input <- load_analysis_data(config)
set.seed(if (is.null(config$random_seed)) 20260928L else config$random_seed)
v <- as.data.frame(input$maf@data)
v$Tumor_Sample_Barcode <- as.character(v$Tumor_Sample_Barcode)
v$paper_gene <- repair_paper_symbol(v$Hugo_Symbol)
d <- repair_group_status(add_source_genomics(input$clinical, v), v)
sets <- repair_gene_sets()
stopifnot(!anyDuplicated(d$patient_id), !any(d$patient_duplicate))
tab <- function(x, f) readr::write_tsv(as.data.frame(x), output_path(config, f), na = "NA")
save_plot <- function(p, f, w = 12, h = 8) ggplot2::ggsave(output_path(config, f), p,
  width = w, height = h, dpi = min(config$plot_dpi, 300), bg = "white", limitsize = FALSE)
theme_set(theme_classic(base_size = 12, base_family = "Pretendard") +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9)))
list_names <- c(HRD = "HRD-related genes", MMR = "MMR genes")
definition <- bind_rows(lapply(names(sets), function(s) data.frame(
  gene_list = s, paper_gene = sets[[s]], n_genes = length(sets[[s]]),
  accepted_alias = ifelse(sets[[s]] == "FAM175A", "ABRAXAS1", NA_character_),
  source = "McIntyre 2020, Methods p3940; doi:10.1002/cncr.33038",
  interpretation = if (s == "HRD") "Gene-list variants; NOT a functional HRD diagnosis" else
    "Gene-list variants; NOT dMMR diagnosis or MSI status")))
tab(definition, "31_Repair_gene_list_definition.tsv")
tab(d[, c("Tumor_Sample_Barcode", "patient_id", "HRD_gene_count", "MMR_gene_count",
  "HRD_variant_status", "MMR_variant_status", "Repair_pattern")], "31_Repair_patient_status_private.tsv")
vf <- v[v$paper_gene %in% unlist(sets), ]
vf$gene_list <- ifelse(vf$paper_gene %in% sets$HRD, "HRD", "MMR")
freq <- bind_rows(lapply(names(sets), function(s) bind_rows(lapply(sets[[s]], function(g) {
  z <- vf[vf$paper_gene == g, ]
  data.frame(gene_list = s, gene = g, patients = n_distinct(z$Tumor_Sample_Barcode),
    total_patients = nrow(d), percent = 100 * n_distinct(z$Tumor_Sample_Barcode) / nrow(d),
    variant_rows = nrow(z), assay_coverage = "Unknown; zero is no retained call, not confirmed negative")
}))))
tab(freq, "31_Repair_gene_patient_frequencies.tsv")
classes <- vf |> count(gene_list, paper_gene, Variant_Classification, name = "variant_rows")
tab(classes, "31_Repair_variant_class_counts.tsv")
summary_rows <- bind_rows(lapply(names(sets), function(s) {
  n <- sum(d[[paste0(s, "_gene_count")]] > 0)
  data.frame(gene_list = s, total_patients = nrow(d), variant_positive_n = n,
    no_retained_variant_n = nrow(d) - n, percent = 100 * n / nrow(d),
    multiple_genes_n = sum(d[[paste0(s, "_gene_count")]] > 1))
}))
tab(summary_rows, "31_Repair_gene_list_summary.tsv")
overlap <- d |> count(Repair_pattern, name = "n", .drop = FALSE) |>
  mutate(total_patients = nrow(d), percent = 100 * n / nrow(d))
tab(overlap, "31_Repair_gene_list_overlap.tsv")
save_plot(ggplot(overlap, aes(Repair_pattern, n, fill = Repair_pattern)) +
  geom_col(show.legend = FALSE) + geom_text(aes(label = sprintf("%d (%.1f%%)", n, percent)), vjust = -.5) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Overlap of HRD-related and MMR gene-list variants", x = NULL, y = "Unique patients",
    caption = "At least one retained MAF variant defines a positive gene-list status. This is not HRD or dMMR diagnosis."),
  "31_HRD_MMR_variant_overlap.png", 11, 6)

# True mutation-class tiles: multiple classes share a cell without being hidden.
colors <- c(Missense = "#399B72", Frameshift = "#D67F37", Nonsense = "#C84659",
  Splice = "#8068A3", "In-frame" = "#D4AF37", Other = "#617C8B")
vf$class <- case_when(vf$Variant_Classification == "Missense_Mutation" ~ "Missense",
  vf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "Frameshift",
  vf$Variant_Classification == "Nonsense_Mutation" ~ "Nonsense",
  vf$Variant_Classification == "Splice_Site" ~ "Splice",
  vf$Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins") ~ "In-frame", TRUE ~ "Other")
draw_repair_oncoplot <- function(set, carriers_only = FALSE) {
  source("R/oncoplot_counts.R")
  gene_order <- freq$gene[freq$gene_list == set][order(-freq$patients[freq$gene_list == set],
    freq$gene[freq$gene_list == set])]
  z <- d[if (carriers_only) d[[paste0(set, "_gene_count")]] > 0 else rep(TRUE, nrow(d)), ]
  if (!nrow(z)) return(invisible(NULL))
  z <- z[order(-z[[paste0(set, "_gene_count")]], as.character(z$Repair_pattern), z$Tumor_Sample_Barcode), ]
  m <- matrix("", length(gene_order), nrow(z), dimnames = list(gene_order, z$Tumor_Sample_Barcode))
  sub <- vf[vf$gene_list == set & vf$Tumor_Sample_Barcode %in% colnames(m), ]
  for (i in seq_len(nrow(sub))) {
    a <- sub$paper_gene[i]; b <- sub$Tumor_Sample_Barcode[i]
    m[a, b] <- paste(unique(c(strsplit(m[a, b], ";", fixed = TRUE)[[1]], sub$class[i]))[
      nzchar(unique(c(strsplit(m[a, b], ";", fixed = TRUE)[[1]], sub$class[i])))], collapse = ";")
  }
  tracks <- c("Age", "Sex", "Neoadjuvant", "Preop_platinum_exposure", "T_stage", "N_stage",
    "M_stage", "Differentiation", "R_status", "KRAS_clinical_group", "MSI")
  ann <- z[, tracks, drop = FALSE]
  for (col in c("Neoadjuvant")) ann[[col]] <- ifelse(is.na(ann[[col]]), NA, ifelse(ann[[col]] == "y", "Yes", "No"))
  for (col in c("T_stage", "N_stage", "M_stage")) ann[[col]] <- paste0(substr(col, 1, 1), as.character(ann[[col]]))
  names(ann) <- tracks
  count_variants <- v
  count_variants$Hugo_Symbol <- count_variants$paper_gene
  counts <- oncoplot_variant_counts(count_variants, colnames(m), rownames(m))
  count_annotations <- oncoplot_count_annotations(counts, 9)
  ha <- c(count_annotations$top, make_annotation(ann, tracks))
  af <- function(x, y, w, h, v) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#F1F3F5", col = NA))
    active <- names(v)[v]
    if (length(active)) for (j in seq_along(active)) grid::grid.rect(x,
      y - h * .45 + (j - .5) * h * .9 / length(active), w * .95, h * .9 / length(active),
      gp = grid::gpar(fill = colors[[active[j]]], col = NA))
  }
  suffix <- if (carriers_only) "variant_positive_patients" else "all_patients"
  filename <- paste0("31_", set, "_oncoplot_", suffix, ".png")
  grDevices::png(output_path(config, filename), width = 15, height = if (set == "HRD") 11 else 9.5,
    units = "in", res = min(config$plot_dpi, 300))
  grid::grid.newpage(); grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = "Pretendard")))
  ht <- ComplexHeatmap::oncoPrint(m, alter_fun = af, alter_fun_is_vectorized = FALSE, col = colors,
    top_annotation = ha, right_annotation = count_annotations$right,
    left_annotation = count_annotations$left, show_pct = FALSE,
    remove_empty_rows = FALSE, remove_empty_columns = FALSE,
    column_order = seq_len(ncol(m)), row_order = seq_len(nrow(m)), show_column_names = FALSE,
    column_title = paste0(list_names[set], " | ",
      if (carriers_only) "patients with a gene-list variant" else "all unique patients",
      " (N=", nrow(z), ")\nGrey mutation cells = no retained variant; grey clinical tracks = missing. Row % uses displayed N."),
    column_title_gp = grid::gpar(fontsize = 13), row_names_gp = grid::gpar(fontsize = 12),
    pct_gp = grid::gpar(fontsize = 10), pct_digits = 1,
    heatmap_legend_param = list(title = "MAF variant class"))
  ComplexHeatmap::draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  invisible(filename)
}
for (s in names(sets)) for (carrier in c(FALSE, TRUE)) draw_repair_oncoplot(s, carrier)

# Baseline/biomarker descriptions: outcomes are not baseline clinical features.
continuous <- c("Age", "BMI", "CA19_9", "CEA", "Tumor_size", "TMB", "MAF_variant_count")
categorical <- setdiff(c(curated_categorical(), "KRAS_clinical_group", "Stage_Group"),
  c("KRAS_subtype", "Death_event", "Recurrence_event", "Recurrence_pattern", "Distant_pattern"))
stopifnot(!length(excluded_analysis_fields(c(continuous, categorical))))
clinical_summaries <- tests <- pairs <- list()
groupings <- c("HRD_variant_status", "MMR_variant_status", "Repair_pattern")
for (grouping in groupings) for (variable in c(continuous, categorical)) {
  value <- d[[variable]]; g <- droplevels(factor(d[[grouping]])); id <- paste(grouping, variable)
  for (level in levels(g)) {
    x <- value[g == level]; N <- length(x); nobs <- sum(!is.na(x))
    if (variable %in% continuous) {
      q <- if (nobs) quantile(x, c(.25, .5, .75), na.rm = TRUE, names = FALSE) else rep(NA_real_, 3)
      result <- data.frame(grouping, variable, display_label = figure_label(variable), group = level,
        category = "Continuous", total_n = N, observed_n = nobs, missing_n = N - nobs,
        n = NA_integer_, percent = NA_real_, q1 = q[1], median = q[2], q3 = q[3])
    } else {
      tt <- table(ifelse(is.na(x), "Missing / unavailable", as.character(x)))
      result <- data.frame(grouping, variable, display_label = figure_label(variable), group = level,
        category = names(tt), total_n = N, observed_n = nobs, missing_n = N - nobs,
        n = as.integer(tt), percent = 100 * as.integer(tt) / N, q1 = NA_real_, median = NA_real_, q3 = NA_real_)
    }
    clinical_summaries[[paste(id, level)]] <- result
  }
  ok <- complete.cases(value, g); x <- value[ok]; gg <- droplevels(g[ok])
  method <- if (variable %in% continuous) {
    if (nlevels(gg) == 2) "Wilcoxon rank-sum" else "Kruskal-Wallis"
  } else "Fisher exact / Monte Carlo 10000"
  test <- tryCatch(if (length(unique(x)) < 2 || nlevels(gg) < 2) stop("Insufficient variation") else
    if (variable %in% continuous) {
      if (nlevels(gg) == 2) suppressWarnings(wilcox.test(x ~ gg, exact = FALSE)) else kruskal.test(x ~ gg)
    } else fisher.test(table(x, gg), simulate.p.value = nrow(table(x, gg)) > 2 || nlevels(gg) > 2, B = 10000),
    error = identity)
  tests[[id]] <- data.frame(grouping, variable, n = sum(ok), missing_n = sum(!ok), method,
    statistic = if (!inherits(test, "error") && !is.null(test$statistic)) unname(test$statistic) else NA_real_,
    p = if (inherits(test, "error")) NA_real_ else test$p.value,
    status = if (inherits(test, "error")) conditionMessage(test) else "Estimated")
  pw <- if (variable %in% continuous) pairwise_numeric(value, g) else pairwise_categorical(value, g)
  if (nrow(pw)) pairs[[id]] <- cbind(data.frame(grouping, variable), pw)
}
summary <- bind_rows(clinical_summaries)
tt <- bind_rows(tests) |> group_by(grouping) |> mutate(q_BH = p.adjust(p, "BH")) |> ungroup()
tab(summary, "32_Repair_group_clinical_distributions.tsv")
tab(tt, "32_Repair_group_clinical_tests.tsv")
tab(bind_rows(pairs), "32_Repair_group_clinical_pairwise_tests.tsv")
for (set in names(sets)) {
  grouping <- paste0(set, "_variant_status")
  show <- c("Sex", "Neoadjuvant", "Preop_platinum_exposure", "T_stage", "N_stage", "M_stage", "Differentiation", "R_status")
  z <- summary[summary$grouping == grouping & summary$variable %in% show, ]
  z$category <- ifelse(z$variable %in% c("T_stage", "N_stage", "M_stage") & z$category != "Missing / unavailable",
    paste0(substr(z$variable, 1, 1), z$category), z$category)
  z$category[z$category == "y"] <- "Yes"; z$category[z$category == "n"] <- "No"
  z$category[z$category == "F"] <- "Female"; z$category[z$category == "M"] <- "Male"
  z$category[z$category == "MD"] <- "Moderate"; z$category[z$category == "WD"] <- "Well"
  z$category[z$category == "PD"] <- "Poor"
  save_plot(ggplot(z, aes(group, percent, fill = category)) + geom_col(show.legend = FALSE) +
    geom_text(aes(label = ifelse(percent >= 6, paste0(category, ifelse(percent < 18, " ", "\n"), sprintf("%.1f%%", percent)), "")),
      position = position_stack(vjust = .5), size = 4.5, show.legend = FALSE) +
    facet_wrap(~display_label, scales = "free", ncol = 4) +
    labs(title = paste(list_names[set], "and clinical characteristics"), x = NULL, y = "% within variant-status group",
      fill = "Recorded category",
      subtitle = paste0("No retained variant: n=", sum(d[[grouping]] == "No retained variant"),
        "; Variant detected: n=", sum(d[[grouping]] == "Variant detected")),
      caption = "Includes missing entries in denominators. Gene-list status is not HRD/dMMR diagnosis.\nSegments <6% are not labeled; all categories and BH-adjusted tests are in 32_Repair_group_clinical_distributions.tsv / tests.tsv.") +
    theme(text = element_text(size = 16), axis.text.x = element_text(size = 12),
      strip.text = element_text(size = 16), plot.caption = element_text(size = 11)),
    paste0("32_", set, "_clinical_distributions.png"), 16, 8)
  num <- d |> select(all_of(c(grouping, continuous))) |>
    pivot_longer(all_of(continuous), names_to = "variable", values_to = "value") |>
    mutate(group = .data[[grouping]])
  num_tests <- bind_rows(pairs)
  num_tests <- num_tests[num_tests$grouping == grouping & num_tests$variable %in% continuous, ]
  panel_labels <- setNames(paste0(figure_label(num_tests$variable), "\nPairwise Holm p ",
    ifelse(is.na(num_tests$p_holm), "not estimable",
      ifelse(num_tests$p_holm < .001, "<0.001", paste0("= ", format.pval(num_tests$p_holm, digits = 2))))),
    num_tests$variable)
  num$variable <- factor(panel_labels[num$variable], levels = unname(panel_labels[continuous]))
  save_plot(ggplot(num[!is.na(num$value), ], aes(group, value, fill = group)) +
    geom_boxplot(outlier.alpha = .3, show.legend = FALSE) +
    facet_wrap(~variable, scales = "free_y", ncol = 3) +
    labs(title = paste(list_names[set], "and continuous clinical measures"), x = NULL, y = NULL,
      caption = "Observed values only. Separate clinical TMB (mut/Mb) and Variant count (MAF rows). Pairwise details: 32_Repair_group_clinical_pairwise_tests.tsv.\nHolm family: one pair per clinical field; BH across fields is reported separately in tests.tsv."),
    paste0("32_", set, "_continuous_clinical_distributions.png"), 14, 10)
}
# All gene-specific clinical profiles: descriptive, retaining sparse groups.
gene_profiles <- list()
for (set in names(sets)) for (gene in sets[[set]]) {
  ids <- unique(vf$Tumor_Sample_Barcode[vf$paper_gene == gene])
  for (variable in c(continuous, categorical)) for (status in c("Variant detected", "No retained variant")) {
    keep <- if (status == "Variant detected") d$Tumor_Sample_Barcode %in% ids else !d$Tumor_Sample_Barcode %in% ids
    x <- d[[variable]][keep]; N <- length(x)
    id <- paste(gene, variable, status)
    if (variable %in% continuous) {
      q <- if (any(!is.na(x))) quantile(x, c(.25, .5, .75), na.rm = TRUE, names = FALSE) else rep(NA_real_, 3)
      gene_profiles[[id]] <- data.frame(gene_list = set, gene, status, variable, category = "Continuous",
        total_n = N, missing_n = sum(is.na(x)), n = NA_integer_, percent = NA_real_, q1 = q[1], median = q[2], q3 = q[3])
    } else if (N) {
      xx <- table(ifelse(is.na(x), "Missing / unavailable", as.character(x)))
      gene_profiles[[id]] <- data.frame(gene_list = set, gene, status, variable, category = names(xx),
        total_n = N, missing_n = sum(is.na(x)), n = as.integer(xx), percent = 100 * as.integer(xx) / N,
        q1 = NA_real_, median = NA_real_, q3 = NA_real_)
    }
  }
}
tab(bind_rows(gene_profiles), "32_Per_gene_clinical_distributions.tsv")
tab(data.frame(excluded = c("LOH", "CNV", "Tumor purity", "Tumor cellularity"),
  scope = "Excluded from model predictors, summary comparisons and plot annotations; raw input unchanged"),
  "31_Excluded_analysis_domains.tsv")
message("HRD-related/MMR oncoplots and clinical profiles complete: ", config$output_dir)
