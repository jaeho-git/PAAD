#!/usr/bin/env Rscript
# Run from the repository root:
# Rscript --vanilla scripts/source_comparisons.R --config=config/local.R
#
# 25: KRAS detailed/grouped agreement and individual discrepancy list
# 26: Clinical- and MAF-derived KRAS analyses, separately
# 27: Clinical TMB and MAF variant-count distributions, pairing and missingness
# 28: All approved clinical variables: distributions and genomic associations
# 29: Source-specific adjusted OS/recorded-recurrence models and diagnostics
# All tables/figures go to the configured, Git-ignored output directory.
source("R/config.R")
source("R/data.R")
source("R/genomics.R")
source("R/pairwise_tests.R")
config <- load_config_from_command_line()
if (!identical(config$clinical_schema, "curated_v3")) {
  message("Source comparisons skipped: curated_v3 workbook required.")
  quit(status = 0)
}
input <- load_analysis_data(config)
set.seed(if (is.null(config$random_seed)) 20260928L else config$random_seed)
d <- add_source_genomics(input$clinical, as.data.frame(input$maf@data))
stopifnot(!anyDuplicated(d$patient_id), !any(d$patient_duplicate))
for (g in c("KRAS", "TP53", "SMAD4", "CDKN2A")) {
  d[[g]] <- as.integer(d$Tumor_Sample_Barcode %in% input$maf@data$Tumor_Sample_Barcode[input$maf@data$Hugo_Symbol == g])
}
d$Age10 <- d$Age / 10
d$TMB_availability <- factor(ifelse(is.na(d$TMB), "Missing", "Observed"), levels = c("Observed", "Missing"))
continuous <- curated_continuous()
categorical <- curated_categorical()
tab <- function(x, f) readr::write_tsv(as.data.frame(x), output_path(config, f), na = "NA")
plot_save <- function(p, f, w = 11, h = 7) ggplot2::ggsave(output_path(config, f), p,
  width = w, height = h, dpi = min(config$plot_dpi, 300L), bg = "white", limitsize = FALSE)
theme_set(theme_classic(base_size = 12) + theme(legend.position = "bottom",
  plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9)))
source_labels <- c(TMB = "Clinical TMB (mutations/Mb)", MAF_variant_count = "Variant count",
  KRAS_clinical_group = "Clinical KRAS subtype", KRAS_MAF_group = "MAF-derived KRAS subtype",
  KRAS_clinical_group_G12D_vs_other = "Clinical KRAS: G12D vs other",
  KRAS_MAF_group_G12D_vs_other = "MAF KRAS: G12D vs other")
label <- function(x) { out <- figure_label(x); i <- x %in% names(source_labels); out[i] <- source_labels[x[i]]; out }
tab(data.frame(variable = c(continuous, categorical, "Age_group", "Stage_Group"),
  display_label = clinical_label(c(continuous, categorical, "Age_group", "Stage_Group")),
  figure_label = label(c(continuous, categorical, "Age_group", "Stage_Group")),
  type = c(rep("continuous", length(continuous)), rep("categorical", length(categorical) + 2L))),
  "28_Applied_variable_mapping.tsv")

# 1. KRAS: retain both original strings; normalize only a separate comparison key.
d$KRAS_clinical_normalized <- normalize_kras(d$KRAS_subtype)
d$KRAS_MAF_normalized <- normalize_kras(d$KRAS_subtype_raw)
d$KRAS_group_agreement <- as.character(d$KRAS_clinical_group) == as.character(d$KRAS_MAF_group)
both_absent <- d$KRAS_clinical_normalized == "No KRAS call" & d$KRAS_MAF_normalized == "No KRAS call"
d$KRAS_exact_comparable <- !is.na(d$KRAS_clinical_normalized) & d$KRAS_clinical_normalized != "No KRAS call" &
  d$KRAS_MAF_exact_resolved & d$KRAS_MAF_normalized != "No KRAS call"
d$KRAS_exact_agreement <- ifelse(d$KRAS_exact_comparable,
  d$KRAS_clinical_normalized == d$KRAS_MAF_normalized, NA)
d$KRAS_comparison_status <- case_when(
  is.na(d$KRAS_clinical_normalized) ~ "Clinical subtype missing",
  both_absent ~ "Both report no detection / no call",
  (d$KRAS_clinical_normalized == "No KRAS call") != (d$KRAS_MAF_normalized == "No KRAS call") ~ "Detection / no-call discordance",
  d$KRAS_exact_comparable & d$KRAS_exact_agreement ~ "Exact subtype concordant",
  d$KRAS_exact_comparable & !d$KRAS_exact_agreement ~ "Exact subtype discordant",
  !d$KRAS_group_agreement ~ "Group discordant; exact MAF unresolved",
  TRUE ~ "Group concordant; exact MAF unresolved")
private_cols <- c("Tumor_Sample_Barcode", "patient_id", "KRAS_subtype", "KRAS_subtype_raw",
  "KRAS_clinical_normalized", "KRAS_MAF_normalized", "KRAS_clinical_group", "KRAS_MAF_group",
  "KRAS_MAF_exact_resolved", "KRAS_exact_comparable", "KRAS_exact_agreement",
  "KRAS_group_agreement", "KRAS_comparison_status")
tab(d[, private_cols], "25_KRAS_source_comparison_private.tsv")
tab(as.data.frame(input$maf@data[input$maf@data$Hugo_Symbol == "KRAS", ]), "25_KRAS_MAF_variant_evidence_private.tsv")
tab(d[!d$KRAS_group_agreement | (!is.na(d$KRAS_exact_agreement) & !d$KRAS_exact_agreement) |
  grepl("unresolved|missing", d$KRAS_comparison_status), private_cols], "25_KRAS_discrepancies_and_unresolved_private.tsv")
for (detail in c("group", "detailed")) {
  a <- if (detail == "group") as.character(d$KRAS_clinical_group) else d$KRAS_clinical_normalized
  b <- if (detail == "group") as.character(d$KRAS_MAF_group) else d$KRAS_MAF_normalized
  z <- as.data.frame(table(Clinical = a, MAF = b), stringsAsFactors = FALSE)
  tab(z, paste0("25_KRAS_", detail, "_cross_table.tsv"))
  plot_save(ggplot(z, aes(MAF, Clinical, fill = Freq)) + geom_tile(color = "white") +
    geom_text(aes(label = ifelse(Freq == 0, "", Freq)), size = 3.1) +
    scale_fill_gradient(low = "#F1F5F8", high = "#3D89B0") +
    labs(title = paste("KRAS source comparison:", detail), x = "MAF-derived call", y = "Clinical report",
      fill = "Patients", caption = "Unique patients. Unresolved MAF annotations are not exact amino-acid matches. No call is not confirmed wild type.") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    paste0("25_KRAS_", detail, "_heatmap.png"), if (detail == "group") 10 else 13, if (detail == "group") 7 else 11)
}
status <- d |> count(KRAS_comparison_status, name = "n") |> mutate(percent = 100 * n / nrow(d))
tab(status, "25_KRAS_comparison_status.tsv")
plot_save(ggplot(status, aes(n, reorder(KRAS_comparison_status, n))) +
  geom_col(fill = "#367EA4") + geom_text(aes(label = n), hjust = -.2) +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "KRAS agreement and unresolved comparisons", x = "Patients", y = NULL),
  "25_KRAS_disagreement_types.png", 12, 6)
group_table <- table(d$KRAS_clinical_group, d$KRAS_MAF_group)
po <- sum(diag(group_table)) / sum(group_table)
pe <- sum(rowSums(group_table) * colSums(group_table)) / sum(group_table)^2
exact_n <- sum(d$KRAS_exact_comparable)
agreement <- data.frame(comparison = c("Collapsed category", "Exact subtype among resolvable mutation-positive pairs"),
  n = c(sum(group_table), exact_n), agreement_n = c(sum(diag(group_table)), sum(d$KRAS_exact_agreement, na.rm = TRUE)),
  agreement_fraction = c(po, if (exact_n) mean(d$KRAS_exact_agreement, na.rm = TRUE) else NA),
  kappa = c(if (pe < 1) (po - pe) / (1 - pe) else NA, NA))
tab(agreement, "25_KRAS_agreement_summary.tsv")
kras_freq <- bind_rows(lapply(c("KRAS_clinical_group", "KRAS_MAF_group"), function(v) {
  data.frame(source = label(v), subtype = as.character(d[[v]])) |> count(source, subtype, name = "n") |>
    mutate(denominator = nrow(d), percent = 100 * n / denominator)
}))
tab(kras_freq, "26_KRAS_group_frequencies.tsv")
for (v in c("KRAS_subtype", "KRAS_subtype_raw")) tab(d |> count(.data[[v]], name = "n") |>
  mutate(denominator = nrow(d), percent = 100 * n / denominator), paste0("26_", v, "_detailed_frequencies.tsv"))
plot_save(ggplot(kras_freq, aes(subtype, percent, fill = source)) + geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(.9), vjust = -.3, size = 3.5) +
  scale_fill_manual(values = c("#30799E", "#C4783B")) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "KRAS subtype distributions by source", x = NULL, y = "% of unique patients", fill = NULL,
    caption = "Other KRAS includes rare/multiple variants. Clinical WT/Not detected and MAF no call share a comparison category only."),
  "26_KRAS_source_distributions.png")

# 2. Kaplan-Meier analyses for each KRAS source. No arbitrary TMB split.
# Recurrence follow-up is undefined in M1 and remains missing, never zero.
km_tests <- km_medians <- km_risks <- km_pairs <- list()
km <- function(data, variable, endpoint, file) {
  time <- if (endpoint == "OS") "OS_months" else "RFS_m"
  event <- if (endpoint == "OS") "survive" else "Recur"
  z <- data[complete.cases(data[, c(variable, time, event)]), ]
  z$group <- droplevels(factor(z[[variable]])); z$time <- z[[time]]; z$event <- z[[event]]
  if (nlevels(z$group) < 2 || sum(z$event) == 0) return(invisible(NULL))
  f <- survival::survfit(survival::Surv(time, event) ~ group, data = z, conf.type = "log-log")
  lr <- survival::survdiff(survival::Surv(time, event) ~ group, data = z)
  p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  id <- paste(variable, endpoint)
  pw <- pairwise_logrank(z$time, z$event, z$group)
  km_pairs[[id]] <<- cbind(data.frame(variable, endpoint), pw)
  km_tests[[id]] <<- data.frame(variable, endpoint, n = nrow(z), events = sum(z$event), p)
  s <- summary(f, censored = TRUE)
  curve <- data.frame(time = s$time, probability = s$surv, lower = s$lower, upper = s$upper,
    censor = s$n.censor, group = sub("^group=", "", as.character(s$strata))) |>
    bind_rows(data.frame(time = 0, probability = 1, lower = 1, upper = 1, censor = 0, group = levels(z$group))) |>
    arrange(group, time)
  med <- as.data.frame(summary(f)$table)
  km_medians[[id]] <<- data.frame(variable, endpoint, group = sub("^group=", "", rownames(med)),
    n = med$records, events = med$events, median = med$median, lower95 = med[["0.95LCL"]], upper95 = med[["0.95UCL"]])
  times <- seq(0, max(z$time), 12)
  risk <- tidyr::expand_grid(group = levels(z$group), time = times) |>
    rowwise() |> mutate(n_risk = sum(as.character(z$group) == group & z$time >= time)) |> ungroup()
  km_risks[[id]] <<- mutate(risk, variable = variable, endpoint = endpoint)
  curve$group <- factor(curve$group, levels = levels(z$group))
  cols <- setNames(c("#777777", "#377DA4", "#C87837", "#4C916B", "#8B68A1")[seq_len(nlevels(z$group))], levels(z$group))
  p1 <- ggplot(curve, aes(time, probability, color = group)) + geom_step(linewidth = .8) +
    geom_point(data = curve[curve$censor > 0, ], shape = 3, size = 1) +
    scale_color_manual(values = cols) + scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = times, limits = c(0, max(z$time))) +
    labs(title = paste(label(variable), if (endpoint == "OS") "overall survival" else "recorded recurrence"), subtitle = sprintf("N=%d; events=%d; log-rank p=%.3g", nrow(z), sum(z$event), p),
      x = "Months from surgery", y = if (endpoint == "OS") "Overall survival" else "Recurrence-only event-free probability",
      color = NULL, caption = if (endpoint == "OS") "Unadjusted association; not a causal effect." else
        "Recorded recurrence only. Deaths are not added as events. This is not CIF or death-inclusive DFS/RFS.")
  if (nlevels(z$group) <= 2) p1 <- p1 + geom_ribbon(aes(ymin = lower, ymax = upper, fill = group),
    color = NA, alpha = .1, show.legend = FALSE) + scale_fill_manual(values = cols)
  p2 <- ggplot(risk, aes(time, factor(group, levels = rev(levels(z$group))), label = n_risk)) + geom_text(size = 3.3) +
    scale_x_continuous(breaks = times, limits = c(0, max(z$time))) +
    labs(title = "Number at risk", x = "Months from surgery", y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  p1 <- p1 + labs(caption = paste(p1$labels$caption, pairwise_caption(pw), sep = "\n"))
  g1 <- ggplotGrob(p1); g2 <- ggplotGrob(p2); w <- grid::unit.pmax(g1$widths, g2$widths); g1$widths <- w; g2$widths <- w
  png(output_path(config, file), width = 12, height = 9, units = "in", res = min(config$plot_dpi, 300L))
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = c(.73, .27))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1)); grid::grid.draw(g1); grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2)); grid::grid.draw(g2); grid::upViewport(2); dev.off()
}
for (v in c("KRAS_clinical_group", "KRAS_MAF_group")) {
  d[[paste0(v, "_G12D_vs_other")]] <- factor(ifelse(is.na(d[[v]]) | d[[v]] == "No KRAS call", NA,
    ifelse(d[[v]] == "G12D", "G12D", "Other KRAS")), levels = c("Other KRAS", "G12D"))
  for (e in c("OS", "Recorded_recurrence")) for (g in c(v, paste0(v, "_G12D_vs_other"))) {
    km(d, g, e, paste0("26_", g, "_", e, ".png"))
  }
}
kt <- bind_rows(km_tests) |> group_by(endpoint) |> mutate(q_BH = p.adjust(p, "BH")) |> ungroup()
tab(kt, "26_KRAS_logrank_tests.tsv")
tab(bind_rows(km_pairs), "26_KRAS_pairwise_logrank.tsv")
tab(bind_rows(km_medians), "26_KRAS_KM_medians.tsv")
tab(bind_rows(km_risks), "26_KRAS_KM_number_at_risk.tsv")

# 3. TMB: no callable Mb denominator is present in the supplied inputs.
# Therefore MAF_variant_count is NOT re-labelled as normalized TMB_MAF.
pair <- d[complete.cases(d[, c("TMB", "MAF_variant_count")]), ]
numeric_summary <- function(x) {
  observed <- x[!is.na(x)]
  q <- if (length(observed)) quantile(observed, c(0, .25, .5, .75, 1), names = FALSE) else rep(NA_real_, 5)
  data.frame(total_n = length(x), observed_n = length(observed), missing_n = sum(is.na(x)),
    mean = if (length(observed)) mean(observed) else NA, sd = sd(observed),
    minimum = q[1], q1 = q[2], median = q[3], q3 = q[4], maximum = q[5])
}
metrics <- c("TMB", "MAF_variant_count")
ns <- bind_rows(lapply(c("All unique patients", "Paired patients"), function(cohort) {
  z <- if (cohort == "Paired patients") pair else d
  bind_rows(lapply(metrics, function(v) cbind(data.frame(cohort, metric = v, unit = label(v)), numeric_summary(z[[v]]))))
}))
tab(ns, "27_TMB_distribution_summary.tsv")
tab(d[, c("Tumor_Sample_Barcode", "patient_id", "TMB", "MAF_variant_count", "TMB_availability",
  "MSI", "NGS_group", "M_stage")], "27_TMB_MAF_pairs_and_missingness_private.tsv")
rho <- if (nrow(pair) >= 3) cor(pair$TMB, pair$MAF_variant_count, method = "spearman") else NA_real_
boot <- if (nrow(pair) >= 3) replicate(2000, {
  ii <- sample.int(nrow(pair), replace = TRUE)
  suppressWarnings(cor(pair$TMB[ii], pair$MAF_variant_count[ii], method = "spearman"))
}) else NA_real_
ci <- if (any(is.finite(boot))) quantile(boot, c(.025, .975), na.rm = TRUE, names = FALSE) else c(NA, NA)
tab(data.frame(n = nrow(pair), spearman_rho = rho, lower95 = ci[1], upper95 = ci[2],
  bootstrap_replicates = 2000, method = "Patient-paired nonparametric percentile bootstrap; association, NOT agreement"),
  "27_TMB_MAF_rank_correlation.tsv")
plot_save(ggplot(pair, aes(MAF_variant_count, TMB)) + geom_point(alpha = .35, color = "#307A9E") +
  labs(title = "Clinical TMB versus Variant count",
    subtitle = sprintf("Paired N=%d; Spearman rho=%.3f (bootstrap 95%% CI %.3f to %.3f)", nrow(pair), rho, ci[1], ci[2]),
    x = label("MAF_variant_count"), y = label("TMB"),
    caption = "Different units: no identity line or Bland-Altman limits. Correlation does not establish interchangeability.\nShared assay/pipeline provenance is unconfirmed. Missing clinical TMB is not filled."),
  "27_TMB_MAF_paired_scatter.png")
long <- d |> select(all_of(metrics)) |> pivot_longer(everything(), names_to = "metric", values_to = "value") |>
  mutate(metric = label(metric))
plot_save(ggplot(long[!is.na(long$value), ], aes(value)) + geom_histogram(bins = 25, fill = "#367FA4", color = "white") +
  facet_wrap(~metric, scales = "free", ncol = 2) +
  labs(title = "Source-specific burden distributions", x = "Value in the panel's stated units", y = "Patients",
    caption = sprintf("Clinical TMB observed %d/%d; Variant count observed %d/%d. Units and denominators differ.",
      sum(!is.na(d$TMB)), nrow(d), sum(!is.na(d$MAF_variant_count)), nrow(d))),
  "27_TMB_source_distributions.png", 12, 6)
paired_long <- pair |> select(all_of(metrics)) |> pivot_longer(everything(), names_to = "metric", values_to = "value") |>
  mutate(metric = label(metric))
plot_save(ggplot(paired_long, aes(value)) + stat_ecdf(geom = "step", color = "#367FA4") +
  facet_wrap(~metric, scales = "free_x") + scale_y_continuous(labels = scales::percent) +
  labs(title = paste("Burden distributions in the same", nrow(pair), "patients"), x = "Value in the panel's stated units",
    y = "Empirical cumulative proportion"), "27_TMB_paired_ECDF.png", 12, 6)
heat <- pair[order(pair$MAF_variant_count, pair$TMB), metrics]
heat$patient_order <- seq_len(nrow(heat))
heat <- heat |> pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
  group_by(metric) |> mutate(z_score = as.numeric(scale(value))) |> ungroup()
plot_save(ggplot(heat, aes(patient_order, label(metric), fill = z_score)) + geom_tile() +
  scale_fill_gradient2(low = "#336C9E", mid = "white", high = "#B85446", midpoint = 0) +
  labs(title = "Paired burden profiles (each source standardized separately)", x = "Patients ordered by Variant count then clinical TMB",
    y = NULL, fill = "Within-source z",
    caption = "Standardization is for visualization only. It does not convert counts to mutations/Mb or prove agreement."),
  "27_TMB_paired_heatmap.png", 12, 4)
tab(data.frame(analysis = c("TMB_MAF (mutations/Mb)", "Bland-Altman / absolute agreement", "Missing clinical TMB replacement"),
  status = c("Not calculated", "Not estimated", "Not performed"),
  reason = c("Validated callable panel Mb and assay filtering/calibration are unavailable. A denominator is never inferred from correlation.",
    "The two supplied measures have different units; normalized comparable measurements are needed.",
    "Clinical TMB missingness is preserved; source choice remains a researcher decision.")), "27_TMB_analysis_limits.tsv")

# 4. New clinical fields are analysed, not merely carried into the dataset.
# Descriptive figures include all 40 fields except identifiers (not meaningful
# clinical plots). Tests never treat outcomes/adjuvant treatment as baseline
# covariates in the adjusted survival models.
desc <- bind_rows(lapply(continuous, function(v) cbind(data.frame(variable = v, display_label = label(v)), numeric_summary(d[[v]]))))
tab(desc, "28_Continuous_clinical_summary.tsv")
cats <- bind_rows(lapply(c(categorical, "Age_group", "Stage_Group"), function(v) {
  data.frame(variable = v, display_label = label(v), category = ifelse(is.na(d[[v]]), "Missing / unavailable", as.character(d[[v]]))) |>
    count(variable, display_label, category, name = "n") |> mutate(denominator = nrow(d), percent = 100 * n / denominator)
}))
tab(cats, "28_Categorical_clinical_summary.tsv")
cl <- d |> select(all_of(continuous)) |> pivot_longer(everything(), names_to = "variable", values_to = "value")
plot_save(ggplot(cl[!is.na(cl$value), ], aes(value)) + geom_histogram(bins = 25, fill = "#367FA4", color = "white") +
  facet_wrap(~variable, scales = "free", ncol = 3, labeller = labeller(variable = figure_label)) +
  labs(title = "Continuous clinical variables", x = "Recorded value", y = "Patients",
    caption = "No arbitrary high/low cutoffs. Missing counts and units are listed in the companion table."),
  "28_Continuous_clinical_distributions.png", 14, 10)
cat_sets <- split(setdiff(c(categorical, "Age_group", "Stage_Group"), "KRAS_subtype"),
  ceiling(seq_along(setdiff(c(categorical, "Age_group", "Stage_Group"), "KRAS_subtype")) / 6))
for (i in seq_along(cat_sets)) {
  z <- cats[cats$variable %in% cat_sets[[i]], ]
  plot_save(ggplot(z, aes(n, category)) + geom_col(fill = "#367FA4") +
    geom_text(aes(label = n), hjust = -.1, size = 3) +
    facet_wrap(~display_label, scales = "free", ncol = 2) +
    scale_x_continuous(expand = expansion(mult = c(0, .2))) +
    labs(title = "Clinical categories and missingness", x = "Patients", y = NULL,
      caption = "Blank clinical entries remain missing, including Adjuvant_RT. Original MSI wording is retained."),
    paste0("28_Categorical_distributions_", i, ".png"), 13, 11)
}
# All comparisons use complete pairs; report the denominator and dropped rows.
association <- function(value, group, variable, grouping) {
  ok <- complete.cases(value, group); x <- value[ok]; g <- droplevels(factor(group[ok]))
  note <- ""; method <- if (is.numeric(value)) "Kruskal-Wallis" else "Fisher Monte Carlo (10000 replicates)"
  test <- tryCatch({
    if (nlevels(g) < 2 || length(unique(x)) < 2) stop("Fewer than two observed categories/values")
    if (is.numeric(value)) kruskal.test(x, g) else fisher.test(table(x, g), simulate.p.value = TRUE, B = 10000)
  }, error = function(e) { note <<- conditionMessage(e); NULL })
  data.frame(variable, grouping, n = sum(ok), missing_n = sum(!ok), p = if (is.null(test)) NA_real_ else test$p.value, method, note)
}
clinical_features <- setdiff(c(continuous, categorical, "Age_group", "Stage_Group"),
  c("TMB", "KRAS_subtype", "OS_months", "DFS_months", "Death_event", "Recurrence_event", "Recurrence_pattern", "Distant_pattern"))
associations <- list(); strata <- list(); clinical_pairs <- list()
for (v in clinical_features) for (g in c("KRAS_clinical_group", "KRAS_MAF_group", "TMB_availability", "KRAS", "TP53", "SMAD4", "CDKN2A")) {
  id <- paste(v, g)
  associations[[id]] <- association(d[[v]], d[[g]], v, g)
  pw <- if (is.numeric(d[[v]])) pairwise_numeric(d[[v]], d[[g]]) else pairwise_categorical(d[[v]], d[[g]])
  if (nrow(pw)) clinical_pairs[[id]] <- cbind(data.frame(variable = v, grouping = g), pw)
  strata[[id]] <- bind_rows(lapply(unique(as.character(d[[g]][!is.na(d[[g]])])), function(level) {
    z <- d[!is.na(d[[g]]) & as.character(d[[g]]) == level, ]
    if (is.numeric(z[[v]])) return(cbind(data.frame(variable = v, grouping = g, group = level), numeric_summary(z[[v]])))
    data.frame(variable = v, grouping = g, group = level, category = ifelse(is.na(z[[v]]), "Missing / unavailable", as.character(z[[v]]))) |>
      count(variable, grouping, group, category, name = "n") |> mutate(denominator = nrow(z), percent = 100 * n / denominator)
  }))
}
assoc <- bind_rows(associations) |> group_by(grouping) |> mutate(q_BH = p.adjust(p, "BH")) |> ungroup()
tab(assoc, "28_Clinical_genomic_association_tests.tsv")
tab(bind_rows(clinical_pairs), "28_Clinical_pairwise_tests.tsv")
tab(bind_rows(strata), "28_Clinical_source_stratified_tables.tsv")
tab(assoc[assoc$grouping == "TMB_availability", ], "27_TMB_missingness_clinical_comparisons.tsv")
plot_save(ggplot(assoc, aes(label(grouping), label(variable), fill = -log10(pmax(q_BH, 1e-12)))) +
  geom_tile(color = "white") + geom_text(aes(label = ifelse(is.na(q_BH), "NA", sprintf("%.2f", q_BH))), size = 2.7) +
  scale_fill_gradient(low = "#F1F5F8", high = "#30799E", na.value = "grey80") +
  labs(title = "Clinical variables and genomic measures", x = NULL, y = NULL, fill = "-log10(BH q)",
    caption = "Cells show BH q within each grouping family. Exploratory associations; not causal treatment effects.") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)),
  "28_Clinical_genomic_associations.png", 13, 13)
# TMB source-specific clinical comparisons, including each newly added field.
burden_tests <- list(); burden_pairs <- list()
burden_categories <- setdiff(clinical_features, c(continuous, "Age_group"))
for (metric in metrics) {
  for (v in clinical_features) {
    if (is.numeric(d[[v]])) {
      ok <- complete.cases(d[, c(metric, v)])
      ct <- tryCatch(cor.test(d[[metric]][ok], d[[v]][ok], method = "spearman", exact = FALSE), error = function(e) NULL)
      burden_tests[[paste(metric, v)]] <- data.frame(variable = v, metric, n = sum(ok),
        rho = if (is.null(ct)) NA_real_ else unname(ct$estimate), p = if (is.null(ct)) NA_real_ else ct$p.value,
        method = "Spearman", note = if (is.null(ct)) "Not estimable" else "")
    } else {
      pw <- pairwise_numeric(d[[metric]], d[[v]])
      if (nrow(pw)) burden_pairs[[paste(metric, v)]] <- cbind(data.frame(metric, variable = v), pw)
      a <- association(d[[metric]], d[[v]], v, metric)
      burden_tests[[paste(metric, v)]] <- transmute(a, variable, metric = grouping, n, rho = NA_real_, p, method, note)
    }
  }
  batches <- split(burden_categories, ceiling(seq_along(burden_categories) / 6))
  for (i in seq_along(batches)) {
    boxes <- bind_rows(lapply(batches[[i]], function(v) data.frame(variable = label(v), category = as.character(d[[v]]), value = d[[metric]])))
    boxes <- boxes[complete.cases(boxes), ]
    plot_save(ggplot(boxes, aes(category, value)) + geom_boxplot(fill = "#BDDAE8", outlier.alpha = .2) +
      facet_wrap(~variable, scales = "free_x", ncol = 2) +
      labs(title = paste(label(metric), "by clinical group"), x = NULL, y = label(metric),
        caption = "Observed values only; missing categories excluded. All pairs and Holm p: 27_Burden_pairwise_tests.tsv.") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1)),
      paste0("27_", metric, "_clinical_groups_", i, ".png"), 13, 11)
  }
}
bt <- bind_rows(burden_tests) |> group_by(metric) |> mutate(q_BH = p.adjust(p, "BH")) |> ungroup()
tab(bt, "27_TMB_clinical_associations.tsv")
tab(bind_rows(burden_pairs), "27_Burden_pairwise_tests.tsv")
missing_maf <- association(d$MAF_variant_count, d$TMB_availability, "MAF_variant_count", "TMB_availability")
tab(missing_maf, "27_TMB_missingness_MAF_count_comparison.tsv")
plot_save(ggplot(d, aes(TMB_availability, MAF_variant_count, fill = TMB_availability)) + geom_boxplot(outlier.alpha = .3) +
  scale_fill_manual(values = c("#367FA4", "#C87837")) +
  labs(title = "Variant count in patients with observed or missing clinical TMB", x = "Clinical TMB", y = label("MAF_variant_count"),
    subtitle = sprintf("Observed n=%d; missing n=%d; Kruskal-Wallis p=%.3g", nrow(pair), nrow(d) - nrow(pair), missing_maf$p),
    fill = NULL), "27_TMB_missingness_MAF_counts.png")

# 5. Continuous TMB and source-specific KRAS Cox models.
# Main adjustment includes M stage because M1 is INCLUDED. Recurrence analyses
# drop invariant M stage after complete-case selection. No adjuvant causal model.
coef_rows <- diagnostic_rows <- ph_rows <- cox_pairs <- list()
fit_source_cox <- function(data, exposure, endpoint, cohort, standardized = FALSE) {
  time <- if (endpoint == "OS") "OS_months" else "RFS_m"
  event <- if (endpoint == "OS") "survive" else "Recur"
  adjust <- c("Age10", "Sex", "T_stage", "Neoadjuvant", "M_stage")
  z <- droplevels(data[complete.cases(data[, c(exposure, time, event, adjust)]), ])
  id <- paste(exposure, endpoint, cohort, if (standardized) "per paired SD" else "native scale", sep = ": ")
  z$time <- z[[time]]; z$event <- z[[event]]; z$exposure <- z[[exposure]]
  scale_value <- if (standardized) sd(pair[[exposure]]) else 1
  if (standardized) z$exposure <- z$exposure / scale_value
  adjust <- adjust[vapply(z[, adjust, drop = FALSE], function(x) length(unique(x)) > 1, logical(1))]
  for (v in intersect(c("Sex", "T_stage", "Neoadjuvant", "M_stage"), adjust)) z[[v]] <- factor(z[[v]])
  form <- reformulate(c("exposure", adjust), response = "survival::Surv(time, event)")
  warnings <- character()
  f <- tryCatch(withCallingHandlers(survival::coxph(form, data = z, x = TRUE, model = TRUE, ties = "efron"),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
  if (inherits(f, "error")) {
    diagnostic_rows[[id]] <<- data.frame(model = id, exposure, endpoint, cohort, n = nrow(z), events = sum(z$event),
      scale_value, formula = paste(deparse(form), collapse = " "), status = conditionMessage(f))
    return(invisible(NULL))
  }
  pw <- pairwise_cox(f, z, "exposure")
  if (nrow(pw)) cox_pairs[[id]] <<- cbind(data.frame(model = id, exposure, endpoint, cohort), pw)
  s <- summary(f); epv <- sum(z$event) / length(coef(f))
  status <- if (length(warnings) || any(!is.finite(s$conf.int)) || anyNA(coef(f))) "Unstable" else if (epv < 10) "Sparse: events/parameter <10" else "Estimated"
  ph <- tryCatch(survival::cox.zph(f), error = function(e) NULL)
  if (!is.null(ph)) ph_rows[[id]] <<- data.frame(model = id, term = rownames(ph$table), p = ph$table[, "p"])
  nonlinear_p <- NA_real_
  if (is.numeric(z$exposure) && length(unique(z$exposure)) >= 5 && sum(z$event) >= 30) {
    nonlinear <- tryCatch(suppressWarnings(survival::coxph(
      reformulate(c("splines::ns(exposure, df = 3)", adjust), response = "survival::Surv(time, event)"), data = z)), error = function(e) NULL)
    if (!is.null(nonlinear)) nonlinear_p <- pchisq(2 * (as.numeric(logLik(nonlinear)) - as.numeric(logLik(f))), df = 2, lower.tail = FALSE)
  }
  diagnostic_rows[[id]] <<- data.frame(model = id, exposure, endpoint, cohort, n = nrow(z), events = sum(z$event),
    scale_value, parameters = length(coef(f)), events_per_parameter = epv, formula = paste(deparse(form), collapse = " "),
    status, warnings = paste(warnings, collapse = "; "), nonlinearity_p = nonlinear_p)
  idx <- f$assign[["exposure"]]
  coef_rows[[id]] <<- data.frame(model = id, exposure, endpoint, cohort, scale = if (standardized) "per paired SD" else "native scale",
    scale_value, term = rownames(s$coefficients)[idx], n = nrow(z), events = sum(z$event),
    HR = s$conf.int[idx, "exp(coef)"], lower95 = s$conf.int[idx, "lower .95"], upper95 = s$conf.int[idx, "upper .95"],
    p = s$coefficients[idx, "Pr(>|z|)"], PH_p = if (!is.null(ph)) ph$table["exposure", "p"] else NA_real_,
    nonlinearity_p = nonlinear_p, status)
}
for (e in c("OS", "Recorded_recurrence")) {
  for (v in c("KRAS_clinical_group", "KRAS_MAF_group",
    "KRAS_clinical_group_G12D_vs_other", "KRAS_MAF_group_G12D_vs_other")) fit_source_cox(d, v, e, "All available")
  for (v in metrics) for (cohort in c("All available", "Same paired patients")) {
    z <- if (cohort == "All available") d else pair
    for (standard in c(FALSE, TRUE)) fit_source_cox(z, v, e, cohort, standard)
  }
}
tab(bind_rows(cox_pairs), "29_Source_Cox_pairwise_contrasts.tsv")
coeff <- bind_rows(coef_rows) |> group_by(endpoint, cohort, scale) |> mutate(q_BH = p.adjust(p, "BH")) |> ungroup()
tab(coeff, "29_Source_specific_Cox_coefficients.tsv")
tab(bind_rows(diagnostic_rows), "29_Source_specific_Cox_diagnostics.tsv")
tab(bind_rows(ph_rows), "29_Source_specific_Cox_PH_tests.tsv")
forest <- coeff[coeff$exposure %in% metrics & coeff$scale == "per paired SD" & coeff$status == "Estimated", ]
forest$display <- paste(label(forest$exposure), forest$cohort, sep = " / ")
plot_save(ggplot(forest, aes(HR, display, xmin = lower95, xmax = upper95, color = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") + geom_errorbar(orientation = "y", width = .15) +
  geom_point(size = 2.5) + facet_wrap(~endpoint, ncol = 1) + scale_x_log10() +
  labs(title = "TMB source sensitivity: adjusted outcome associations", x = "HR per source-specific paired-cohort SD (95% CI)",
    y = NULL, color = NULL, caption = "Adjusted for age/10, sex, T category, neoadjuvant treatment and variable M stage.\nCheck PH/nonlinearity diagnostics before interpreting a single linear HR. Pairing removes cohort-composition differences, not measurement differences."),
  "29_TMB_source_outcome_sensitivity.png", 14, 8)
kforest <- coeff[grepl("^KRAS", coeff$exposure) & coeff$status == "Estimated", ]
kforest$display <- paste(label(kforest$exposure), sub("^exposure", "", kforest$term), sep = ": ")
plot_save(ggplot(kforest, aes(HR, display, xmin = lower95, xmax = upper95)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") + geom_errorbar(orientation = "y", width = .15) +
  geom_point(size = 2, color = "#367FA4") + facet_wrap(~endpoint, ncol = 1) + scale_x_log10() +
  labs(title = "KRAS source-specific adjusted outcome associations", x = "Hazard ratio (95% CI)", y = NULL,
    caption = "Subtype models reference no call; G12D models reference other KRAS. Sources are analysed separately.\nUnstable/sparse fits are not plotted; all attempted models and PH tests are retained in tables."),
  "29_KRAS_source_adjusted_outcomes.png", 14, 12)
tab(d, "29_Integrated_analysis_dataset_private.tsv")

# 6. Reading guide is last, after every result has been written.
files <- sort(list.files(config$output_dir, pattern = "[.](png|tsv|xlsx)$"))
guide <- c("# 새 임상정보·확정 매핑 적용 결과", "",
  sprintf("- 주 분석: 중복 환자 ID의 모든 행을 제외한 %d명. M0=%d, M1=%d. M1도 OS와 주 분석에 포함.", nrow(d), sum(d$M_stage == "0"), sum(d$M_stage == "1")),
  sprintf("- OS 사망 %d건. 기록된 재발 분석 %d명, 재발 %d건. M1 재발은 정의되지 않아 이 endpoint에서만 제외.", sum(d$survive), sum(complete.cases(d[, c("RFS_m", "Recur")])), sum(d$Recur, na.rm = TRUE)),
  sprintf("- 임상 TMB 관측 %d명, 결측 %d명; MAF와 동일 환자 비교 %d쌍.", sum(!is.na(d$TMB)), sum(is.na(d$TMB)), nrow(pair)),
  sprintf("- KRAS 축약 범주 불일치 %d명. 정확한 아미노산 비교 가능 %d쌍; 비교 가능 쌍 중 불일치 %d쌍.", sum(!d$KRAS_group_agreement, na.rm = TRUE), exact_n, sum(!d$KRAS_exact_agreement, na.rm = TRUE)), "",
  "## 먼저 볼 결과", "",
  "- 25_KRAS_group_heatmap.png / 25_KRAS_detailed_heatmap.png: 축약 범주와 세부 문자열 비교.",
  "- 25_KRAS_discrepancies_and_unresolved_private.tsv: 불일치 또는 상세 판정 불가 환자별 검토.",
  "- 26_*.png: 임상 KRAS 및 MAF KRAS 각각의 빈도, OS, 기록된 재발 KM.",
  "- 27_TMB_MAF_paired_scatter.png, 27_TMB_distribution_summary.tsv: 서로 다른 단위의 두 지표.",
  "- Pairwise tables: all pairs, raw p and within-family Holm p; significant Holm pairs appear on KM curves.",
  "- Variant count = retained MAF rows per sample, not mutations/Mb.",
  "- 27_TMB_missingness_*.tsv/png: 임상 TMB 결측군의 특성과 MAF count 비교.",
  "- 28_*.png/tsv: 신규 변수 포함 임상 분포, 유전자·KRAS 출처별 연관성, 적용 자료형/약어.",
  "- 29_*.png/tsv: KRAS 출처별 보정 모형, TMB 전체 관측군 및 동일 paired cohort 민감도, PH·비선형성 진단.",
  "- 18–24: 확장 Table 1, Neoadjuvant 층화 CMH, driver/TP53/HRR 등 참고논문 분석.",
  "- 1–17: 기존 분석 방법을 새 임상정보·전체 고유 환자에 적용. TMB/RFS라는 옛 파일명은 호환용이며 제목에 실제 지표를 표시.", "",
  "## 해석 및 보류 항목", "",
  "- KRAS_subtype는 임상 원문, KRAS_subtype_raw는 MAF 파생 세부 문자열. 분석용 축약 범주는 별도 열이며 원문을 덮어쓰지 않음.",
  "- Q61 등 좌표만으로 정확히 확정하지 못한 MAF 변이는 unresolved로 표시. no call을 생물학적 WT로 단정하지 않음.",
  "- TMB는 병리과 mutations/Mb, MAF_variant_count는 retained MAF 행 수. 검증된 Mb/필터 정보가 없어 TMB_MAF, Bland–Altman 및 절대 일치도는 계산하지 않음.",
  "- 높은 상관관계만으로 두 지표의 교환 가능성이나 독립 검증을 주장할 수 없음. 결측 TMB 자동 대체 없음.",
  "- 생존모형은 탐색적 연관성. PH/비선형성 또는 희소성 경고를 확인. 치료효과·인과효과를 의미하지 않음.",
  "- 재발 전 사망 포함 DFS/RFS, CIF/Gray/Fine–Gray, LOH/CNA/생식세포/기능적 HRD 분석은 필요한 자료가 없어 보류.",
  "- Adjuvant_RT의 빈칸은 No가 아님. MSI 원문 범주 유지. Stage_Group에 IV 포함. 자료형과 표기 약어는 28_Applied_variable_mapping.tsv 참고.", "",
  "## 다시 실행", "", "```sh", "Rscript --vanilla run_all.R --config=config/local.R",
  "Rscript --vanilla tests/verify_curated_outputs.R --config=config/local.R", "```", "",
  "결과·원본자료·환자별 표는 Git 제외 대상입니다. 원본 XLSX/MAF와 작성하신 매핑 XLSX는 변경하지 않았습니다.", "",
  "## 파일 목록", "", paste0("- [", files, "](", files, ")"))
writeLines(enc2utf8(guide), output_path(config, "00_README_results.md"), useBytes = TRUE)
message("Source comparisons complete: ", config$output_dir)
