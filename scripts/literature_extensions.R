#!/usr/bin/env Rscript
# Additional, directly readable analyses inspired by McIntyre (2020) and
# Campbell (2025). Run from PAAD root with --config=config/local.R.
# Main cohort: unique-patient M0+M1 records. Full source rows are also described.
# Methods and departures from the papers: docs/updated_analysis.md.
# Tables are UTF-8 TSV so the statistical calculations remain reproducible in R.
source("R/config.R")
source("R/data.R")
source("R/genomics.R")
source("R/pairwise_tests.R")
source("R/repair_gene_sets.R")
config <- load_config_from_command_line()
if (!identical(config$clinical_schema, "curated_v3")) {
  message("Literature extensions skipped: curated_v3 workbook required.")
  quit(status = 0)
}
analysis_data <- load_analysis_data(config)
set.seed(if (is.null(config$random_seed)) 20260928L else config$random_seed)
all_clinical <- analysis_data$clinical_all
full_maf <- maftools::read.maf(config$maf_file, clinicalData = all_clinical, verbose = FALSE)
variants <- as.data.frame(full_maf@data)
genes <- c("KRAS", "TP53", "SMAD4", "CDKN2A")
save_table <- function(x, name) readr::write_tsv(as.data.frame(x), output_path(config, name), na = "NA")
save_plot <- function(p, name, width = 10, height = 7) {
  ggplot2::ggsave(output_path(config, name), p, width = width, height = height,
    dpi = min(config$plot_dpi, 400L), bg = "white", limitsize = FALSE)
}
theme_set(theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold"),
  legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9)))
palette <- c("#276F9E", "#CA6B31", "#498C66", "#8E66A3", "#858585")

# 1. MAF-derived status. No call does not imply confirmed biological wild type.
annotate_genomics <- function(d, v) {
  for (g in genes) d[[g]] <- as.integer(d$Tumor_Sample_Barcode %in% v$Tumor_Sample_Barcode[v$Hugo_Symbol == g])
  d$driver_count <- rowSums(d[, genes])
  d$driver_group <- factor(ifelse(d$driver_count >= 3, "3+", as.character(d$driver_count)), levels = c("1", "0", "2", "3+"))
  d$Any_driver <- as.integer(d$driver_count > 0)
  d <- add_source_genomics(d, v)
  d$KRAS_MAF <- d$KRAS_MAF_group
  d$KRAS_report_group <- as.character(d$KRAS_clinical_group)
  tp <- v[v$Hugo_Symbol == "TP53", ] |> group_by(Tumor_Sample_Barcode) |>
    summarise(TP53_variant_rows = n(), TP53_class = if (n() != 1) "Multiple calls (excluded)" else
      if (first(Variant_Classification) %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins")) "Truncating" else
      if (first(Variant_Classification) %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins")) "Nontruncating" else "Other / splice (excluded)", .groups = "drop")
  d <- left_join(d, tp, by = "Tumor_Sample_Barcode")
  d$TP53_class[is.na(d$TP53_class)] <- "No TP53 call"
  d$TP53_truncation <- factor(ifelse(d$TP53_class %in% c("Nontruncating", "Truncating"), d$TP53_class, NA), levels = c("Nontruncating", "Truncating"))
  d$KRAS_TP53 <- factor(ifelse(d$KRAS + d$TP53 == 0, "Neither", ifelse(d$KRAS + d$TP53 == 1, "One gene", "Both genes")), levels = c("Neither", "One gene", "Both genes"))
  d$G12D_vs_other <- factor(ifelse(d$KRAS == 0, NA, ifelse(d$KRAS_MAF == "G12D", "G12D", "Other KRAS")), levels = c("Other KRAS", "G12D"))
  d$Age10 <- d$Age / 10
  for (col in c("Sex", "T_stage", "N_stage", "Neoadjuvant", "LVI", "PNI", "R_status", "Differentiation", "M_stage")) d[[col]] <- factor(d[[col]])
  for (col in c("LVI", "PNI")) d[[col]] <- relevel(d[[col]], ref = "Negative")
  d$Differentiation <- factor(d$Differentiation, levels = c("WD", "MD", "PD"))
  d
}
all_data <- annotate_genomics(all_clinical, variants)
d <- all_data[all_data$Tumor_Sample_Barcode %in% analysis_data$clinical$Tumor_Sample_Barcode, ]
save_table(d, "18_Primary_analysis_dataset_private.tsv")

# 2. Cohort accounting, missingness and extended clinical Table 1.
flow <- data.frame(step = c("Curated workbook samples", "Unique patient identifiers", "Duplicate-ID rows excluded", "M1 patients INCLUDED after duplicates", "Primary M0+M1 unique patients", "Valid recorded recurrence follow-up"),
  n = c(nrow(all_data), n_distinct(all_data$patient_id), sum(all_data$patient_duplicate),
    sum(!all_data$patient_duplicate & all_data$M_stage == "1"), nrow(d), sum(!is.na(d$RFS_m) & !is.na(d$Recur))))
save_table(flow, "18_Cohort_flow.tsv")
flow_plot <- flow[c(1, 3, 4, 5, 6), ]; flow_plot$step <- factor(flow_plot$step, levels = rev(flow_plot$step))
save_plot(ggplot(flow_plot, aes(n, step)) + geom_col(fill = palette[1], width = .65) +
  geom_text(aes(label = n), hjust = -.15) + scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Updated PDAC cohort accounting", x = "Records / patients", y = NULL,
    caption = "Excluded counts are shown separately; these bars are not additive."), "18_Cohort_flow.png", 10, 5)
continuous <- setdiff(curated_continuous(), c("OS_months", "DFS_months"))
categorical <- curated_categorical()
describe <- function(x, cohort) {
  bind_rows(lapply(c(continuous, categorical), function(v) {
    value <- x[[v]]; observed <- sum(!is.na(value))
    if (v %in% continuous) {
      q <- if (observed) quantile(value, c(.25, .5, .75), na.rm = TRUE) else rep(NA_real_, 3)
      return(data.frame(cohort, variable = v, category = "Continuous", total_n = nrow(x), observed_n = observed,
        missing_n = sum(is.na(value)), n = NA_integer_, percent = NA_real_, median = q[2], q1 = q[1], q3 = q[3],
        mean = if (observed) mean(value, na.rm = TRUE) else NA_real_, sd = sd(value, na.rm = TRUE)))
    }
    t <- table(ifelse(is.na(value), "Missing / unavailable", as.character(value)))
    data.frame(cohort, variable = v, category = names(t), total_n = nrow(x), observed_n = observed, missing_n = sum(is.na(value)),
      n = as.integer(t), percent = as.integer(t) / nrow(x) * 100, median = NA_real_, q1 = NA_real_, q3 = NA_real_, mean = NA_real_, sd = NA_real_)
  }))
}
baseline <- bind_rows(describe(all_data, "All source samples (not unique patients)"), describe(d, "Primary M0+M1 unique patients"),
  describe(d[d$Neoadjuvant == "n", ], "Primary: no neoadjuvant"), describe(d[d$Neoadjuvant == "y", ], "Primary: neoadjuvant"))
save_table(baseline, "18_Extended_clinical_Table1.tsv")
display_table <- baseline |>
  mutate(characteristic = ifelse(category == "Continuous", variable, paste(variable, category, sep = ": ")),
    value = ifelse(category == "Continuous", sprintf("%.2f [%.2f, %.2f]; missing %d", median, q1, q3, missing_n),
      sprintf("%d (%.1f%%)", n, percent))) |>
  select(characteristic, cohort, value) |>
  tidyr::pivot_wider(names_from = cohort, values_from = value)
save_table(display_table, "18_Clinical_Table1_display.tsv")

# 3. Pathology vs four gene calls, stratified by neoadjuvant treatment (CMH).
# Generalized CMH handles multicategory T/N/grade without inventing cutoffs.
cmh_results <- cmh_counts <- cmh_pairs <- list()
for (g in genes) for (v in c("T_stage", "N_stage", "Differentiation", "LVI", "PNI", "R_status")) {
  z <- d[complete.cases(d[, c(g, v, "Neoadjuvant")]), ]
  pw <- pairwise_categorical(z[[g]], z[[v]], z$Neoadjuvant)
  if (nrow(pw)) cmh_pairs[[paste(g, v)]] <- cbind(data.frame(gene = g, characteristic = v), pw)
  arr <- table(droplevels(factor(z[[v]])), factor(z[[g]], levels = 0:1), droplevels(z$Neoadjuvant))
  test <- tryCatch(mantelhaen.test(arr, correct = FALSE), error = function(e) e)
  note <- if (inherits(test, "error")) conditionMessage(test) else "Generalized CMH; strata = neoadjuvant Yes/No"
  cmh_results[[paste(g, v)]] <- data.frame(gene = g, characteristic = v, n = nrow(z), missing_n = nrow(d) - nrow(z),
    p = if (inherits(test, "error")) NA_real_ else test$p.value, method = note)
  cmh_counts[[paste(g, v)]] <- z |> transmute(gene = g, characteristic = v, category = as.character(.data[[v]]),
    gene_status = ifelse(.data[[g]] == 1, "Retained variant", "No retained call")) |>
    count(gene, characteristic, category, gene_status, name = "n") |>
    group_by(gene, characteristic, gene_status) |> mutate(denominator = sum(n), percent = 100 * n / denominator) |> ungroup()
}
cmh <- bind_rows(cmh_results); cmh$q_BH <- p.adjust(cmh$p, "BH")
save_table(cmh, "19_Neoadjuvant_adjusted_pathology_CMH.tsv")
save_table(bind_rows(cmh_pairs), "19_Pathology_pairwise_CMH.tsv")
save_table(bind_rows(cmh_counts), "19_Pathology_by_driver_counts.tsv")
format_p <- function(x) ifelse(is.na(x), "NA", ifelse(x < .001, "<0.001", sprintf("%.3f", x)))
save_plot(ggplot(cmh, aes(gene, figure_label(characteristic), fill = -log10(pmax(q_BH, 1e-16)))) + geom_tile(color = "white") +
  geom_text(aes(label = paste0("p ", format_p(p), "\nq ", format_p(q_BH))), size = 3.2) +
  scale_fill_gradient(low = "#F3F6F9", high = "#73A6C7", na.value = "grey90") +
  labs(title = "Driver calls and pathology, adjusted for neoadjuvant treatment", x = NULL, y = NULL,
    fill = "-log10(BH q)", caption = "CMH omnibus tests. BH correction across the 24 gene-characteristic tests."), "19_Adjusted_pathology_associations.png", 10, 6)

# 4. Cox models. Prespecified covariates; no p-value-based variable selection.
cox_rows <- ph_rows <- model_rows <- cox_pairs <- list()
fit_cox <- function(data, exposure, adjust, family) {
  adjust <- setdiff(adjust, exposure)
  terms <- c(exposure, adjust)
  z <- droplevels(data[complete.cases(data[, c("OS_m", "survive", terms)]), ])
  id <- paste(family, exposure, sep = ": ")
  form <- reformulate(terms, response = "survival::Surv(OS_m, survive)")
  warn <- character()
  fit <- tryCatch(withCallingHandlers(survival::coxph(form, data = z, x = TRUE, model = TRUE, ties = "efron"),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
  if (inherits(fit, "error")) {
    model_rows[[id]] <<- data.frame(model = id, n = nrow(z), events = sum(z$survive), parameters = NA, events_per_parameter = NA, status = conditionMessage(fit), formula = deparse(form))
    return(invisible(NULL))
  }
  pw <- pairwise_cox(fit, z, exposure)
  if (nrow(pw)) cox_pairs[[id]] <<- cbind(data.frame(model = id, exposure, family), pw)
  s <- summary(fit); pcount <- length(coef(fit))
  invalid <- length(warn) || any(!is.finite(coef(fit))) || any(!is.finite(s$conf.int))
  epv <- sum(z$survive) / pcount
  status <- if (invalid) paste("Unstable;", paste(warn, collapse = "; ")) else if (epv < 10) "Sparse model: <10 events/parameter" else "Estimated"
  model_rows[[id]] <<- data.frame(model = id, n = nrow(z), events = sum(z$survive), parameters = pcount, events_per_parameter = epv, status, formula = deparse(form))
  rows <- data.frame(model = id, family, exposure, term = rownames(s$coefficients), n = nrow(z), events = sum(z$survive),
    HR = s$conf.int[, "exp(coef)"], lower95 = s$conf.int[, "lower .95"], upper95 = s$conf.int[, "upper .95"],
    p = s$coefficients[, "Pr(>|z|)"], status, row.names = NULL)
  # Determine exposure terms from the design-matrix term assignment, not prefixes.
  assignment <- fit$assign[[exposure]]
  rows$is_exposure <- seq_len(nrow(rows)) %in% assignment
  cox_rows[[id]] <<- rows
  ph <- tryCatch(survival::cox.zph(fit), error = identity)
  if (!inherits(ph, "error")) ph_rows[[id]] <<- data.frame(model = id, term = rownames(ph$table), chisq = ph$table[, 1], df = ph$table[, 2], p = ph$table[, 3])
  invisible(fit)
}
basic <- c("Age10", "Sex", "T_stage", "Neoadjuvant", "M_stage")
extended <- c(basic, "N_stage", "LVI", "PNI", "R_status")
exposures <- c(genes, "KRAS_MAF", "driver_group", "Any_driver", "TP53_truncation", "KRAS_TP53", "G12D_vs_other")
for (v in exposures) {
  fit_cox(d, v, character(), "Univariable genomic")
  fit_cox(d, v, basic, "Adjusted genomic")
  fit_cox(d, v, extended, "Extended genomic sensitivity")
}
for (v in c("Age10", "Sex", "Neoadjuvant", "T_stage", "N_stage", "LVI", "PNI", "R_status", "Differentiation", "Operation", "Diabetes")) fit_cox(d, v, character(), "Univariable clinical")
fit_cox(d, "N_stage", c("Age10", "Sex", "T_stage", "Neoadjuvant", "M_stage", "LVI", "PNI", "R_status", "Differentiation"), "Joint clinical model")
# Localized-only sensitivity complements the primary M0+M1 analysis.
for (v in genes) fit_cox(d[d$M_stage == "0", ], v, setdiff(basic, "M_stage"), "M0-only sensitivity")
concordant <- d[as.character(d$KRAS_MAF) == d$KRAS_report_group, ]
for (v in c("KRAS", "KRAS_MAF")) fit_cox(concordant, v, basic, "KRAS concordant-only sensitivity")
cox <- bind_rows(cox_rows)
ph_table <- bind_rows(ph_rows)
cox$PH_exposure_p <- vapply(seq_len(nrow(cox)), function(i) {
  x <- ph_table$p[ph_table$model == cox$model[i] & ph_table$term == cox$exposure[i]]
  if (length(x)) x[1] else NA_real_
}, numeric(1))
cox$PH_global_p <- vapply(cox$model, function(id) {
  x <- ph_table$p[ph_table$model == id & ph_table$term == "GLOBAL"]
  if (length(x)) x[1] else NA_real_
}, numeric(1))
cox$q_BH <- NA_real_
for (f in unique(cox$family)) {
  idx <- which(cox$family == f & cox$is_exposure)
  cox$q_BH[idx] <- p.adjust(cox$p[idx], "BH")
}
save_table(cox, "20_OS_Cox_coefficients.tsv")
save_table(bind_rows(model_rows), "20_OS_Cox_model_diagnostics.tsv")
save_table(ph_table, "20_OS_Cox_PH_tests.tsv")
forest <- cox[cox$family == "Adjusted genomic" & cox$is_exposure & cox$status == "Estimated", ]
forest$label <- vapply(seq_len(nrow(forest)), function(i) {
  e <- forest$exposure[i]; term <- forest$term[i]
  label <- if (e %in% genes) paste0(e, ": variant vs no call") else
    if (e == "KRAS_MAF") paste0(sub("KRAS_MAF", "KRAS ", term), " vs no call") else
    if (e == "driver_group") paste0("Driver genes: ", sub("driver_group", "", term), " vs 1") else
    if (e == "Any_driver") "Any driver call vs none" else
    if (e == "TP53_truncation") "TP53: truncating vs nontruncating" else
    if (e == "KRAS_TP53") paste0("KRAS/TP53: ", sub("KRAS_TP53", "", term), " vs neither") else
    "KRAS: G12D vs other variants"
  paste0(label, if (!is.na(forest$PH_exposure_p[i]) && forest$PH_exposure_p[i] < .05) " *" else "", " (n=", forest$n[i], ")")
}, character(1))
forest$label <- factor(forest$label, levels = rev(unique(forest$label)))
save_plot(ggplot(forest, aes(HR, label, xmin = lower95, xmax = upper95)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") + geom_errorbar(orientation = "y", width = .15) +
  geom_point(color = palette[1], size = 2.5) + scale_x_log10() +
  labs(title = "Overall survival: adjusted genomic associations", x = "Hazard ratio (95% CI; log scale)", y = NULL,
    caption = "Separate models adjusted for age (per 10 years), sex, T category, neoadjuvant treatment and M stage.\n* Exposure proportional-hazards test p<0.05: a constant HR may be inadequate; see 36-month RMST sensitivity.\nRetained MAF calls do not establish pathogenicity or biological wild type. Exploratory, not causal."),
  "20_Adjusted_OS_forest.png", 12, max(7, .35 * nrow(forest) + 2))

# 5. Kaplan-Meier estimates: time zero, censor marks, confidence intervals,
# shared time axis, and number at risk. Endpoint is OS from surgery.
save_table(bind_rows(cox_pairs), "20_Cox_pairwise_contrasts.tsv")
km_summaries <- risk_tables <- logrank_rows <- km_pairs <- list()
km_plot <- function(data, variable, title, file) {
  z <- data[complete.cases(data[, c(variable, "OS_m", "survive")]) & data$OS_m > 0, ]
  z$group <- droplevels(factor(z[[variable]]))
  if (variable %in% c(genes, "Any_driver")) z$group <- factor(z[[variable]], levels = 0:1, labels = c("No retained call", "Retained variant"))
  if (nlevels(z$group) < 2) return(invisible(NULL))
  fit <- survival::survfit(survival::Surv(OS_m, survive) ~ group, data = z, conf.type = "log-log")
  test <- survival::survdiff(survival::Surv(OS_m, survive) ~ group, data = z)
  pv <- pchisq(test$chisq, length(test$n) - 1, lower.tail = FALSE)
  pw <- pairwise_logrank(z$OS_m, z$survive, z$group)
  km_pairs[[variable]] <<- cbind(data.frame(variable), pw)
  logrank_rows[[variable]] <<- data.frame(variable, n = nrow(z), p = pv)
  ss <- summary(fit, censored = TRUE)
  curve <- data.frame(time = ss$time, survival = ss$surv, lower = ss$lower, upper = ss$upper, censor = ss$n.censor,
    group = sub("^group=", "", as.character(ss$strata))) |>
    bind_rows(data.frame(time = 0, survival = 1, lower = 1, upper = 1, censor = 0, group = levels(z$group))) |>
    arrange(group, time)
  curve$group <- factor(curve$group, levels = levels(z$group))
  times <- seq(0, floor(max(z$OS_m) / 12) * 12, 12)
  risk <- tidyr::expand_grid(group = levels(z$group), time = times) |>
    rowwise() |> mutate(n_risk = sum(as.character(z$group) == group & z$OS_m >= time)) |> ungroup()
  risk_tables[[variable]] <<- mutate(risk, variable = variable)
  med <- as.data.frame(summary(fit)$table)
  med$group <- sub("^group=", "", rownames(med))
  for (g in levels(z$group)) {
    gz <- z[z$group == g, ]; gf <- survival::survfit(survival::Surv(OS_m, survive) ~ 1, data = gz, conf.type = "log-log")
    for (h in c(12, 36, 60)) {
      hs <- if (max(gz$OS_m) >= h) summary(gf, times = h, extend = FALSE) else NULL
      mr <- med[med$group == g, ]
      km_summaries[[paste(variable, g, h)]] <<- data.frame(variable, group = g, n = nrow(gz), events = sum(gz$survive),
        median_months = mr$median, median_lower95 = mr$`0.95LCL`, median_upper95 = mr$`0.95UCL`, time_months = h,
        survival = if (!is.null(hs) && length(hs$surv)) hs$surv else NA_real_,
        lower95 = if (!is.null(hs) && length(hs$lower)) hs$lower else NA_real_, upper95 = if (!is.null(hs) && length(hs$upper)) hs$upper else NA_real_,
        n_risk = sum(gz$OS_m >= h))
    }
  }
  cols <- setNames(palette[seq_len(nlevels(z$group))], levels(z$group))
  p <- ggplot(curve, aes(time, survival, color = group, fill = group)) +
    geom_step(linewidth = .8) + geom_point(data = curve[curve$censor > 0, ], shape = 3, size = 1.2) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = times, limits = c(0, max(z$OS_m))) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title = title, subtitle = sprintf("M0+M1 unique patients: N=%d; log-rank p=%.3g", nrow(z), pv),
      x = "Months from surgery", y = "Overall survival", color = NULL, fill = NULL,
      caption = pairwise_caption(pw))
  if (nlevels(z$group) <= 2) p <- p + geom_ribbon(data = curve[is.finite(curve$lower) & is.finite(curve$upper), ],
    aes(ymin = lower, ymax = upper), alpha = .10, color = NA)
  rp <- ggplot(risk, aes(time, factor(group, levels = rev(levels(z$group))), label = n_risk)) + geom_text(size = 3.5) +
    scale_x_continuous(breaks = times, limits = c(0, max(z$OS_m))) +
    labs(title = "Number at risk", x = "Months from surgery", y = NULL) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  grDevices::png(output_path(config, file), width = 11, height = 8, units = "in", res = min(config$plot_dpi, 400L))
  g1 <- ggplotGrob(p); g2 <- ggplotGrob(rp)
  widths <- grid::unit.pmax(g1$widths, g2$widths); g1$widths <- widths; g2$widths <- widths
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = c(.73, .27))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1)); grid::grid.draw(g1); grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2)); grid::grid.draw(g2); grid::upViewport(2)
  grDevices::dev.off()
}
for (g in genes) km_plot(d, g, paste("OS by", g, "retained variant status"), paste0("21_OS_", g, ".png"))
km_plot(d, "driver_group", "OS by number of genes with retained variants", "21_OS_driver_count_0_1_2_3plus.png")
km_plot(d, "TP53_truncation", "OS by TP53 truncating variant class", "21_OS_TP53_truncation.png")
km_plot(d, "KRAS_TP53", "OS by KRAS and TP53 co-occurrence", "21_OS_KRAS_TP53_cooccurrence.png")
km_plot(d, "G12D_vs_other", "OS: KRAS G12D vs other KRAS calls", "21_OS_KRAS_G12D_vs_other.png")
save_table(bind_rows(km_pairs), "21_OS_pairwise_logrank.tsv")
save_table(bind_rows(km_summaries), "21_OS_medians_and_landmark_estimates.tsv")
save_table(bind_rows(risk_tables), "21_OS_number_at_risk.tsv")
lr <- bind_rows(logrank_rows); lr$q_BH <- p.adjust(lr$p, "BH"); save_table(lr, "21_OS_logrank_tests.tsv")

# Restricted mean survival to 36 months does not require proportional hazards.
# These are unadjusted contrasts, not replacements for the adjusted Cox models.
rmst_rows <- list(); rmst_pairs <- list()
for (v in c(genes, "Any_driver", "G12D_vs_other", "driver_group")) {
  z <- d[complete.cases(d[, c(v, "OS_m", "survive")]), ]; z$group <- droplevels(factor(z[[v]]))
  if (any(tapply(z$OS_m, z$group, max) < 36)) next
  f <- survival::survfit(survival::Surv(OS_m, survive) ~ group, data = z)
  m <- as.data.frame(summary(f, rmean = 36)$table)
  pairs <- t(combn(seq_len(nrow(m)), 2))
  pp <- bind_rows(lapply(seq_len(nrow(pairs)), function(k) {
    a <- pairs[k, 1]; b <- pairs[k, 2]; delta <- m$rmean[a] - m$rmean[b]
    se <- sqrt(m[["se(rmean)"]][a]^2 + m[["se(rmean)"]][b]^2)
    data.frame(variable = v, group1 = sub("^group=", "", rownames(m)[a]), group2 = sub("^group=", "", rownames(m)[b]),
      n1 = m$records[a], n2 = m$records[b], horizon_months = 36, difference_months = delta,
      lower95 = delta - qnorm(.975) * se, upper95 = delta + qnorm(.975) * se,
      statistic = delta / se, p = 2 * pnorm(-abs(delta / se)))
  }))
  pp$p_holm <- p.adjust(pp$p, "holm"); rmst_pairs[[v]] <- pp
  ref <- 1L
  for (i in seq_len(nrow(m))) {
    delta <- m$rmean[i] - m$rmean[ref]
    se <- if (i == ref) NA_real_ else sqrt(m$`se(rmean)`[i]^2 + m$`se(rmean)`[ref]^2)
    rmst_rows[[paste(v, i)]] <- data.frame(variable = v, group = sub("^group=", "", rownames(m)[i]),
      reference = sub("^group=", "", rownames(m)[ref]), horizon_months = 36, n = m$records[i],
      rmst_months = m$rmean[i], rmst_se = m$`se(rmean)`[i], difference_months = delta,
      difference_lower95 = delta - qnorm(.975) * se, difference_upper95 = delta + qnorm(.975) * se,
      p = if (is.finite(se) && se > 0) 2 * pnorm(-abs(delta / se)) else NA_real_)
  }
}
rmst <- bind_rows(rmst_rows); rmst$q_BH <- p.adjust(rmst$p, "BH")
save_table(rmst, "21_Unadjusted_RMST_36months.tsv")
save_table(bind_rows(rmst_pairs), "21_RMST_pairwise_contrasts.tsv")

# 6. Observed recurrence patterns (NOT cumulative incidence or Fine-Gray).
rec <- d[d$Recur == 1 & !is.na(d$Recur), ]
recurrence <- rec |> count(Recurrence_pattern, name = "n") |> mutate(denominator_recurrent = nrow(rec), percent = n / nrow(rec) * 100)
distant <- rec[!is.na(rec$Recurrence_pattern) & rec$Recurrence_pattern == "Distant only", ] |>
  count(Distant_pattern, name = "n") |> mutate(denominator_distant_only = sum(n), percent = n / sum(n) * 100)
save_table(recurrence, "22_Observed_recurrence_patterns.tsv")
save_table(distant, "22_Distant_only_recurrence_sites.tsv")
save_plot(ggplot(recurrence, aes(Recurrence_pattern, n, fill = Recurrence_pattern)) + geom_col(width = .65) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, percent)), vjust = -.4) +
  scale_fill_manual(values = palette) + scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Observed recurrence distribution", x = NULL, y = "Patients with documented recurrence",
    caption = "Percentages among patients with defined recorded recurrence; not a time-to-event cumulative incidence.") + theme(legend.position = "none"), "22_Observed_recurrence_patterns.png", 9, 6)
save_plot(ggplot(distant, aes(reorder(Distant_pattern, n), n)) + geom_col(fill = palette[1]) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, percent)), hjust = -.1) + coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  labs(title = "Sites among distant-only recurrence", x = NULL, y = "Patients",
    caption = "Excludes local + distant recurrence; the source site variable applies only to distant-only cases."), "22_Distant_only_recurrence_sites.png", 9, 6)
by_gene <- bind_rows(lapply(genes, function(g) rec |> mutate(gene = g, status = ifelse(.data[[g]] == 1, "Retained variant", "No retained call")) |>
  count(gene, status, Recurrence_pattern, name = "n") |> group_by(gene, status) |> mutate(denominator_recurrent = sum(n), percent = n / sum(n) * 100) |> ungroup()))
save_table(by_gene, "22_Recurrence_patterns_by_driver.tsv")

# 7. Cross-check report KRAS against MAF groups; actual TMB is kept separate.
concordance <- d |> transmute(Tumor_Sample_Barcode, patient_id, KRAS_report, KRAS_report_group, KRAS_MAF,
  comparable_group_agreement = KRAS_report_group == as.character(KRAS_MAF))
save_table(concordance, "23_KRAS_report_MAF_concordance_private.tsv")
tab <- concordance |> count(KRAS_report_group, KRAS_MAF, name = "n")
tab <- tidyr::complete(tab, KRAS_report_group = levels(d$KRAS_MAF), KRAS_MAF, fill = list(n = 0L))
tab$KRAS_report_group <- factor(tab$KRAS_report_group, levels = rev(levels(all_data$KRAS_MAF)))
save_table(tab, "23_KRAS_report_MAF_concordance_counts.tsv")
save_plot(ggplot(tab, aes(KRAS_MAF, KRAS_report_group, fill = n)) + geom_tile(color = "white") + geom_text(aes(label = n)) +
  scale_fill_gradient(low = "#EEF4F7", high = "#63A2C9") +
  labs(title = "KRAS: report-derived vs MAF-derived categories", x = "MAF category", y = "Report category",
    caption = "All unique patients. Other KRAS collapses rare/multiple calls; agreement does not validate exact rare subtypes.\nReport WT/Not detected is mapped only for this comparison; neither label proves biological wild type."), "23_KRAS_concordance.png", 10, 7)
counts <- variants |> count(Tumor_Sample_Barcode, name = "retained_variant_count")
tmb <- left_join(d, counts, by = "Tumor_Sample_Barcode")
rho <- cor(tmb$TMB_report, tmb$retained_variant_count, method = "spearman", use = "complete.obs")
save_table(tmb[, c("Tumor_Sample_Barcode", "TMB_report", "retained_variant_count", "MSI")], "23_TMB_and_variant_counts_private.tsv")
save_plot(ggplot(tmb[!is.na(tmb$TMB_report), ], aes(retained_variant_count, TMB_report)) + geom_point(alpha = .4, color = palette[1]) +
  labs(title = "Reported TMB and Variant count", subtitle = sprintf("Complete pairs: %d; Spearman rho=%.3f", sum(!is.na(tmb$TMB_report)), rho),
    x = "Variant count", y = "Reported TMB (mutations/Mb)",
    caption = "Distinct measures: no panel-size normalization is inferred from the MAF row count."), "23_Reported_TMB_vs_variant_count.png", 9, 6)

# 8. DNA-repair gene calls (exploratory; not germline HRD or a functional HRD test).
hrr <- repair_gene_sets()$HRD
mmr <- repair_gene_sets()$MMR
repair <- data.frame(gene = c(hrr, mmr), pathway = c(rep("HRD-related gene list", length(hrr)), rep("MMR gene list", length(mmr))))
repair$n <- vapply(repair$gene, function(g) n_distinct(variants$Tumor_Sample_Barcode[variants$Hugo_Symbol == g & variants$Tumor_Sample_Barcode %in% d$Tumor_Sample_Barcode]), integer(1))
repair$denominator <- nrow(d); repair$percent <- repair$n / nrow(d) * 100
save_table(repair, "24_DNA_repair_gene_variant_counts.tsv")
save_plot(ggplot(repair, aes(reorder(gene, n), percent, fill = pathway)) + geom_col() + coord_flip() +
  scale_fill_manual(values = palette[c(1, 2)]) + labs(title = "Patients with variants in HRD-related and MMR genes", x = NULL, y = "% of primary cohort", fill = NULL,
    caption = "McIntyre 2020 gene lists. Retained MAF variants, without pathogenicity or panel-coverage adjudication.\nZero means no retained call in this file; not proven assay-negative. Not a germline/functional HRD result."), "24_DNA_repair_gene_variants.png", 10, 8)
save_table(data.frame(item = c("Primary analysis cohort", "Baseline table source rows", "Recurrence-only endpoint", "Unavailable paper analyses", "Cox adjustment departures", "NGS ascertainment", "TMB", "Multiple testing"),
  detail = c(paste("M0+M1 and unique patient IDs; N=", nrow(d)), paste("All source samples N=", nrow(all_data), "; duplicate samples retained only in descriptive source table"),
    "User confirmed recorded recurrence event. Legacy RFS filenames are retained, but death-inclusive DFS/RFS is not estimated.",
    "No LOH, copy-number loss, germline variants, pathogenicity classification, or validated competing-risk recurrence times.",
    "Age/10, sex, T category, neoadjuvant y/n and M stage; no Clavien-Dindo, anatomic resectability or detailed therapy modalities. Not an exact paper replication.",
    "Sequencing/entry dates are absent. Survivor/ascertainment bias cannot be corrected by delayed entry.",
    "Reported mutations/Mb is distinct from legacy retained MAF variant count. Tumor-size and assay units remain unconfirmed.",
    "BH correction for omnibus families; Holm for all pairs within each variable/model/endpoint. Pointwise CIs are not multiplicity-adjusted. Exploratory study.")), "24_Analysis_definitions_and_limits.tsv")
save_table(data.frame(line = capture.output(sessionInfo())), "24_R_session_info.tsv")

message("Literature extensions complete: ", config$output_dir)
