# Additional manuscript-facing analyses requested for the curated clinical update.
# This module never exports patient identifiers. It is called after the core
# genomic and survival models have been fitted in scripts/manuscript_analysis.R.

run_manuscript_story_extensions <- function(d, out, savefig, tab, pal, extended,
                                            heat, fg, kk, kf, core, cpair) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("Package gridExtra is required.")

  fmt_p <- function(x) ifelse(is.na(x), "NE", ifelse(x < .001, "<0.001", sprintf("%.3f", x)))
  clean_label <- function(x) {
    labs <- c(Age="Age", Sex="Sex", Neoadjuvant="Neoadjuvant treatment",
      Preop_platinum_exposure="Preoperative platinum", Adjuvant_treatment="Adjuvant treatment",
      T_stage="T category", N_stage="N category (0/1/2)",
      N_category="N category (0/1)", M_stage="M category", Differentiation="Differentiation",
      LVI="LVI", PNI="PNI", R_status="Margin", KRAS_group="KRAS subtype",
      KRAS="KRAS variant", TP53="TP53 variant", SMAD4="SMAD4 variant",
      CDKN2A="CDKN2A variant", TMB="Reported TMB", Variant_count="Variant count",
      Tumor_size="Tumor size", CA19_9="CA19-9")
    unname(ifelse(x %in% names(labs), labs[x], x))
  }

  save_table_png <- function(df, title, filename, widths = NULL, font_size = 8.5,
                             width = 16, height = NULL) {
    shown <- as.data.frame(df, stringsAsFactors = FALSE)
    stripe <- rep(c("#FFFFFF", "#F4F4F4"), length.out = nrow(shown))
    theme <- gridExtra::ttheme_minimal(
      base_size = font_size, base_family = "Pretendard",
      core = list(bg_params = list(fill = stripe, col = NA),
        fg_params = list(hjust = 0, x = .03)),
      colhead = list(bg_params = list(fill = "#E6E6E6", col = NA),
        fg_params = list(fontface = "bold", hjust = 0, x = .03)))
    grob <- gridExtra::tableGrob(shown, rows = NULL, theme = theme)
    if (!is.null(widths)) grob$widths <- grid::unit(widths, "null")
    if (is.null(height)) height <- max(3.5,
      grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE) + .75)
    canvas <- gridExtra::arrangeGrob(grob, top = grid::textGrob(title,
      gp = grid::gpar(fontfamily = "Pretendard", fontsize = 16, fontface = "bold"),
      x = .01, hjust = 0))
    ragg::agg_png(file.path(out, "figures", paste0(filename, ".png")),
      width = width, height = height, units = "in", res = 300, background = "white")
    grid::grid.draw(canvas); grDevices::dev.off()
  }

  omnibus_categorical <- function(outcome, group) {
    ok <- complete.cases(outcome, group)
    a <- table(droplevels(factor(outcome[ok])), droplevels(factor(group[ok])))
    if (nrow(a) < 2 || ncol(a) < 2) return(list(p = NA_real_, method = "Not estimable", effect = NA_real_))
    chi <- suppressWarnings(chisq.test(a, correct = FALSE))
    sparse <- any(chi$expected < 5)
    fit <- if (sparse) fisher.test(a, simulate.p.value = any(dim(a) > 2), B = 20000) else chi
    v <- sqrt(unname(chi$statistic) / (sum(a) * min(nrow(a) - 1, ncol(a) - 1)))
    list(p = fit$p.value,
      method = if (sparse) "Fisher exact; Monte Carlo 20,000 for tables larger than 2x2" else "Pearson chi-square",
      effect = v)
  }

  pathology_status <- function(gene, variable) {
    z <- d[complete.cases(d[, c(gene, variable)]), c(gene, variable), drop = FALSE]
    z$status <- factor(ifelse(z[[gene]] == 1, "Variant detected", "Not detected"),
      levels = c("Not detected", "Variant detected"))
    z$category <- droplevels(factor(z[[variable]]))
    global <- omnibus_categorical(z$status, z$category)
    counts <- z |> dplyr::count(category, status, name = "n") |>
      dplyr::group_by(category) |> dplyr::mutate(denominator = sum(n), percent = 100 * n / denominator) |>
      dplyr::ungroup()
    pairs <- pairwise_categorical(z$status, z$category)
    pairs$gene <- gene; pairs$characteristic <- variable
    annotation <- paste0("Overall ", global$method, " p ", fmt_p(global$p),
      "; Cramer's V ", sprintf("%.2f", global$effect), ".\n", pairwise_caption(pairs, 85))
    plot <- ggplot2::ggplot(counts, ggplot2::aes(category, percent, fill = status)) +
      ggplot2::geom_col(position = "stack", width = .68) +
      ggplot2::geom_text(ggplot2::aes(label = paste0(n, "/", denominator)),
        position = ggplot2::position_stack(vjust = .5), size = 3.5) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
      ggplot2::scale_fill_manual(values = c("Not detected" = "#D4D4D4", "Variant detected" = pal[2])) +
      ggplot2::labs(x = clean_label(variable), y = "Patients within category (%)", fill = NULL,
        title = paste(gene, "variant detection by", tolower(clean_label(variable))), caption = annotation) +
      ggplot2::theme(legend.position = "top")
    list(plot = plot, counts = counts, pairs = pairs,
      global = data.frame(gene = gene, characteristic = variable, n = nrow(z),
        p = global$p, method = global$method, cramers_v = global$effect))
  }

  # Figure 2: make the comparison direction explicit. Percentages are variant
  # detection within pathology categories, not pathology composition within a
  # variant-status group.
  tp_grade <- pathology_status("TP53", "Differentiation")
  cd_nstage <- pathology_status("CDKN2A", "N_stage")
  cd_ncat <- pathology_status("CDKN2A", "N_category")
  tab(dplyr::bind_rows(tp_grade$global, cd_nstage$global, cd_ncat$global),
    "Input_Figure2_global_tests.tsv")
  tab(dplyr::bind_rows(tp_grade$pairs, cd_nstage$pairs, cd_ncat$pairs),
    "Input_Figure2_pairwise_tests.tsv")
  tab(dplyr::bind_rows(
    dplyr::mutate(tp_grade$counts, gene = "TP53", characteristic = "Differentiation"),
    dplyr::mutate(cd_nstage$counts, gene = "CDKN2A", characteristic = "N_stage"),
    dplyr::mutate(cd_ncat$counts, gene = "CDKN2A", characteristic = "N_category")),
    "Input_Figure2_counts.tsv")
  savefig(tp_grade$plot, "Figure2_TP53_differentiation", 7.2, 6.2)
  savefig(cd_nstage$plot, "Figure2_CDKN2A_N_stage", 7.2, 6.2)
  savefig(cd_ncat$plot, "Figure2_CDKN2A_N_category", 7.2, 6.2)
  fig2 <- cowplot::plot_grid(heat, tp_grade$plot, cd_nstage$plot, cd_ncat$plot,
    ncol = 2, rel_heights = c(.95, 1.05))
  savefig(fig2, "Figure2_Genes_and_pathology_revised", 15, 12)

  compare_groups <- function(z, group_var, variables, family_label) {
    rows <- pairs <- list()
    for (nm in variables) {
      ok <- complete.cases(z[, c(group_var, nm)])
      zz <- z[ok, , drop = FALSE]
      g <- droplevels(factor(zz[[group_var]]))
      x <- zz[[nm]]
      if (length(unique(g)) < 2) next
      is_continuous <- is.numeric(x) && length(unique(x[!is.na(x)])) > 5
      if (is_continuous) {
        global <- if (nlevels(g) == 2) suppressWarnings(wilcox.test(x ~ g, exact = FALSE)) else kruskal.test(x, g)
        eff <- if (nlevels(g) == 2) {
          w <- unname(global$statistic); n1 <- sum(g == levels(g)[1]); n2 <- sum(g == levels(g)[2])
          2 * w / (n1 * n2) - 1
        } else NA_real_
        rows[[nm]] <- data.frame(family = family_label, variable = nm, type = "Continuous",
          n = nrow(zz), p = global$p.value, method = global$method,
          effect = eff, effect_definition = if (nlevels(g) == 2) paste0("Rank-biserial correlation; positive values indicate higher values in ", levels(g)[1]) else "Not defined for omnibus Kruskal-Wallis")
        pw <- pairwise_numeric(x, g)
      } else {
        global <- omnibus_categorical(x, g)
        rows[[nm]] <- data.frame(family = family_label, variable = nm, type = "Categorical",
          n = nrow(zz), p = global$p, method = global$method, effect = global$effect,
          effect_definition = "Cramer's V; magnitude only")
        pw <- pairwise_categorical(x, g)
      }
      if (nrow(pw)) { pw$family <- family_label; pw$variable <- nm; pairs[[nm]] <- pw }
    }
    result <- dplyr::bind_rows(rows)
    result$p_holm <- p.adjust(result$p, "holm")
    list(global = result, pairwise = dplyr::bind_rows(pairs))
  }

  describe_by_group <- function(z, group_var, variables) {
    dplyr::bind_rows(lapply(variables, function(nm) {
      x <- z[[nm]]; g <- z[[group_var]]; ok <- complete.cases(x, g)
      x <- x[ok]; g <- droplevels(factor(g[ok]))
      if (is.numeric(x) && length(unique(x)) > 5) {
        return(dplyr::bind_rows(lapply(levels(g), function(lev) {
          y <- x[g == lev]; q <- stats::quantile(y, c(.25, .5, .75), na.rm = TRUE)
          data.frame(variable = nm, group = lev, category = "Median [IQR]", n = length(y),
            value = sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3]), percent = NA_real_)
        })))
      }
      dplyr::bind_rows(lapply(levels(g), function(lev) {
        y <- as.character(x[g == lev]); tt <- table(y)
        data.frame(variable = nm, group = lev, category = names(tt), n = as.integer(tt),
          value = sprintf("%d (%.1f%%)", as.integer(tt), 100 * as.integer(tt) / length(y)),
          percent = 100 * as.integer(tt) / length(y))
      }))
    }))
  }

  association_plot <- function(x, title, subtitle) {
    x$label <- factor(clean_label(x$variable), levels = rev(clean_label(x$variable)))
    ggplot2::ggplot(x, ggplot2::aes(-log10(pmax(p_holm, 1e-12)), label)) +
      ggplot2::geom_vline(xintercept = -log10(.05), linetype = 2, color = "#777777") +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = -log10(pmax(p_holm, 1e-12)), yend = label), color = "#BDBDBD") +
      ggplot2::geom_point(ggplot2::aes(color = type), size = 3) +
      ggplot2::geom_text(ggplot2::aes(label = paste0("p ", fmt_p(p), "; Holm p ", fmt_p(p_holm))),
        hjust = 0, nudge_x = .08, size = 3.1) +
      ggplot2::scale_color_manual(values = c(Categorical = pal[2], Continuous = pal[3])) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, .45))) +
      ggplot2::labs(x = expression(-log[10]("Holm-adjusted p")), y = NULL, color = NULL,
        title = title, subtitle = subtitle, caption = "All prespecified variables are shown; the dashed line marks adjusted p = 0.05.") +
      ggplot2::theme(legend.position = "top")
  }

  km_analysis <- function(z, group_var, time_var, event_var, title, y_title, family_name) {
    z <- z[complete.cases(z[, c(group_var, time_var, event_var)]), , drop = FALSE]
    z <- z[is.finite(z[[time_var]]) & z[[time_var]] > 0 & z[[event_var]] %in% 0:1, , drop = FALSE]
    z$group <- droplevels(factor(z[[group_var]])); z$time <- z[[time_var]]; z$event <- z[[event_var]]
    fit <- survival::survfit(survival::Surv(time, event) ~ group, data = z, conf.type = "log-log")
    lr <- survival::survdiff(survival::Surv(time, event) ~ group, data = z)
    p <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
    pw <- pairwise_logrank(z$time, z$event, z$group)
    st <- summary(fit, censored = TRUE)
    curve <- data.frame(time = st$time, prob = st$surv, censor = st$n.censor,
      group = sub("^group=", "", as.character(st$strata))) |>
      dplyr::bind_rows(data.frame(time = 0, prob = 1, censor = 0, group = levels(z$group))) |>
      dplyr::arrange(group, time)
    max_time <- max(z$time); breaks <- seq(0, floor(max_time / 12) * 12, 12)
    if (length(breaks) > 9) breaks <- pretty(c(0, max_time), n = 7)
    risk <- tidyr::expand_grid(group = levels(z$group), time = breaks) |>
      dplyr::rowwise() |> dplyr::mutate(n_risk = sum(z$group == group & z$time >= time)) |> dplyr::ungroup()
    cols <- setNames(pal[seq_len(nlevels(z$group))], levels(z$group))
    p1 <- ggplot2::ggplot(curve, ggplot2::aes(time, prob, color = group)) +
      ggplot2::geom_step(linewidth = .85) +
      ggplot2::geom_point(data = curve[curve$censor > 0, ], shape = 3, size = .7) +
      ggplot2::scale_color_manual(values = cols) + ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::scale_x_continuous(breaks = breaks) +
      ggplot2::labs(title = title, subtitle = sprintf("N=%d, events=%d; global log-rank p %s", nrow(z), sum(z$event), fmt_p(p)),
        x = NULL, y = y_title, color = NULL, caption = pairwise_caption(pw, 95)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    risk$group <- factor(risk$group, levels = rev(levels(z$group)))
    p2 <- ggplot2::ggplot(risk, ggplot2::aes(time, group, label = n_risk, color = group)) +
      ggplot2::geom_text(size = 3) + ggplot2::scale_color_manual(values = cols) +
      ggplot2::scale_x_continuous(breaks = breaks) + ggplot2::labs(x = "Months", y = NULL) +
      ggplot2::theme_classic(base_family = "Pretendard", base_size = 10) +
      ggplot2::theme(legend.position = "none", axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    tab(data.frame(family = family_name, variable = group_var, endpoint = y_title,
      n = nrow(z), events = sum(z$event), chisq = lr$chisq, df = length(lr$n) - 1, p = p,
      method = "Global log-rank"), paste0("Input_", family_name, "_logrank.tsv"))
    tab(pw, paste0("Input_", family_name, "_pairwise_logrank.tsv"))
    list(plot = cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(.77, .23), align = "v", axis = "lr"),
      global_p = p, pairwise = pw, n = nrow(z), events = sum(z$event), data = z)
  }

  cox_row <- function(z, group_var, time_var, event_var, covariates, model_name, endpoint) {
    needed <- unique(c(group_var, time_var, event_var, all.vars(reformulate(covariates))))
    z <- droplevels(z[complete.cases(z[, needed]), , drop = FALSE])
    z$time <- z[[time_var]]; z$event <- z[[event_var]]
    covariates <- covariates[vapply(covariates, function(a)
      all(vapply(all.vars(reformulate(a)), function(v) length(unique(z[[v]])) > 1, logical(1))), logical(1))]
    form <- reformulate(c(group_var, covariates), response = "survival::Surv(time,event)")
    fit <- survival::coxph(form, data = z, x = TRUE, model = TRUE, ties = "efron")
    base_form <- reformulate(covariates, response = "survival::Surv(time,event)")
    global0 <- survival::coxph(base_form, data = z, x = TRUE, model = TRUE, ties = "efron")
    lrt <- anova(global0, fit, test = "LRT")
    pairs <- if (nlevels(factor(z[[group_var]])) > 2) pairwise_cox(fit, z, group_var) else {
      s <- summary(fit); idx <- fit$assign[[group_var]]; rr <- s$coefficients[idx, , drop = FALSE]; cc <- s$conf.int[idx, , drop = FALSE]
      lev <- levels(factor(z[[group_var]]))
      data.frame(group1 = lev[2], group2 = lev[1], n1 = sum(z[[group_var]] == lev[2]), n2 = sum(z[[group_var]] == lev[1]),
        statistic = rr[1, "z"], df = 1, p = rr[1, "Pr(>|z|)"], status = "Estimated",
        HR = cc[1, "exp(coef)"], lower95 = cc[1, "lower .95"], upper95 = cc[1, "upper .95"],
        p_holm = rr[1, "Pr(>|z|)"], method = "Cox Wald Z; pointwise 95% CI", family_size = 1)
    }
    ph <- tryCatch(survival::cox.zph(fit)$table, error = function(e) NULL)
    pairs$model <- model_name; pairs$endpoint <- endpoint; pairs$n <- nrow(z); pairs$events <- sum(z$event)
    pairs$global_LRT_p <- lrt$`Pr(>|Chi|)`[2]
    pairs$PH_global_p <- if (is.null(ph)) NA_real_ else ph["GLOBAL", "p"]
    pairs$comparison_group <- pairs$group1; pairs$reference_group <- pairs$group2
    pairs$ratio_definition <- "Hazard in comparison_group / hazard in reference_group"
    pairs$direction <- ratio_direction(pairs$HR, pairs$group1, pairs$group2,
      ifelse(endpoint == "OS", "all-cause death", endpoint), pairs$p_holm < .05)
    pairs$formula <- paste(deparse(form), collapse = " ")
    pairs
  }

  # Figure 3: platinum exposure. Unknown is excluded, and the within-neoadjuvant
  # model is shown because platinum exposure is structurally nested in treatment.
  platinum <- d[d$Preop_platinum_exposure %in% c("No", "Yes"), , drop = FALSE]
  platinum$Preop_platinum <- factor(platinum$Preop_platinum_exposure, levels = c("No", "Yes"))
  platinum_vars <- c("Age", "Sex", "T_stage", "N_category", "M_stage", "Differentiation",
    "LVI", "PNI", "R_status", "KRAS_group", "TP53", "SMAD4", "CDKN2A", "TMB", "Variant_count")
  pc <- compare_groups(platinum, "Preop_platinum", platinum_vars, "Preop_platinum")
  tab(pc$global, "Input_Preop_platinum_global_comparisons.tsv")
  tab(pc$pairwise, "Input_Preop_platinum_pairwise_comparisons.tsv")
  tab(describe_by_group(platinum, "Preop_platinum", platinum_vars),
    "Input_Preop_platinum_descriptive.tsv")
  platinum_summary <- platinum |> dplyr::count(Preop_platinum, Neoadjuvant, name = "n")
  tab(platinum_summary, "Input_Preop_platinum_population.tsv")
  pp <- association_plot(pc$global, "Clinical and genomic context of preoperative platinum",
    sprintf("Yes n=%d; No n=%d; Unknown excluded n=%d", sum(platinum$Preop_platinum == "Yes"),
      sum(platinum$Preop_platinum == "No"), sum(d$Preop_platinum_exposure == "Unknown")))
  pkm <- km_analysis(platinum, "Preop_platinum", "OS_m", "survive",
    "Overall survival by preoperative platinum exposure", "Overall survival", "Preop_platinum_OS")
  neo <- platinum[platinum$Neoadjuvant == "y", , drop = FALSE]
  covs <- c("splines::ns(Age, df = 3)", "Sex", "T_stage", "N_stage", "M_stage",
    "strata(Differentiation)", "LVI", "PNI", "R_status")
  pcx <- dplyr::bind_rows(
    cox_row(platinum, "Preop_platinum", "OS_m", "survive", character(), "All evaluable; unadjusted", "OS"),
    cox_row(neo, "Preop_platinum", "OS_m", "survive", covs, "Neoadjuvant-treated; adjusted", "OS"),
    cox_row(platinum, "Preop_platinum", "RFS_m", "Recur", character(), "All evaluable; unadjusted", "Recorded recurrence"),
    cox_row(neo, "Preop_platinum", "RFS_m", "Recur", covs, "Neoadjuvant-treated; adjusted", "Recorded recurrence"))
  pcx$p_holm_all_models <- p.adjust(pcx$p, "holm")
  tab(pcx, "Input_Preop_platinum_Cox_models.tsv")
  pforest <- forest_result_plot(pcx, "Preoperative platinum and outcomes",
    paste0(pcx$endpoint, ": Yes\nReference: No; ", pcx$model), p_column = "p_holm_all_models", p_heading = "Holm p") +
    ggplot2::labs(caption = "HR = Yes / No. HR > 1 indicates a higher event hazard in the platinum-exposed group.\nThe adjusted models are restricted to neoadjuvant-treated patients; observational associations are not treatment effects.")
  savefig(pp, "Figure3_Platinum_context", 10, 8)
  savefig(pkm$plot, "Figure3_Platinum_OS", 10, 7.5)
  savefig(pforest, "Figure3_Platinum_forest", 12, 5.5)
  savefig(cowplot::plot_grid(pp, pkm$plot, pforest, ncol = 1, rel_heights = c(.37, .38, .25)),
    "Figure3_Preoperative_platinum", 13, 18)

  # Figure 4: recurrence pattern is defined only after recurrence. Blank values
  # among non-recurrent patients are excluded, not treated as a fourth category.
  recur <- d[d$Recur == 1 & !is.na(d$Recurrence_pattern), , drop = FALSE]
  recur$Recurrence_group <- factor(recur$Recurrence_pattern,
    levels = c("Local only", "Distant only", "Local + distant"))
  recur_vars <- c("Age", "Sex", "T_stage", "N_category", "Differentiation", "LVI", "PNI",
    "R_status", "KRAS_group", "TP53", "SMAD4", "CDKN2A", "TMB", "Variant_count")
  rc <- compare_groups(recur, "Recurrence_group", recur_vars, "Recurrence_pattern")
  tab(rc$global, "Input_Recurrence_pattern_global_comparisons.tsv")
  tab(rc$pairwise, "Input_Recurrence_pattern_pairwise_comparisons.tsv")
  tab(describe_by_group(recur, "Recurrence_group", recur_vars),
    "Input_Recurrence_pattern_descriptive.tsv")
  rdist <- recur |> dplyr::count(Recurrence_group, name = "n") |>
    dplyr::mutate(percent = 100 * n / sum(n))
  tab(rdist, "Input_Recurrence_pattern_distribution.tsv")
  rbar <- ggplot2::ggplot(rdist, ggplot2::aes(Recurrence_group, percent, fill = Recurrence_group)) +
    ggplot2::geom_col(width = .66) + ggplot2::geom_text(ggplot2::aes(label = sprintf("%d (%.1f%%)", n, percent)), vjust = -.35) +
    ggplot2::scale_fill_manual(values = setNames(pal[2:4], levels(recur$Recurrence_group))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .14))) +
    ggplot2::labs(x = NULL, y = "Patients among recorded recurrences (%)", title = "Recorded recurrence patterns",
      subtitle = sprintf("Recurrent patients with a recorded pattern: N=%d", nrow(recur)),
      caption = "Non-recurrent blanks and patients without an evaluable recurrence endpoint are not a fourth group.") +
    ggplot2::theme(legend.position = "none")
  rp <- association_plot(rc$global, "Clinicopathologic and genomic differences across recurrence patterns",
    "Global tests across Local only, Distant only, and Local + distant; all pairwise tests use Holm correction")
  recur$post_recurrence_months <- recur$OS_m - recur$RFS_m
  nonpositive <- sum(!is.na(recur$post_recurrence_months) & recur$post_recurrence_months <= 0)
  post <- recur[is.na(recur$post_recurrence_months) | recur$post_recurrence_months > 0, , drop = FALSE]
  rkm <- km_analysis(post, "Recurrence_group", "post_recurrence_months", "survive",
    "Survival after recorded recurrence", "Post-recurrence overall survival", "Recurrence_pattern_post_OS")
  rcx <- cox_row(post, "Recurrence_group", "post_recurrence_months", "survive",
    c("splines::ns(Age, df = 3)", "Sex", "T_stage", "N_stage", "M_stage",
      "strata(Differentiation)", "LVI", "PNI", "R_status"),
    "Adjusted post-recurrence survival", "death after recorded recurrence")
  tab(rcx, "Input_Recurrence_pattern_Cox_models.tsv")
  rforest <- forest_result_plot(rcx, "Adjusted post-recurrence comparisons",
    paste0(rcx$group1, "\nReference: ", rcx$group2), p_column = "p_holm", p_heading = "Holm p") +
    ggplot2::labs(caption = paste0("HR = first pattern / reference pattern; HR > 1 indicates higher death hazard after recorded recurrence.\n",
      "Exploratory post-baseline analysis; nonpositive post-recurrence intervals excluded n=", nonpositive, "."))
  savefig(rbar, "Figure4_Recurrence_distribution", 8, 5)
  savefig(rp, "Figure4_Recurrence_context", 10, 8)
  savefig(rkm$plot, "Figure4_Post_recurrence_OS", 10, 7.5)
  savefig(rforest, "Figure4_Post_recurrence_forest", 12, 5)
  fig4 <- cowplot::plot_grid(cowplot::plot_grid(rbar, rp, ncol = 2, rel_widths = c(.38, .62)),
    rkm$plot, rforest, ncol = 1, rel_heights = c(.34, .39, .27))
  savefig(fig4, "Figure4_Recurrence_phenotype", 15, 19)

  # Figure 5 consolidates the prognostic evidence instead of presenting every
  # available plot as a main figure.
  fig5 <- cowplot::plot_grid(fg, kk, kf, ncol = 1, rel_heights = c(.24, .42, .34))
  savefig(fig5, "Figure5_Genomic_prognosis", 14, 18)

  # Main publication tables are also rendered as PNG. Detailed calculation
  # files remain Input_*.tsv and are not labelled supplementary tables.
  t1 <- readr::read_tsv(file.path(out, "tables", "Table1_Cohort_characteristics.tsv"), show_col_types = FALSE)
  t1shown <- t1 |> dplyr::transmute(Characteristic = clean_label(variable),
    Category = dplyr::case_when(category == "y" ~ "Yes", category == "n" ~ "No", TRUE ~ category),
    `Value, n (%) or median [IQR]` = value, `Observed n` = observed, `Missing n` = missing)
  save_table_png(t1shown, "Table 1. Cohort characteristics", "Table1_Cohort_characteristics",
    widths = c(2.1, 2.1, 2.5, 1, 1), font_size = 7.8, width = 15)

  gene_rows <- core |> dplyr::filter(exposure %in% c("KRAS", "TP53", "SMAD4", "CDKN2A")) |>
    dplyr::transmute(endpoint = "OS", variable = exposure, comparison_group = "Variant detected",
      reference_group = "Not detected", n, events, HR, lower95, upper95, raw_p = p,
      adjusted_p = q_BH, adjustment = "BH across exposure coefficients", method = "Adjusted Cox",
      direction = interpretation)
  kras_context <- core |> dplyr::filter(exposure == "KRAS_group") |> dplyr::slice(1)
  kras_rows <- cpair |> dplyr::filter(exposure == "KRAS_group", family == "Extended genomic", endpoint == "OS") |>
    dplyr::transmute(endpoint = "OS", variable = "KRAS subtype", comparison_group = group1,
      reference_group = group2, n = kras_context$n[[1]], events = kras_context$events[[1]], HR, lower95, upper95,
      raw_p = p, adjusted_p = p_holm, adjustment = "Holm across 10 subtype pairs",
      method = "Adjusted Cox contrast", direction = interpretation)
  platinum_rows <- pcx |> dplyr::transmute(endpoint, variable = "Preoperative platinum", comparison_group,
    reference_group, n, events, HR, lower95, upper95, raw_p = p, adjusted_p = p_holm_all_models,
    adjustment = "Holm across 4 prespecified models", method = model, direction)
  recurrence_rows <- rcx |> dplyr::transmute(endpoint = "Post-recurrence OS", variable = "Recurrence pattern",
    comparison_group, reference_group, n, events, HR, lower95, upper95, raw_p = p,
    adjusted_p = p_holm, adjustment = "Holm across 3 pattern pairs", method = model, direction)
  main_survival <- dplyr::bind_rows(gene_rows, kras_rows, platinum_rows, recurrence_rows)
  tab(main_survival, "Table2_Main_survival_models.tsv")
  display <- main_survival |> dplyr::transmute(Endpoint = endpoint, Variable = variable,
    Comparison = paste0(comparison_group, " / ", reference_group), `n / events` = ifelse(is.na(n), "See model input", paste0(n, " / ", events)),
    `HR [95% CI]` = sprintf("%.2f [%.2f, %.2f]", HR, lower95, upper95),
    `Raw p` = fmt_p(raw_p), `Adjusted p` = fmt_p(adjusted_p), Adjustment = adjustment)
  save_table_png(display, "Table 2. Main survival models and pairwise contrasts", "Table2_Main_survival_models",
    widths = c(1.2, 1.5, 2.2, 1.2, 1.7, .9, 1.1, 2.2), font_size = 7.3, width = 18)

  invisible(list(platinum = pc, platinum_cox = pcx, recurrence = rc,
    recurrence_cox = rcx, table2 = main_survival))
}
