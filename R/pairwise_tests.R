# Explicit all-pairs tests. Holm families are one plot/variable/model, not the
# whole study. All planned pairs (including unestimable pairs) remain in tables.
# No filtering on an omnibus p-value. No patient identifiers are exported.
pair_frame <- function(group) {
  g <- droplevels(factor(group))
  if (nlevels(g) < 2) return(data.frame())
  pairs <- t(combn(levels(g), 2))
  data.frame(group1 = pairs[, 1], group2 = pairs[, 2],
    n1 = as.integer(table(g)[pairs[, 1]]), n2 = as.integer(table(g)[pairs[, 2]]),
    statistic = NA_real_, df = NA_real_, p = NA_real_, status = "Not estimable")
}
finish_pairs <- function(x, method) {
  if (!nrow(x)) return(x)
  x$p_holm <- p.adjust(x$p, "holm", n = nrow(x))
  x$method <- method
  x$family_size <- nrow(x)
  x
}
pairwise_logrank <- function(time, event, group) {
  z <- data.frame(time, event, group)
  z <- z[complete.cases(z) & is.finite(time) & time > 0 & event %in% 0:1, ]
  out <- pair_frame(z$group)
  if (!nrow(out)) return(out)
  out$events1 <- out$events2 <- NA_integer_
  for (i in seq_len(nrow(out))) {
    a <- out$group1[i]; b <- out$group2[i]
    zz <- z[z$group %in% c(a, b), ]; zz$group <- droplevels(factor(zz$group))
    out$events1[i] <- sum(zz$event[zz$group == a])
    out$events2[i] <- sum(zz$event[zz$group == b])
    if (!sum(zz$event)) { out$status[i] <- "No events"; next }
    f <- tryCatch(survival::survdiff(survival::Surv(time, event) ~ group, data = zz), error = identity)
    if (inherits(f, "error")) { out$status[i] <- conditionMessage(f); next }
    out$statistic[i] <- f$chisq; out$df[i] <- 1
    out$p[i] <- pchisq(f$chisq, 1, lower.tail = FALSE)
    out$status[i] <- if (is.finite(out$p[i])) "Estimated" else "Degenerate variance"
  }
  finish_pairs(out, "Pairwise log-rank; chi-square, df=1")
}
pairwise_numeric <- function(value, group) {
  ok <- complete.cases(value, group) & is.finite(value)
  x <- value[ok]; g <- droplevels(factor(group[ok])); out <- pair_frame(g)
  if (!nrow(out)) return(out)
  if (nlevels(g) == 2) {
    f <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))
    out$statistic <- unname(f$statistic); out$p <- f$p.value
    method <- "Wilcoxon rank-sum W; two-sided, continuity correction"
  } else {
    # Dunn Z based on ranks pooled across ALL groups, with tie correction.
    n <- length(x); ties <- as.numeric(table(x))
    variance <- n * (n + 1) / 12 - sum(ties^3 - ties) / (12 * (n - 1))
    ranks <- tapply(rank(x), g, mean)
    out$statistic <- as.numeric((ranks[out$group1] - ranks[out$group2]) /
      sqrt(variance * (1 / out$n1 + 1 / out$n2)))
    out$p <- 2 * pnorm(-abs(out$statistic))
    method <- "Dunn Z; two-sided, pooled ranks, ties corrected"
  }
  out$status <- ifelse(is.finite(out$p), "Estimated", "Degenerate variance")
  finish_pairs(out, method)
}
pairwise_categorical <- function(value, group, strata = NULL) {
  ok <- complete.cases(value, group)
  if (!is.null(strata)) ok <- ok & !is.na(strata)
  x <- as.character(value[ok]); g <- droplevels(factor(group[ok]))
  ss <- if (!is.null(strata)) droplevels(factor(strata[ok])) else NULL
  out <- pair_frame(g)
  if (!nrow(out)) return(out)
  out$odds_ratio <- out$lower95 <- out$upper95 <- NA_real_
  out$outcome_first_level <- out$outcome_second_level <- NA_character_
  for (i in seq_len(nrow(out))) {
    keep <- g %in% c(out$group1[i], out$group2[i])
    xx <- droplevels(factor(x[keep])); gg <- factor(g[keep], levels = c(out$group1[i], out$group2[i]))
    if (nlevels(xx) < 2) { out$status[i] <- "Only one outcome category"; next }
    arr <- if (is.null(ss)) table(xx, gg) else table(xx, gg, droplevels(ss[keep]))
    f <- tryCatch(if (is.null(ss)) fisher.test(arr, simulate.p.value = nrow(arr) > 2, B = 10000)
      else mantelhaen.test(arr, correct = FALSE), error = identity)
    if (inherits(f, "error")) { out$status[i] <- conditionMessage(f); next }
    out$p[i] <- f$p.value; out$status[i] <- "Estimated"
    if (!is.null(f$statistic)) out$statistic[i] <- unname(f$statistic)
    if (!is.null(f$parameter)) out$df[i] <- unname(f$parameter)
    if (!is.null(f$estimate)) {
      out$odds_ratio[i] <- unname(f$estimate)
      out$lower95[i] <- f$conf.int[1]; out$upper95[i] <- f$conf.int[2]
      out$outcome_first_level[i] <- levels(xx)[1]; out$outcome_second_level[i] <- levels(xx)[2]
    }
  }
  finish_pairs(out, if (is.null(ss)) "Fisher exact (2x2); Monte Carlo 10000 (>2 outcome levels)"
    else "Pairwise CMH, stratified; no continuity correction")
}
pairwise_cox <- function(fit, data, variable) {
  if (is.numeric(data[[variable]]) || length(unique(data[[variable]])) < 3) return(data.frame())
  out <- pair_frame(data[[variable]])
  out$HR <- out$lower95 <- out$upper95 <- NA_real_
  nd <- data[rep(1L, nlevels(factor(data[[variable]]))), , drop = FALSE]
  lev <- levels(droplevels(factor(data[[variable]])))
  nd[[variable]] <- factor(lev, levels = lev)
  # Cox method preserves training spline knots, factor levels and strata.
  mm <- model.matrix(fit, data = nd)
  mm <- mm[, names(coef(fit)), drop = FALSE]
  rownames(mm) <- lev
  for (i in seq_len(nrow(out))) {
    delta <- mm[out$group1[i], ] - mm[out$group2[i], ]
    take <- which(delta != 0)
    if (!length(take) || any(!is.finite(coef(fit)[take]))) next
    estimate <- sum(delta[take] * coef(fit)[take])
    se <- sqrt(as.numeric(t(delta[take]) %*% vcov(fit)[take, take, drop = FALSE] %*% delta[take]))
    if (!is.finite(se) || se <= 0) next
    out$statistic[i] <- estimate / se
    out$p[i] <- 2 * pnorm(-abs(out$statistic[i]))
    out$HR[i] <- exp(estimate)
    out$lower95[i] <- exp(estimate - qnorm(.975) * se)
    out$upper95[i] <- exp(estimate + qnorm(.975) * se)
    out$status[i] <- "Estimated; HR is group1 relative to group2"
  }
  finish_pairs(out, "Cox Wald Z contrast; pointwise 95% CI, Holm-adjusted p")
}
pairwise_caption <- function(x, width = 110) {
  if (!nrow(x)) return("Pairwise comparisons not estimable.")
  sig <- x[is.finite(x$p_holm) & x$p_holm < .05, ]
  if (!nrow(sig)) return("Pairwise comparisons: none significant after Holm correction.")
  s <- paste0(sig$group1, " vs ", sig$group2, ": p=", formatC(sig$p_holm, digits = 3, format = "g"))
  paste(strwrap(paste("Significant pairs (Holm-adjusted p):", paste(s, collapse = "; ")), width = width), collapse = "\n")
}
