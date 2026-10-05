source("R/config.R")
source("R/pairwise_tests.R")
config <- load_config_from_command_line()
read <- function(f) readr::read_tsv(output_path(config, f), show_col_types = FALSE)
families <- list(
  "26_KRAS_pairwise_logrank.tsv" = c("variable", "endpoint"),
  "21_OS_pairwise_logrank.tsv" = "variable",
  "20_Cox_pairwise_contrasts.tsv" = "model",
  "29_Source_Cox_pairwise_contrasts.tsv" = "model",
  "19_Pathology_pairwise_CMH.tsv" = c("gene", "characteristic"),
  "13_Driver_pathology_pairwise_Fisher.tsv" = c("gene", "characteristic"),
  "27_Burden_pairwise_tests.tsv" = c("metric", "variable"),
  "28_Clinical_pairwise_tests.tsv" = c("variable", "grouping"),
  "21_RMST_pairwise_contrasts.tsv" = "variable",
  "32_Repair_group_clinical_pairwise_tests.tsv" = c("grouping", "variable"))
for (name in names(families)) {
  d <- read(name)
  id <- do.call(paste, c(d[families[[name]]], sep = "|"))
  for (z in split(d, id)) {
    stopifnot(isTRUE(all.equal(z$p_holm, p.adjust(z$p, "holm", n = nrow(z)))))
    k <- length(unique(c(z$group1, z$group2)))
    stopifnot(nrow(z) == choose(k, 2), !anyDuplicated(paste(z$group1, z$group2)))
  }
}
d <- read("29_Integrated_analysis_dataset_private.tsv")
p <- read("26_KRAS_pairwise_logrank.tsv")
for (v in unique(p$variable)) for (e in unique(p$endpoint)) {
  time <- if (e == "OS") d$OS_months else d$RFS_m
  event <- if (e == "OS") d$survive else d$Recur
  got <- p[p$variable == v & p$endpoint == e, ]
  expected <- pairwise_logrank(time, event, factor(d[[v]], levels = unique(c(got$group1, got$group2))))
  key <- paste(expected$group1, expected$group2)
  idx <- match(paste(got$group1, got$group2), key)
  stopifnot(!anyNA(idx), all(abs(got$p - expected$p[idx]) < 1e-12),
    all(got$n1 == expected$n1[idx]), all(got$events2 == expected$events2[idx]))
}
for (f in list.files(config$output_dir, pattern = "^4_.*pairwise_.*[.]tsv$")) {
  d <- read(f)
  for (z in if ("gene" %in% names(d)) split(d, d$gene) else list(d)) {
    stopifnot(isTRUE(all.equal(z$p_holm, p.adjust(z$p, "holm", n = nrow(z)))))
  }
}
cat("All pairwise output families validated: complete pairs, Holm, actual-source log-rank n/events/p.\n")
