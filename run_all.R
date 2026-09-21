#!/usr/bin/env Rscript

# Run every PAAD analysis in a clear, documented order.
#
# Each file under scripts/ is a complete standalone analysis. This runner starts
# a fresh R session for each one so the combined run uses exactly the same code
# as an individual run and never relies on objects left by another analysis.
#
# Run from the repository root:
#   Rscript --vanilla run_all.R --config=config/local.R
#   Rscript --vanilla run_all.R --config=config/config.synthetic.R

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

analysis_scripts <- c(
  "scripts/oncoplots.R",
  "scripts/clinical_summary.R",
  "scripts/clinical_group_comparisons.R",
  "scripts/mutation_heatmaps.R",
  "scripts/somatic_interactions.R",
  "scripts/stage_distribution.R",
  "scripts/driver_kras_survival.R"
)

missing_scripts <- analysis_scripts[!file.exists(file.path(project_root, analysis_scripts))]
if (length(missing_scripts)) {
  stop(
    "Run this file from the PAAD repository root. Missing: ",
    paste(missing_scripts, collapse = ", "),
    call. = FALSE
  )
}

arguments <- commandArgs(trailingOnly = TRUE)
config_arguments <- grep("^--config=", arguments, value = TRUE)
if (length(config_arguments) > 1L) {
  stop("Specify --config only once.", call. = FALSE)
}

config_file <- if (length(config_arguments)) {
  sub("^--config=", "", config_arguments[[1]])
} else {
  "config/local.R"
}
if (!grepl("^(/|[A-Za-z]:[/\\])", config_file)) {
  config_file <- file.path(project_root, config_file)
}
config_file <- normalizePath(config_file, winslash = "/", mustWork = TRUE)

rscript <- file.path(R.home("bin"), "Rscript")

for (index in seq_along(analysis_scripts)) {
  script <- analysis_scripts[[index]]
  message(
    "\n[", index, "/", length(analysis_scripts), "] Running ", script
  )

  status <- system2(
    command = rscript,
    args = c(
      "--vanilla",
      shQuote(file.path(project_root, script)),
      shQuote(paste0("--config=", config_file))
    )
  )

  if (!identical(as.integer(status), 0L)) {
    stop("Analysis failed: ", script, " (exit status ", status, ")", call. = FALSE)
  }
}

message("\nAll PAAD analyses completed successfully.")
