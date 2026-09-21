#!/usr/bin/env Rscript

# Dependency-light integration checks for the public PAAD repository. This file
# uses base R assertions so testthat is not required.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
setwd(project_root)

failures <- character()
run_test <- function(name, code) {
  tryCatch(
    {
      force(code)
      message("PASS: ", name)
    },
    error = function(error) {
      failures <<- c(failures, paste0(name, ": ", conditionMessage(error)))
      message("FAIL: ", name, " — ", conditionMessage(error))
    }
  )
}

expect_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
  invisible(TRUE)
}

expect_error <- function(code, pattern = NULL) {
  error <- tryCatch({
    force(code)
    NULL
  }, error = identity)
  if (is.null(error)) stop("Expected an error, but the expression succeeded.", call. = FALSE)
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(error), fixed = TRUE)) {
    stop("Unexpected error: ", conditionMessage(error), call. = FALSE)
  }
  invisible(error)
}

command_status <- function(command, args) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = output)
}

required_files <- c(
  "R/bootstrap.R", "R/config.R", "R/data.R", "R/plot_helpers.R",
  "R/plot_genomics.R", "R/driver_survival.R", "R/workflows.R",
  "scripts/01_oncoplots.R", "scripts/02_clinical_summary.R",
  "scripts/03_comparison_plots.R", "scripts/04_heatmaps.R",
  "scripts/05_somatic_interactions.R", "scripts/06_stage_counts.R",
  "scripts/07_driver_kras_survival.R", "run_all.R",
  "config/config.example.R", "config/config.synthetic.R",
  "data/example/PDAC_clinical_info.example.xlsx",
  "data/example/result_527명.example.txt",
  "data/example/PDAC_oncopanel.example.maf",
  "data/example/generate_examples.py"
)

run_test("required public files exist", {
  missing <- required_files[!file.exists(required_files)]
  expect_true(!length(missing), paste("Missing:", paste(missing, collapse = ", ")))
})

run_test("all R source and entrypoint files parse", {
  files <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
             list.files("scripts", pattern = "[.]R$", full.names = TRUE),
             "run_all.R")
  invisible(lapply(files, parse))
})

run_test("synthetic files match their deterministic generator", {
  result <- command_status("python3", c("data/example/generate_examples.py", "--check"))
  expect_true(result$status == 0L, paste(result$output, collapse = "\n"))
})

run_test("local/private paths are ignored but synthetic fixtures are not", {
  ignored <- c("data/raw/ignore-probe.txt", "config/local.R", "outputs/ignore-probe.txt", "local_notes/ignore-probe.txt")
  for (path in ignored) {
    result <- command_status("git", c("check-ignore", "-q", path))
    expect_true(result$status == 0L, paste(path, "must be ignored"))
  }
  public_fixture <- command_status("git", c("check-ignore", "-q", "data/example/PDAC_oncopanel.example.maf"))
  expect_true(public_fixture$status != 0L, "Synthetic fixtures must remain eligible for Git tracking.")
})

run_test("public text files contain no personal absolute home path", {
  roots <- c("R", "scripts", "config", "docs", "tests", "data/example")
  files <- unlist(lapply(roots, function(path) {
    list.files(path, pattern = "([.]R|[.]md|[.]py)$", recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  files <- c(files, "README.md", "data/README.md", ".gitignore", ".gitattributes", "PAAD.Rproj", "run_all.R")
  files <- unique(files[file.exists(files)])
  offenders <- vapply(files, function(path) {
    text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    grepl("/(Users|home)/[A-Za-z0-9._-]+/", text, perl = TRUE)
  }, logical(1))
  expect_true(!any(offenders), paste("Absolute home path found in:", paste(files[offenders], collapse = ", ")))
})

source("R/bootstrap.R")
source_project_modules(project_root)

run_test("configuration rejects output outside repository outputs", {
  safe <- load_config(file.path(project_root, "config", "config.synthetic.R"), project_root)
  unsafe <- unclass(safe)
  unsafe$output_dir <- tempdir()
  expect_error(validate_config(unsafe, project_root), "output_dir must be outputs/")
  expect_error(output_path(safe, "../escape.txt"), "directory separator")
})

synthetic_config <- NULL
synthetic_data <- NULL
run_test("synthetic XLSX, CP949 cohort list, and MAF load together", {
  synthetic_config <<- load_config(file.path(project_root, "config", "config.synthetic.R"), project_root)
  synthetic_data <<- load_analysis_data(synthetic_config)
  expect_true(nrow(synthetic_data$clinical) == 20L, "Expected 20 synthetic clinical rows after PDAC filtering.")
  expect_true(nrow(synthetic_data$maf@data) == 300L, "Expected 300 synthetic retained MAF rows.")
  expect_true(all(grepl("^SYN-TUMOR-", synthetic_data$clinical$Tumor_Sample_Barcode)), "Clinical IDs are not wholly synthetic.")
  expect_true(all(grepl("^SYN-TUMOR-", synthetic_data$maf@data$Tumor_Sample_Barcode)), "MAF IDs are not wholly synthetic.")
  required <- c("Age", "Age_group", "Stage_Group", "OS_m", "RFS_m", "survive", "Recur", "LVI", "PNI", "RM")
  expect_true(all(required %in% names(synthetic_data$clinical)), "Expected derived clinical columns are missing.")
  expect_true(all(stats::na.omit(synthetic_data$clinical$survive) %in% c(0, 1)), "OS event coding must be binary in the fixture.")
  expect_true(all(stats::na.omit(synthetic_data$clinical$Recur) %in% c(0, 1)), "RFS event coding must be binary in the fixture.")
})

run_test("core grouping and driver summaries work on synthetic data", {
  expect_true(identical(get_expected_levels("N", c("N1", "N0")), c("N0", "N1", "N2")), "N level order changed.")
  driver <- make_driver_summary(synthetic_data$maf@data, synthetic_data$clinical)
  expect_true(nrow(driver) == nrow(synthetic_data$clinical), "Driver summary changed the cohort row count.")
  expect_true(all(c("KRAS", "TP53", "SMAD4", "CDKN2A", "KRAS_subtype", "driver_mutation_count") %in% names(driver)), "Driver summary columns are incomplete.")
  expect_true(all(driver$driver_mutation_count >= 0L & driver$driver_mutation_count <= 4L), "Driver count is outside 0–4.")
})

run_test("lightweight workflows write only under outputs", {
  smoke_dir <- file.path(project_root, "outputs", "test_smoke")
  dir.create(smoke_dir, recursive = TRUE, showWarnings = FALSE)
  smoke_config <- synthetic_config
  smoke_config$output_dir <- normalizePath(smoke_dir, winslash = "/", mustWork = TRUE)
  workflow_clinical_summary(smoke_config, synthetic_data)
  workflow_stage_counts(smoke_config, synthetic_data)
  expected <- file.path(smoke_dir, c("2_Clinical_Table_PDAC.xlsx", "10_Stage_Counts_PDAC.png"))
  expect_true(all(file.exists(expected)), "Lightweight workflow outputs were not created.")
  expect_true(all(file.info(expected)$size > 0), "A workflow output is empty.")
})

run_test("full CLI workflow completes with the synthetic config", {
  rscript <- file.path(R.home("bin"), "Rscript")
  result <- command_status(rscript, c("--vanilla", "run_all.R", "--config=config/config.synthetic.R"))
  expect_true(result$status == 0L, paste(result$output, collapse = "\n"))
  output_files <- list.files("outputs/synthetic_run", full.names = TRUE)
  representatives <- file.path("outputs/synthetic_run", c(
    "1_Oncoplot_Main_PDAC.png", "2_Clinical_Table_PDAC.xlsx",
    "4_FreqPlot_PDAC_Stage.png", "5_Heatmap_Top20_PDAC.png",
    "9_Somatic_Interactions_Results_PDAC.xlsx", "10_Stage_Counts_PDAC.png",
    "11_Driver_KRAS_patient_level_summary.xlsx", "17_RFS_by_driver_mutation_count_exact.png"
  ))
  expect_true(length(output_files) == 42L, "Full synthetic run did not produce the expected 42 files.")
  expect_true(all(file.exists(representatives)), "One or more analysis workflow outputs are missing.")
})

if (length(failures)) {
  cat("\nTest failures:\n", paste0("- ", failures, collapse = "\n"), "\n", sep = "")
  quit(status = 1L)
}

cat("\nAll PAAD repository tests passed.\n")
