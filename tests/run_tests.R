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

standalone_outputs <- list(
  "scripts/oncoplots.R" = c(
    "1_Oncoplot_Main_PDAC.png",
    "1_2_Oncoplot_TN_Filtered_PDAC.png"
  ),
  "scripts/clinical_summary.R" = "2_Clinical_Table_PDAC.xlsx",
  "scripts/clinical_group_comparisons.R" = c(
    paste0("4_FreqPlot_PDAC_", c("Differentiation", "NAC", "T", "N", "N_status", "Stage", "Stage_Group"), ".png"),
    paste0("4_TMB_BoxPlot_PDAC_", c("Differentiation", "NAC", "T", "N", "N_status", "Stage", "Stage_Group"), ".png")
  ),
  "scripts/mutation_heatmaps.R" = c(
    "5_Heatmap_Top20_PDAC.png",
    "6_Heatmap_Top5_PDAC.png",
    "7_TMB_Heatmap_PDAC.png"
  ),
  "scripts/somatic_interactions.R" = c(
    "9_Somatic_Interactions_PDAC.png",
    "9_Somatic_Interactions_Results_PDAC.xlsx"
  ),
  "scripts/stage_distribution.R" = "10_Stage_Counts_PDAC.png",
  "scripts/driver_kras_survival.R" = c(
    "11_Driver_KRAS_patient_level_summary.xlsx",
    "11_KRAS_subtype_pie.png",
    "12_Driver_mutation_status_pie.png",
    "13_Driver_gene_mutation_clinicopath_table.tsv",
    "13_Driver_gene_mutation_clinicopath_table.xlsx",
    "13_Driver_gene_mutation_clinicopath_table_display.xlsx",
    "14_17_PDAC_survival_input_patient_level.xlsx",
    "14_OS_by_KRAS_subtype.png",
    "14_OS_by_KRAS_subtype_number_at_risk.tsv",
    "15_OS_by_driver_mutation_count_exact.png",
    "15_OS_by_driver_mutation_count_exact_number_at_risk.tsv",
    "15_2_OS_by_driver_mutation_count_0_12_34.png",
    "15_2_OS_by_driver_mutation_count_0_12_34_number_at_risk.tsv",
    "16_RFS_by_KRAS_subtype.png",
    "16_RFS_by_KRAS_subtype_number_at_risk.tsv",
    "17_RFS_by_driver_mutation_count_exact.png",
    "17_RFS_by_driver_mutation_count_exact_number_at_risk.tsv",
    "17_2_RFS_by_driver_mutation_count_0_12_34.png",
    "17_2_RFS_by_driver_mutation_count_0_12_34_number_at_risk.tsv"
  )
)
expected_outputs <- unname(unlist(standalone_outputs, use.names = FALSE))

required_files <- c(
  "R/config.R", "R/data.R", "R/plot_helpers.R",
  names(standalone_outputs), "scripts/README.md", "run_all.R",
  "config/config.example.R", "config/config.synthetic.R",
  "data/example/PDAC_clinical_info.example.xlsx",
  "data/example/result_527명.example.txt",
  "data/example/PDAC_oncopanel.example.maf",
  "data/example/generate_examples.py",
  "docs/analysis_map.md", "docs/refactoring_notes.md"
)

retired_files <- c(
  "R/bootstrap.R", "R/workflows.R", "R/plot_genomics.R", "R/driver_survival.R",
  "scripts/01_oncoplots.R", "scripts/02_clinical_summary.R",
  "scripts/03_comparison_plots.R", "scripts/04_heatmaps.R",
  "scripts/05_somatic_interactions.R", "scripts/06_stage_counts.R",
  "scripts/07_driver_kras_survival.R"
)

run_test("required public files exist", {
  missing <- required_files[!file.exists(required_files)]
  expect_true(!length(missing), paste("Missing:", paste(missing, collapse = ", ")))
  expect_true(length(expected_outputs) == 42L, "The documented output contract must contain 42 files.")
  expect_true(!anyDuplicated(expected_outputs), "The output contract contains duplicate filenames.")
})

run_test("retired dispatch files and numbered entrypoints are absent", {
  remaining <- retired_files[file.exists(retired_files)]
  numbered <- list.files("scripts", pattern = "^[0-9]{2}_.*[.]R$", full.names = TRUE)
  expect_true(!length(c(remaining, numbered)), paste("Remove retired files:", paste(c(remaining, numbered), collapse = ", ")))
})

run_test("standalone files contain their analysis code directly", {
  for (path in names(standalone_outputs)) {
    text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    expect_true(grepl("load_config", text, fixed = TRUE), paste(path, "does not load its configuration directly."))
    expect_true(grepl("load_analysis_data", text, fixed = TRUE), paste(path, "does not load its analysis data directly."))
    expect_true(grepl("function[[:space:]]*[(]", text), paste(path, "does not expose an analysis function."))
    expect_true(!grepl("run_entrypoint", text, fixed = TRUE), paste(path, "still forwards through the retired dispatcher."))
  }
})

run_test("all public R files parse", {
  files <- c(
    list.files("R", pattern = "[.]R$", full.names = TRUE),
    list.files("scripts", pattern = "[.]R$", full.names = TRUE),
    "run_all.R", "tests/run_tests.R"
  )
  invisible(lapply(files, parse))
})

run_test("public documentation has no retired file references", {
  files <- c("README.md", "scripts/README.md", list.files("docs", pattern = "[.]md$", full.names = TRUE))
  retired_references <- c(
    "R/bootstrap.R", "R/workflows.R",
    "scripts/01_", "scripts/02_", "scripts/03_", "scripts/04_",
    "scripts/05_", "scripts/06_", "scripts/07_"
  )
  offenders <- character()
  for (path in files) {
    text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (any(vapply(retired_references, grepl, logical(1), x = text, fixed = TRUE))) offenders <- c(offenders, path)
  }
  expect_true(!length(offenders), paste("Retired paths remain in:", paste(unique(offenders), collapse = ", ")))
})

run_test("synthetic files match their deterministic generator", {
  result <- command_status("python3", c("data/example/generate_examples.py", "--check"))
  expect_true(result$status == 0L, paste(result$output, collapse = "\n"))
})

run_test("private inputs and outputs are ignored but synthetic fixtures are public", {
  ignored <- c(
    "data/raw/ignore-probe.txt", "data/processed/ignore-probe.txt", "data/private/ignore-probe.txt",
    "config/local.R", "outputs/ignore-probe.txt", "local_notes/ignore-probe.txt"
  )
  for (path in ignored) {
    result <- command_status("git", c("check-ignore", "-q", path))
    expect_true(result$status == 0L, paste(path, "must be ignored"))
  }

  public_fixture <- command_status("git", c("check-ignore", "-q", "data/example/PDAC_oncopanel.example.maf"))
  expect_true(public_fixture$status != 0L, "Synthetic fixtures must remain eligible for Git tracking.")

  tracked <- command_status("git", "ls-files")
  expect_true(tracked$status == 0L, paste(tracked$output, collapse = "\n"))
  private_prefixes <- c("data/raw/", "data/processed/", "data/private/", "outputs/", "local_notes/")
  tracked_private <- tracked$output[vapply(tracked$output, function(path) {
    path == "config/local.R" || any(startsWith(path, private_prefixes))
  }, logical(1))]
  expect_true(!length(tracked_private), paste("Private paths are tracked:", paste(tracked_private, collapse = ", ")))
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

source("R/config.R")
source("R/data.R")
source("R/plot_helpers.R")

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

runtime_root <- file.path(
  project_root,
  "outputs",
  paste0("test_runtime_", Sys.getpid())
)
if (dir.exists(runtime_root)) unlink(runtime_root, recursive = TRUE, force = TRUE)
dir.create(runtime_root, recursive = TRUE, showWarnings = FALSE)

write_test_config <- function(filename, result_directory) {
  path <- file.path(runtime_root, filename)
  runtime_directory <- basename(runtime_root)
  lines <- c(
    "paad_config <- list(",
    '  clinical_file = "data/example/PDAC_clinical_info.example.xlsx",',
    '  target_patients_file = "data/example/result_527명.example.txt",',
    '  maf_file = "data/example/PDAC_oncopanel.example.maf",',
    sprintf(
      '  output_dir = file.path(PROJECT_ROOT, "outputs", "%s", "%s"),',
      runtime_directory,
      result_directory
    ),
    '  target_encoding = "CP949",',
    "  clinical_sheet = 1,",
    '  compatibility_mode = "v19",',
    '  sex_column_policy = "target_patients",',
    "  random_seed = 20250921L,",
    "  top_n = 20L,",
    "  plot_dpi = 150L",
    ")"
  )
  writeLines(lines, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

standalone_config <- write_test_config("config.standalone.R", "standalone")
run_all_config <- write_test_config("config.run_all.R", "run_all")
rscript <- file.path(R.home("bin"), "Rscript")

for (script in names(standalone_outputs)) {
  local({
    script_path <- script
    expected <- standalone_outputs[[script]]
    run_test(paste("standalone CLI", basename(script_path)), {
      result <- command_status(
        rscript,
        c("--vanilla", script_path, shQuote(paste0("--config=", standalone_config)))
      )
      expect_true(result$status == 0L, paste(result$output, collapse = "\n"))
      paths <- file.path(runtime_root, "standalone", expected)
      expect_true(all(file.exists(paths)), paste("Missing outputs:", paste(expected[!file.exists(paths)], collapse = ", ")))
      expect_true(all(file.info(paths)$size > 0), paste(basename(script_path), "created an empty output."))
    })
  })
}

run_test("standalone CLIs collectively create exactly the 42 documented outputs", {
  actual <- list.files(file.path(runtime_root, "standalone"), all.files = FALSE, no.. = TRUE)
  expect_true(
    setequal(actual, expected_outputs),
    paste(
      "Output mismatch.",
      "Missing:", paste(setdiff(expected_outputs, actual), collapse = ", "),
      "Unexpected:", paste(setdiff(actual, expected_outputs), collapse = ", ")
    )
  )
})

run_test("run_all completes and creates exactly the same 42 outputs", {
  result <- command_status(
    rscript,
    c("--vanilla", "run_all.R", shQuote(paste0("--config=", run_all_config)))
  )
  expect_true(result$status == 0L, paste(result$output, collapse = "\n"))
  result_dir <- file.path(runtime_root, "run_all")
  actual <- list.files(result_dir, all.files = FALSE, no.. = TRUE)
  expect_true(length(actual) == 42L, paste("Expected 42 files, found", length(actual)))
  expect_true(
    setequal(actual, expected_outputs),
    paste(
      "Output mismatch.",
      "Missing:", paste(setdiff(expected_outputs, actual), collapse = ", "),
      "Unexpected:", paste(setdiff(actual, expected_outputs), collapse = ", ")
    )
  )
  paths <- file.path(result_dir, expected_outputs)
  expect_true(all(file.info(paths)$size > 0), "run_all created one or more empty outputs.")
})

if (length(failures)) {
  cat("\nTest failures:\n", paste0("- ", failures, collapse = "\n"), "\n", sep = "")
  cat("Failed-run artifacts were kept under ", runtime_root, ".\n", sep = "")
  quit(status = 1L)
}

unlink(runtime_root, recursive = TRUE, force = TRUE)
cat("\nAll PAAD repository tests passed.\n")
