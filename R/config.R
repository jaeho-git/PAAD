# Shared configuration helpers.
#
# Analysis code lives in scripts/*.R. This file only reads config files,
# validates input/output paths, and keeps generated files inside outputs/.

load_config <- function(path, project_root = getwd()) {
  if (!file.exists(path)) {
    stop(
      "Configuration file not found: ", path,
      "\nCopy config/config.example.R to config/local.R and edit local paths.",
      call. = FALSE
    )
  }
  config_env <- new.env(parent = baseenv())
  config_env$PROJECT_ROOT <- project_root
  sys.source(path, envir = config_env)
  if (!exists("paad_config", envir = config_env, inherits = FALSE)) {
    stop("The config file must define a list named paad_config.", call. = FALSE)
  }
  config <- get("paad_config", envir = config_env, inherits = FALSE)
  validate_config(config, project_root)
}

load_config_from_command_line <- function(project_root = getwd(), args = commandArgs(trailingOnly = TRUE)) {
  config_arguments <- grep("^--config=", args, value = TRUE)
  if (length(config_arguments) > 1L) {
    stop("Specify --config only once.", call. = FALSE)
  }

  config_file <- if (length(config_arguments)) {
    sub("^--config=", "", config_arguments[[1]])
  } else {
    "config/local.R"
  }
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", config_file)) {
    config_file <- file.path(project_root, config_file)
  }

  load_config(config_file, project_root)
}

normalize_config_path <- function(path, project_root, must_work = FALSE) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("Each configured path must be one non-empty character value.", call. = FALSE)
  }
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(project_root, path)
  normalizePath(path, winslash = "/", mustWork = must_work)
}

path_is_within <- function(path, parent) {
  path <- paste0(normalizePath(path, winslash = "/", mustWork = FALSE), "/")
  parent <- paste0(normalizePath(parent, winslash = "/", mustWork = FALSE), "/")
  startsWith(path, parent)
}

validate_config <- function(config, project_root) {
  if (!is.list(config)) stop("paad_config must be a list.", call. = FALSE)
  required <- c("clinical_file", "target_patients_file", "maf_file", "output_dir")
  missing <- setdiff(required, names(config))
  if (length(missing)) stop("Missing config fields: ", paste(missing, collapse = ", "), call. = FALSE)

  for (field in c("clinical_file", "target_patients_file", "maf_file")) {
    config[[field]] <- normalize_config_path(config[[field]], project_root, must_work = TRUE)
    if (dir.exists(config[[field]])) stop(field, " must point to a file.", call. = FALSE)
  }

  config$output_dir <- normalize_config_path(config$output_dir, project_root, must_work = FALSE)
  allowed_root <- file.path(project_root, "outputs")
  if (!path_is_within(config$output_dir, allowed_root) &&
      !identical(config$output_dir, normalizePath(allowed_root, winslash = "/", mustWork = FALSE))) {
    stop("output_dir must be outputs/ or one of its subdirectories inside this repository.", call. = FALSE)
  }
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  defaults <- list(
    target_encoding = "CP949",
    clinical_sheet = 1,
    compatibility_mode = "v19",
    sex_column_policy = "legacy",
    random_seed = NULL,
    top_n = 20L,
    plot_dpi = 1000L
  )
  for (field in names(defaults)) if (is.null(config[[field]])) config[[field]] <- defaults[[field]]
  if (!identical(config$compatibility_mode, "v19")) {
    stop("Only compatibility_mode = 'v19' is currently implemented.", call. = FALSE)
  }
  allowed_sex <- c("legacy", "target_patients", "clinical_korean")
  if (!config$sex_column_policy %in% allowed_sex) {
    stop("sex_column_policy must be one of: ", paste(allowed_sex, collapse = ", "), call. = FALSE)
  }
  if (!is.null(config$random_seed)) {
    config$random_seed <- as.integer(config$random_seed)
    if (is.na(config$random_seed)) stop("random_seed must be NULL or an integer.", call. = FALSE)
  }
  config$plot_dpi <- as.integer(config$plot_dpi)
  if (is.na(config$plot_dpi) || config$plot_dpi < 72L) stop("plot_dpi must be an integer of at least 72.", call. = FALSE)
  config$top_n <- as.integer(config$top_n)
  if (is.na(config$top_n) || config$top_n < 1L) stop("top_n must be a positive integer.", call. = FALSE)
  config$project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  class(config) <- c("paad_config", "list")
  config
}

output_path <- function(config, filename) {
  if (grepl("[/\\\\]", filename)) stop("Output filename must not contain a directory separator.", call. = FALSE)
  path <- file.path(config$output_dir, filename)
  if (!path_is_within(path, config$output_dir)) stop("Unsafe output path.", call. = FALSE)
  path
}
