# Project bootstrap helpers. This file has no package dependencies.

paad_project_root <- function(from = getwd()) {
  current <- normalizePath(from, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(current, "R", "bootstrap.R"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("PAAD project root could not be located.", call. = FALSE)
    current <- parent
  }
}

paad_script_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
    candidate <- dirname(script)
    if (basename(candidate) == "scripts") candidate <- dirname(candidate)
    return(paad_project_root(candidate))
  }
  paad_project_root()
}

source_project_modules <- function(project_root = paad_project_root()) {
  modules <- c(
    "config.R", "data.R", "plot_helpers.R", "plot_genomics.R",
    "driver_survival.R", "workflows.R"
  )
  for (module in modules) sys.source(file.path(project_root, "R", module), envir = .GlobalEnv)
  invisible(project_root)
}

parse_config_argument <- function(project_root, args = commandArgs(trailingOnly = TRUE)) {
  config_arg <- grep("^--config=", args, value = TRUE)
  if (length(config_arg) > 1L) stop("Specify --config only once.", call. = FALSE)
  path <- if (length(config_arg)) sub("^--config=", "", config_arg[[1]]) else "config/local.R"
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(project_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

run_entrypoint <- function(workflow) {
  project_root <- paad_script_root()
  source_project_modules(project_root)
  config_path <- parse_config_argument(project_root)
  config <- load_config(config_path, project_root)
  if (is.character(workflow)) workflow <- get(workflow, envir = .GlobalEnv, inherits = TRUE)
  workflow(config)
}
