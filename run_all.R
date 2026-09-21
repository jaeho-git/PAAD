#!/usr/bin/env Rscript
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_path), "R", "bootstrap.R"))
run_entrypoint("workflow_all")
