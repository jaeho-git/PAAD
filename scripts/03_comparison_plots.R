#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "..", "R", "bootstrap.R"))
run_entrypoint("workflow_comparison_plots")
