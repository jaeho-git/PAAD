# Figure lettering belongs to editable presentation/document labels, not pixels.
is_panel_title <- function(x) {
  is.character(x) && length(x) == 1L &&
    grepl("^([a-z]|[0-9]+[a-z])[[:space:].:)]+", x)
}
check_call <- function(expr) {
  if (is.call(expr)) {
    fn <- paste(deparse(expr[[1L]]), collapse = "")
    args <- as.list(expr)[-1L]
    named <- names(args)
    if (!is.null(named) && "title" %in% named) {
      stopifnot(!is_panel_title(args[["title"]]))
    }
    if (fn == "draw_onco" && length(args) >= 3L) {
      stopifnot(!is_panel_title(args[[3L]]))
    }
    if (fn == "cowplot::plot_grid") {
      stopifnot(!"labels" %in% named)
    }
    invisible(lapply(args, check_call))
  } else if (is.expression(expr) || is.pairlist(expr)) {
    invisible(lapply(expr, check_call))
  }
}
for (file in c("R/cohort_flow.R", "scripts/manuscript_analysis.R")) {
  check_call(parse(file))
}
cat("PASS: figure titles have no panel prefixes or embedded panel letters.\n")
