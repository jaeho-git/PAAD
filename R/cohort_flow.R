# One aggregate specification drives the R figure and editable PPTX flowchart.
cohort_flow_spec <- function(flow) {
  n <- stats::setNames(flow$n, flow$item)
  stopifnot(n[["Source records"]] - n[["Duplicate-ID records excluded"]] == n[["Analysis patients"]],
    n[["M0"]] + n[["M1"]] == n[["Analysis patients"]])
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  nodes <- data.frame(
    id = c("source", "excluded", "analysis", "m0", "m1"),
    x = c(265, 655, 265, 110, 510), y = c(20, 145, 250, 400, 400),
    w = c(350, 335, 350, 270, 270), h = c(85, 105, 90, 75, 75),
    label = c(paste0("Updated clinical records\nn = ", fmt(n[["Source records"]])),
      paste0("Excluded duplicate-ID records\nn = ", fmt(n[["Duplicate-ID records excluded"]]),
        " (", fmt(n[["Duplicate patient IDs excluded"]]), " patients)\nAll records for these IDs excluded"),
      paste0("Analysis cohort\n", fmt(n[["Analysis patients"]]), " unique patients"),
      paste0("M0 included\nn = ", fmt(n[["M0"]])),
      paste0("M1 included\nn = ", fmt(n[["M1"]]))),
    fill = c("#FFFFFF", "#FFFFFF", "#EAF2F6", "#FFFFFF", "#FFFFFF"),
    font_size = c(22, 18, 23, 22, 22))
  edges <- data.frame(x = c(440,440,440,440,245,645), y = c(105,197.5,340,370,370,370),
    xend = c(440,655,440,645,245,645), yend = c(250,197.5,370,370,400,400),
    arrow = c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE))
  # The horizontal branch extends left as well as right from the cohort.
  edges$x[4] <- 245
  list(width = 1000, height = 500, nodes = nodes, edges = edges)
}

plot_cohort_flow <- function(spec) {
  nodes <- spec$nodes; edges <- spec$edges
  p <- ggplot2::ggplot() + ggplot2::theme_void(base_family = "Pretendard") +
    ggplot2::coord_fixed(xlim = c(0,spec$width), ylim = c(-spec$height,0), expand = FALSE) +
    ggplot2::labs(title = "Cohort selection") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",size = 16),
      plot.margin = ggplot2::margin(10,10,10,10))
  for (i in seq_len(nrow(edges))) {
    e <- edges[i,]
    args <- list("segment", x = e$x, y = -e$y, xend = e$xend, yend = -e$yend,
      linewidth = .5, colour = "#555555")
    # Supplying arrow=NULL to annotate() drops non-arrow segments in ggplot2.
    if(e$arrow) args$arrow <- grid::arrow(length = grid::unit(2.2,"mm"),type = "closed")
    p <- p + do.call(ggplot2::annotate,args)
  }
  for (i in seq_len(nrow(nodes))) {
    z <- nodes[i,]
    p <- p + ggplot2::annotate("rect",xmin = z$x,xmax = z$x+z$w,ymin = -z$y-z$h,ymax = -z$y,
      fill = z$fill,colour = "#555555",linewidth = .5) +
      ggplot2::annotate("text",x = z$x+z$w/2,y = -z$y-z$h/2,label = z$label,
        family = "Pretendard",size = z$font_size/4,lineheight = 1.1)
  }
  p
}
