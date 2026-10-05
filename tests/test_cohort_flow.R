source("R/cohort_flow.R")
flow <- data.frame(item=c("Source records","Duplicate-ID records excluded","Analysis patients","M0","M1","Duplicate patient IDs excluded"),
  n=c(1015,4,1011,972,39,2))
spec <- cohort_flow_spec(flow)
stopifnot(nrow(spec$nodes)==5, nrow(spec$edges)==6, sum(spec$edges$arrow)==4,
  grepl("1,015",spec$nodes$label[1]), grepl("4 (2 patients)",spec$nodes$label[2],fixed=TRUE),
  grepl("1,011",spec$nodes$label[3]), grepl("972",spec$nodes$label[4]), grepl("39",spec$nodes$label[5]))
plot <- plot_cohort_flow(spec)
stopifnot(inherits(plot,"ggplot"), identical(plot$labels$title,"Cohort selection"))
# Regression: connector segments without arrowheads must not disappear.
built <- ggplot2::ggplot_build(plot)
stopifnot(all(vapply(built$data[seq_len(nrow(spec$edges))],nrow,integer(1))==1))
bad <- flow; bad$n[3] <- 1012
stopifnot(inherits(try(cohort_flow_spec(bad),silent=TRUE),"try-error"))
cat("PASS: cohort flow counts, exclusions, stage inclusion and plotting.\n")
