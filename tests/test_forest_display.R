source("R/forest_display.R")
cowplot::set_null_device(function(width,height)ragg::agg_png(tempfile(fileext=".png"),width=width,height=height,units="in",res=96))
x <- data.frame(term=c("KRAS","TP53","Unavailable"),HR=c(1.2,.8,NA),lower95=c(1.01,.6,NA),
  upper95=c(1.5,1.1,NA),p=c(.0001,.08,NA),q_BH=c(.0012,.16,NA),p_holm=c(.0003,.24,NA))
a <- forest_display_data(x)
stopifnot(identical(a$comparison,x$term),identical(a$p,x$p),identical(a$q_BH,x$q_BH),
  identical(a$p_text,c("<0.001","0.080","NE")),a$hr_text[3]=="Not estimable")
pair <- forest_display_data(x,p_column="p_holm",p_heading="Holm p")
stopifnot(identical(pair$p,x$p_holm),attr(pair,"p_heading")=="Holm p",identical(x$p,c(.0001,.08,NA_real_)))
p <- forest_result_plot(x[1:2,],"Test forest")
stopifnot(inherits(p,"ggplot"),identical(attr(p,"forest_display")$comparison,x$term[1:2]),
  is.null(p$labels$y),is.null(p$labels$x))
cat("PASS: statistics separate from comparison labels; unchanged raw/BH/Holm values; missing estimates explicit.\n")
