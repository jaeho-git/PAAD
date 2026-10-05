# Forest plots with a separate, row-aligned statistics panel. No HR/p/q strings
# are embedded in axis labels. All numbers are existing model outputs.
forest_display_data <- function(x, labels=NULL, p_column="p", p_heading="p") {
  if(is.null(labels)) labels <- x$term
  stopifnot(nrow(x)==length(labels),p_column %in% names(x),all(c("HR","lower95","upper95") %in% names(x)))
  fmt_p <- function(p) ifelse(is.na(p),"NE",ifelse(p<.001,"<0.001",sprintf("%.3f",p)))
  out <- data.frame(comparison=as.character(labels),row=rev(seq_len(nrow(x))),
    HR=x$HR,lower95=x$lower95,upper95=x$upper95,p=x[[p_column]])
  out$hr_text <- ifelse(is.finite(x$HR)&is.finite(x$lower95)&is.finite(x$upper95),
    sprintf("%.2f [%.2f, %.2f]",x$HR,x$lower95,x$upper95),"Not estimable")
  out$p_text <- fmt_p(out$p)
  if("q_BH" %in% names(x)) {out$q_BH <- x$q_BH;out$q_text <- fmt_p(x$q_BH)}
  attr(out,"p_heading") <- p_heading
  out
}

forest_result_plot <- function(x,title,labels=NULL,p_column="p",p_heading="p") {
  a <- forest_display_data(x,labels,p_column,p_heading)
  n <- nrow(a); has_q <- "q_text" %in% names(a); header <- n+.8
  a$comparison <- vapply(a$comparison,function(s)paste(strwrap(s,32),collapse="\n"),character(1))
  base <- ggplot2::ggplot(a,ggplot2::aes(y=row)) +
    ggplot2::scale_y_continuous(limits=c(.4,n+1.3),expand=c(0,0),breaks=NULL) +
    ggplot2::labs(y=NULL,x=NULL) + ggplot2::theme_classic(base_family="Pretendard",base_size=12) +
    ggplot2::theme(axis.line.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),plot.margin=ggplot2::margin(4,8,4,8))
  band <- a[a$row%%2==0,,drop=FALSE]
  base <- base + ggplot2::geom_rect(data=band,ggplot2::aes(ymin=row-.48,ymax=row+.48),
    xmin=-Inf,xmax=Inf,fill="#F5F5F5",colour=NA,inherit.aes=FALSE)
  table_style <- ggplot2::theme(axis.line.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank())
  left <- base + ggplot2::scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
    ggplot2::geom_text(ggplot2::aes(x=.01,label=comparison),hjust=0,size=3.7,family="Pretendard",lineheight=.95) +
    ggplot2::annotate("text",x=.01,y=header,label="Comparison",hjust=0,fontface="bold",size=3.9,family="Pretendard") + table_style
  middle <- base + ggplot2::geom_vline(xintercept=1,linetype=2,colour="#888888") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin=lower95,xmax=upper95,x=HR),orientation="y",width=.15,na.rm=TRUE) +
    ggplot2::geom_point(ggplot2::aes(x=HR),size=2.7,colour="#246A8A",na.rm=TRUE) +
    ggplot2::scale_x_log10(breaks=scales::breaks_log(n=4)) + ggplot2::labs(x="Hazard ratio (log scale)")
  p_pos <- if(has_q) .66 else .86
  right <- base + ggplot2::scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
    ggplot2::geom_text(ggplot2::aes(x=.02,label=hr_text),hjust=0,size=3.7,family="Pretendard") +
    ggplot2::geom_text(ggplot2::aes(x=p_pos,label=p_text),hjust=.5,size=3.7,family="Pretendard") +
    ggplot2::annotate("text",x=.02,y=header,label="HR [95% CI]",hjust=0,fontface="bold",size=3.9,family="Pretendard") +
    ggplot2::annotate("text",x=p_pos,y=header,label=p_heading,hjust=.5,fontface="bold",size=3.9,family="Pretendard") + table_style
  if(has_q) right <- right + ggplot2::geom_text(ggplot2::aes(x=.91,label=q_text),hjust=.5,size=3.7,family="Pretendard") +
    ggplot2::annotate("text",x=.91,y=header,label="q (BH)",hjust=.5,fontface="bold",size=3.9,family="Pretendard")
  plot <- cowplot::ggdraw(cowplot::plot_grid(left,middle,right,nrow=1,align="h",axis="tb",rel_widths=c(.30,.32,.38))) +
    ggplot2::labs(title=title,caption="HR = comparison-group event hazard / reference-group event hazard.\nHR > 1: higher hazard in the comparison group; HR < 1: lower hazard. Not an absolute probability ratio.") + ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=15,family="Pretendard"),
      plot.caption=ggplot2::element_text(hjust=0,size=10,family="Pretendard"),plot.margin=ggplot2::margin(8,12,8,12))
  attr(plot,"forest_display") <- a
  plot
}
