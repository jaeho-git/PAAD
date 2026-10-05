# Explicit numerator/reference metadata. HR concerns event hazard, not absolute
# probability, disease incidence, or a causal treatment effect.
ratio_direction <- function(ratio, numerator, reference, event="death", significant=NA) {
  mapply(function(r,a,b,e,s) {
    if(!is.finite(r)) return("Not estimable")
    direction <- if(r>=1) "higher" else "lower"
    change <- abs(r-1)*100
    paste0(a," has an estimated ",sprintf("%.1f",change),"% ",direction," ",e,
      " hazard than ",b," (reference)",if(is.na(s)) "." else if(s) "; multiplicity-adjusted evidence supports a difference." else "; a difference is not established after multiplicity adjustment.")
  },ratio,numerator,reference,event,significant,USE.NAMES=FALSE)
}
coefficient_direction <- function(x, data) {
  x$comparison_group <- x$reference_group <- NA_character_
  for(i in seq_len(nrow(x))) {
    ex <- x$exposure[i]; term <- x$term[i]
    if(!x$is_exposure[i]) next
    if(is.factor(data[[ex]])) {
      x$comparison_group[i] <- substring(term,nchar(ex)+1L)
      x$reference_group[i] <- levels(data[[ex]])[1]
    } else if(ex=="TMB_SD") {
      x$comparison_group[i] <- "Reported TMB higher by 1 SD";x$reference_group[i] <- "Lower value"
    } else {
      x$comparison_group[i] <- "Variant detected";x$reference_group[i] <- "Not detected"
    }
  }
  x$event_definition <- ifelse(x$endpoint=="OS","all-cause death","recorded recurrence")
  x$ratio_definition <- "Hazard in comparison_group / hazard in reference_group"
  x$interpretation <- ratio_direction(x$HR,x$comparison_group,x$reference_group,x$event_definition,x$q_BH<.05)
  x$interpretation[!x$is_exposure] <- "Covariate coefficient; spline basis coefficients are not standalone clinical effects."
  x
}
