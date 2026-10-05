source("R/prognostic_validation.R")
source("R/effect_direction.R")
source("R/pairwise_tests.R")
set.seed(7)
n<-240;d<-data.frame(Age=runif(n,35,85),group=factor(sample(c("A","B","C"),n,TRUE)))
death<-rexp(n,exp(.02*(d$Age-60)+.5*(d$group=="B"))/.03^-1);cens<-rexp(n,.012)
d$OS_m<-pmin(death,cens);d$survive<-as.integer(death<=cens)
f<-safe_cox(survival::Surv(OS_m,survive)~splines::ns(Age,df=3)+group,d)
stopifnot(!is.null(f))
ours<-survival_probabilities(f,d[1:4,],c(12,24,36))
native<-summary(survival::survfit(f,newdata=d[1:4,]),times=c(12,24,36),extend=TRUE)$surv
stopifnot(max(abs(ours-t(native)))<1e-8)
library(survival)
d$Differentiation<-factor(sample(c("WD","MD","PD"),n,TRUE))
fs<-safe_cox(Surv(OS_m,survive)~splines::ns(Age,df=3)+group+strata(Differentiation),d)
stopifnot(nrow(pairwise_cox(fs,d,"group"))==3)
for(i in 1:4) {
  ss<-summary(survfit(fs,newdata=d[i,,drop=FALSE]),times=c(12,24,36),extend=TRUE)$surv
  stopifnot(max(abs(survival_probabilities(fs,d[i,,drop=FALSE],c(12,24,36))-ss))<1e-8)
}
pw<-pairwise_cox(f,d,"group")
stopifnot(nrow(pw)==3,all(pw$p_holm>=pw$p))
ri<-ratio_direction(c(1.43,.7),c("A","A"),c("B","B"),c("death","recurrence"),c(TRUE,FALSE))
stopifnot(is.character(ri),length(ri)==2,grepl("43.0% higher death",ri[1]),grepl("30.0% lower recurrence",ri[2]))
# Brier score without censoring is ordinary mean squared survival error.
unc<-d;unc$survive<-1L
pr<-runif(n);bs<-ipcw_brier(unc$OS_m,unc$survive,pr,12,censor_survival(unc))
stopifnot(abs(bs-mean((as.integer(unc$OS_m>12)-pr)^2))<1e-10)
v<-validate_survival_models(d,list(Clinical=survival::Surv(OS_m,survive)~splines::ns(Age,df=3),
  Genomic=survival::Surv(OS_m,survive)~splines::ns(Age,df=3)+group),B=20)
stopifnot(all(v$performance$successful==20),all(v$coverage$min_oob_predictions>0),all(is.finite(v$performance$corrected)))
cat("PASS: stratified and unstratified predictions match survfit; spline contrasts; effect direction; IPCW identity; bootstrap refitting and OOB coverage\n")
