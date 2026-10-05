# Survival prediction utilities. No cutpoint search or outcome-driven selection.
# Each bootstrap refits the complete model, including spline knots. All models
# use the same resampled patients. Internal validation is not external validation.
safe_cox <- function(formula,data) {
  warned <- character()
  fit <- tryCatch(withCallingHandlers(survival::coxph(formula,data=data,x=TRUE,y=TRUE,model=TRUE,ties="efron"),
    warning=function(w){warned<<-c(warned,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)NULL)
  if(is.null(fit)||any(!is.finite(coef(fit)))||any(grepl("infinite|converg|singular",warned))) return(NULL)
  fit
}
survival_probabilities <- function(fit,newdata,times) {
  bh <- survival::basehaz(fit,centered=FALSE)
  lp <- as.numeric(predict(fit,newdata=newdata,type="lp",reference="zero"))
  if(!"strata" %in% names(bh)) {
    h <- c(0,bh$hazard)[findInterval(times,bh$time)+1L]
    return(exp(-outer(exp(lp),h,"*")))
  }
  # The manuscript's sole baseline-hazard stratum is pathological grade.
  # Non-PH grade is retained without forcing a constant grade coefficient.
  lev <- unique(as.character(bh$strata))
  grades <- sub("^.*=","",lev)
  if(!all(as.character(newdata$Differentiation) %in% grades)) stop("Unknown prediction stratum")
  ans <- matrix(NA_real_,nrow(newdata),length(times))
  for(k in seq_along(lev)) {
    rows <- as.character(newdata$Differentiation)==grades[k];b<-bh[as.character(bh$strata)==lev[k],]
    h <- c(0,b$hazard)[findInterval(times,b$time)+1L]
    ans[rows,] <- exp(-outer(exp(lp[rows]),h,"*"))
  }
  ans
}
censor_survival <- function(data) survival::survfit(survival::Surv(OS_m,1-survive)~1,data=data)
censor_at <- function(fit,times,left=FALSE) {
  # G(t-) at deaths, G(t) for event-free controls.
  ii <- findInterval(times,fit$time)
  if(left) ii <- ii-as.integer(ii>0L & fit$time[pmax(ii,1L)]==times)
  c(1,fit$surv)[ii+1L]
}
ipcw_brier <- function(time,event,surv_prob,horizon,censor_fit) {
  gd <- censor_at(censor_fit,time,TRUE); gt <- censor_at(censor_fit,horizon)
  cases <- event==1 & time<=horizon; controls <- time>horizon
  if(gt<.05 || any(gd[cases]<.05)) return(NA_real_)
  err <- numeric(length(time))
  err[cases] <- surv_prob[cases]^2/gd[cases]
  err[controls] <- (1-surv_prob[controls])^2/gt
  mean(err)
}
prediction_metrics <- function(fit,data,times=c(12,24,36)) {
  lp <- as.numeric(predict(fit,newdata=data,type="lp",reference="sample"))
  pred <- survival_probabilities(fit,data,times)
  risk36 <- 1-survival_probabilities(fit,data,36)[,1]
  cc <- survival::concordance(survival::Surv(pmin(data$OS_m,36),data$survive*(data$OS_m<=36))~risk36,reverse=TRUE)$concordance
  slope <- tryCatch(if(!is.null(fit$strata)) unname(coef(survival::coxph(survival::Surv(data$OS_m,data$survive)~lp+strata(data$Differentiation)))) else
    unname(coef(survival::coxph(survival::Surv(data$OS_m,data$survive)~lp))),error=function(e)NA_real_)
  g <- censor_survival(data)
  bs <- vapply(seq_along(times),function(j)ipcw_brier(data$OS_m,data$survive,pred[,j],times[j],g),numeric(1))
  setNames(c(cc,slope,bs),c("C_index","Calibration_slope",paste0("Brier_",times)))
}
validate_survival_models <- function(data,formulas,B=500L,times=c(12,24,36),seed=20260929L) {
  set.seed(seed); n<-nrow(data); names_m<-names(formulas)
  fitted<-lapply(formulas,safe_cox,data=data)
  if(any(vapply(fitted,is.null,logical(1)))) stop("Primary prediction model not estimable")
  apparent<-sapply(fitted,prediction_metrics,data=data,times=times)
  metrics<-rownames(apparent)
  optimism<-array(NA_real_,c(B,length(metrics),length(formulas)),dimnames=list(NULL,metrics,names_m))
  oob_metric<-optimism
  oob_sum<-lapply(formulas,function(x)matrix(0,n,length(times)))
  oob_n<-matrix(0,n,length(formulas),dimnames=list(NULL,names_m))
  failures<-list()
  for(b in seq_len(B)) {
    idx<-sample.int(n,n,replace=TRUE);train<-data[idx,,drop=FALSE];oob<-setdiff(seq_len(n),unique(idx))
    for(k in seq_along(formulas)) {
      f<-safe_cox(formulas[[k]],train)
      if(is.null(f)){failures[[length(failures)+1L]]<-data.frame(iteration=b,model=names_m[k],reason="Fit failed or unstable");next}
      a<-tryCatch(prediction_metrics(f,train,times),error=function(e)NULL)
      v<-tryCatch(prediction_metrics(f,data,times),error=function(e)NULL)
      if(is.null(a)||is.null(v)){failures[[length(failures)+1L]]<-data.frame(iteration=b,model=names_m[k],reason="Evaluation failed");next}
      optimism[b,,k]<-a-v
      if(length(oob)>50) {
        oo<-tryCatch(prediction_metrics(f,data[oob,,drop=FALSE],times),error=function(e)NULL)
        if(!is.null(oo)) {
          oob_metric[b,,k]<-oo
          oob_sum[[k]][oob,]<-oob_sum[[k]][oob,]+survival_probabilities(f,data[oob,,drop=FALSE],times)
          oob_n[oob,k]<-oob_n[oob,k]+1
        }
      }
    }
    if(b%%50L==0L) message("Prediction validation bootstrap ",b,"/",B)
  }
  results<-do.call(rbind,lapply(seq_along(formulas),function(k) do.call(rbind,lapply(seq_along(metrics),function(j) {
    opt<-optimism[,j,k];oo<-oob_metric[,j,k]
    data.frame(model=names_m[k],metric=metrics[j],n=n,events=sum(data$survive),apparent=apparent[j,k],
      optimism=mean(opt,na.rm=TRUE),corrected=apparent[j,k]-mean(opt,na.rm=TRUE),
      successful=sum(is.finite(opt)),requested=B,oob_median=median(oo,na.rm=TRUE),
      oob_p025=unname(quantile(oo,.025,na.rm=TRUE)),oob_p975=unname(quantile(oo,.975,na.rm=TRUE)))
  }))))
  delta<-do.call(rbind,lapply(2:length(formulas),function(k) do.call(rbind,lapply(seq_along(metrics),function(j) {
    oo<-oob_metric[,j,k]-oob_metric[,j,1]
    a<-apparent[j,k]-apparent[j,1]
    opt<-optimism[,j,k]-optimism[,j,1]
    data.frame(model=names_m[k],reference=names_m[1],metric=metrics[j],apparent_delta=a,
      corrected_delta=a-mean(opt,na.rm=TRUE),oob_median_delta=median(oo,na.rm=TRUE),
      oob_p025=unname(quantile(oo,.025,na.rm=TRUE)),oob_p975=unname(quantile(oo,.975,na.rm=TRUE)),
      paired_iterations=sum(is.finite(oo)),interval_definition="Central 95% OOB bootstrap distribution; not a confidence interval for corrected delta")
  }))))
  calibration<-list()
  for(k in seq_along(formulas)) for(j in seq_along(times)) {
    pred<-oob_sum[[k]][,j]/oob_n[,k]
    if(any(!is.finite(pred))) stop("Insufficient OOB prediction coverage")
    bin<-pmin(5L,ceiling(rank(pred,ties.method="first")/n*5))
    for(q in 1:5) {
      zz<-data[bin==q,,drop=FALSE];sf<-survival::survfit(survival::Surv(OS_m,survive)~1,data=zz)
      ss<-summary(sf,times=times[j],extend=TRUE)
      calibration[[length(calibration)+1L]]<-data.frame(model=names_m[k],month=times[j],bin=q,n=nrow(zz),
        n_at_risk=sum(zz$OS_m>=times[j]),predicted_survival=mean(pred[bin==q]),observed_survival=ss$surv,
        lower95=ss$lower,upper95=ss$upper,method="KM within quintiles of bagged out-of-bag survival predictions")
    }
  }
  list(fits=fitted,performance=results,delta=delta,calibration=do.call(rbind,calibration),
    failures=if(length(failures))do.call(rbind,failures) else data.frame(iteration=integer(),model=character(),reason=character()),
    coverage=data.frame(model=names_m,min_oob_predictions=apply(oob_n,2,min),median_oob_predictions=apply(oob_n,2,median)))
}
standardized_survival <- function(fit,data,variable,times) {
  lev<-levels(data[[variable]])
  sapply(lev,function(g){nd<-data;nd[[variable]]<-factor(g,levels=lev);colMeans(survival_probabilities(fit,nd,times))})
}
bootstrap_standardized_survival <- function(data,formula,variable="KRAS_group",B=500L,times=c(12,24,36),seed=20260930L) {
  set.seed(seed);f<-safe_cox(formula,data);lev<-levels(data[[variable]])
  point<-standardized_survival(f,data,variable,times)
  samples<-array(NA_real_,c(B,length(times),length(lev)))
  for(b in seq_len(B)) {
    z<-data[sample.int(nrow(data),replace=TRUE),,drop=FALSE];fb<-safe_cox(formula,z)
    if(!is.null(fb)) samples[b,,]<-tryCatch(standardized_survival(fb,z,variable,times),error=function(e)NA_real_)
  }
  est<-do.call(rbind,lapply(seq_along(lev),function(k) do.call(rbind,lapply(seq_along(times),function(j)
    data.frame(group=lev[k],month=times[j],survival=point[j,k],lower95=quantile(samples[,j,k],.025,na.rm=TRUE),
      upper95=quantile(samples[,j,k],.975,na.rm=TRUE),successful=sum(is.finite(samples[,j,k])),reference_population_n=nrow(data))))))
  pairs<-t(combn(seq_along(lev),2));res<-list()
  for(i in seq_len(nrow(pairs))) for(j in seq_along(times)) {
    a<-pairs[i,1];b<-pairs[i,2];dd<-samples[,j,a]-samples[,j,b];delta<-point[j,a]-point[j,b];se<-sd(dd,na.rm=TRUE)
    res[[length(res)+1L]]<-data.frame(group1=lev[a],reference_group=lev[b],month=times[j],survival_difference=delta,
      lower95=unname(quantile(dd,.025,na.rm=TRUE)),upper95=unname(quantile(dd,.975,na.rm=TRUE)),
      p=2*pnorm(-abs(delta/se)),definition="Standardized survival probability in group1 minus reference; positive favors group1")
  }
  dif<-do.call(rbind,res);dif$p_holm<-p.adjust(dif$p,"holm")
  dif$method<-"Bootstrap pointwise percentile CI; approximate normal test using bootstrap SE; Holm across 30 differences"
  curves<-standardized_survival(f,data,variable,0:36)
  curve<-data.frame(month=rep(0:36,length(lev)),group=rep(lev,each=37),survival=as.vector(curves))
  list(estimates=est,differences=dif,curve=curve)
}
