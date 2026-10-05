# Clinical extensions invoked explicitly from scripts/manuscript_analysis.R.
# All outputs are aggregates; patient-level resampling objects stay in private/.
run_clinical_extensions <- function(d,config,out,savefig,tab,pal,extended) {
  source("R/prognostic_validation.R")
  B<-config$validation_bootstraps
  base_formula<-reformulate(extended,response="survival::Surv(OS_m,survive)")
  formulas<-list(Clinical=base_formula,
    Clinical_KRAS=update(base_formula,.~.+KRAS_group),
    Clinical_driver_count=update(base_formula,.~.+driver_group),
    Clinical_genes=update(base_formula,.~.+KRAS_group+TP53+SMAD4+CDKN2A))
  needed<-unique(unlist(lapply(formulas,all.vars)))
  z<-d[complete.cases(d[,needed]),,drop=FALSE]
  model_labels<-c(Clinical="Clinical pathology",Clinical_KRAS="Clinical + KRAS subtype",
    Clinical_driver_count="Clinical + driver count",Clinical_genes="Clinical + KRAS and 3 genes")
  tab(data.frame(model=names(formulas),formula=vapply(formulas,function(f)paste(deparse(f),collapse=" "),character(1)),
    n=nrow(z),events=sum(z$survive),excluded_missing=nrow(d)-nrow(z),
    information_time="Pretreatment clinical/molecular measurements per investigator; surgical pathology used for postoperative prognosis"),"Input_Prediction_models.tsv")
  message("Clinical prediction comparison: ",nrow(z)," patients, ",sum(z$survive)," deaths")
  validation<-validate_survival_models(z,formulas,B=B)
  tab(validation$performance,"Supplementary_Table9_Prediction_performance.tsv")
  tab(validation$delta,"Supplementary_Table9_Prediction_differences.tsv")
  tab(validation$calibration,"Supplementary_Table9_OOB_calibration.tsv")
  tab(validation$failures,"Input_Validation_failures.tsv")
  tab(validation$coverage,"Input_Validation_coverage.tsv")
  # Global tests compare nested models on precisely the same complete cases.
  global<-do.call(rbind,lapply(2:4,function(k) {
    a<-validation$fits[[1]];b<-validation$fits[[k]];df<-length(coef(b))-length(coef(a));lr<-2*(b$loglik[2]-a$loglik[2])
    data.frame(model=names(formulas)[k],reference="Clinical",n=nrow(z),df=df,LR=lr,p=pchisq(lr,df,lower.tail=FALSE))
  }));global$p_holm<-p.adjust(global$p,"holm");tab(global,"Supplementary_Table9_Model_global_tests.tsv")
  perf<-validation$performance;perf$label<-factor(model_labels[perf$model],levels=rev(model_labels))
  cp<-ggplot2::ggplot(subset(perf,metric=="C_index"),ggplot2::aes(corrected,label))+
    ggplot2::geom_point(size=3,color=pal[2])+ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f",corrected)),hjust=-.3,size=4)+
    ggplot2::scale_x_continuous(limits=c(.5,.8),breaks=seq(.5,.8,.05))+
    ggplot2::labs(x="Optimism-corrected 36-month Harrell C-index (higher is better)",y=NULL,title="Additional discrimination from genomic information",
      caption=paste0(B," patient bootstraps; identical patients across all models. Internal validation only.\nPoint estimates shown; bootstrap OOB variability intervals are reported separately, not as corrected-value CIs."))
  bp<-ggplot2::ggplot(subset(perf,metric=="Brier_36"),ggplot2::aes(corrected,label))+
    ggplot2::geom_point(size=3,color=pal[4])+ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f",corrected)),hjust=-.3,size=4)+
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(.1,.25)))+
    ggplot2::labs(x="Optimism-corrected 36-month IPCW Brier score (lower is better)",y=NULL,title="Prediction error at 36 months",caption="Censoring-adjusted prediction error; an association does not necessarily improve prediction.")
  cal<-subset(validation$calibration,month==36);cal$label<-model_labels[cal$model]
  cap<-ggplot2::ggplot(cal,ggplot2::aes(predicted_survival,observed_survival))+
    ggplot2::geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lower95,ymax=upper95),width=.015,color=pal[2])+
    ggplot2::geom_point(size=2,color=pal[2])+ggplot2::facet_wrap(~label,ncol=2)+
    ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+ggplot2::labs(x="Predicted 36-month survival",y="Observed 36-month survival",
      title="Out-of-bag calibration",caption="Quintiles of bagged OOB predictions; bars are pointwise KM 95% CIs.\nEach patient's prediction uses only bootstrap fits that omitted that patient. Quintiles are not clinical risk thresholds.")
  savefig(cp,"Figure5A_Incremental_discrimination",11,4.5);savefig(bp,"Figure5B_Prediction_error",11,4.5)
  savefig(cap,"Figure5C_OOB_calibration",10,8)
  savefig(cowplot::plot_grid(cp,bp,cap,ncol=1,rel_heights=c(.25,.25,.5)),"Figure5_Added_prognostic_value",13,16)

  message("Standardized KRAS survival and absolute differences")
  std<-bootstrap_standardized_survival(z,formulas$Clinical_KRAS,B=B)
  tab(std$estimates,"Supplementary_Table10_Adjusted_survival.tsv")
  tab(std$differences,"Supplementary_Table10_Survival_differences.tsv")
  tab(std$curve,"Input_Standardized_survival_curves.tsv")
  sp<-ggplot2::ggplot(std$curve,ggplot2::aes(month,survival,color=group))+ggplot2::geom_line(linewidth=.9)+
    ggplot2::scale_color_manual(values=setNames(pal,levels(z$KRAS_group)))+
    ggplot2::scale_y_continuous(limits=c(0,1),labels=scales::percent)+
    ggplot2::labs(x="Months from surgery",y="Standardized overall survival",color=NULL,title="KRAS subtype survival after clinical pathology adjustment",
      caption="Predictions standardized to the same complete-case population. Observational associations, not causal effects.\n12-, 24- and 36-month estimates and pointwise bootstrap CIs are in Supplementary Table 10.")
  dd<-subset(std$differences,month==36);dd$contrast<-paste0(dd$group1,"\nReference: ",dd$reference_group)
  dp<-ggplot2::ggplot(dd,ggplot2::aes(survival_difference*100,reorder(contrast,survival_difference)))+
    ggplot2::geom_vline(xintercept=0,linetype=2,color="grey50")+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=lower95*100,xmax=upper95*100),orientation="y",width=.2)+
    ggplot2::geom_point(color=pal[2],size=2.5)+
    ggplot2::geom_text(ggplot2::aes(x=38,label=ifelse(p_holm<.001,"<0.001",sprintf("%.3f",p_holm))),size=3.6)+
    ggplot2::scale_x_continuous(limits=c(-25,43),breaks=seq(-20,30,10))+
    ggplot2::labs(x="36-month survival difference (percentage points)",y=NULL,
      title="Absolute survival differences between KRAS subtypes",subtitle="Right column: Holm p (all 30 pairs-by-time contrasts)",caption="Difference = first group minus reference. Positive: higher survival in the first group.\nPointwise bootstrap CIs; approximate p values use bootstrap SE. Full results: Supplementary Table 10.")
  savefig(sp,"Figure4C_Standardized_survival",12,7);savefig(dp,"Figure4D_Absolute_survival_difference",12,8)

  # Focused pathology effects; outcomes are prespecified, not data-mined cutoffs.
  set.seed(20260931L);pathres<-pathprobs<-list()
  for(spec in list(c("CDKN2A","Nodal involvement"),c("TP53","Poor differentiation"))) {
    g<-spec[1];label<-spec[2];dd<-d
    dd$outcome<-if(g=="CDKN2A")as.integer(as.character(dd$N_stage)!="0") else ifelse(is.na(dd$Differentiation),NA,as.integer(dd$Differentiation=="PD"))
    covs<-c("splines::ns(Age, df = 3)","Sex","T_stage","M_stage","Neoadjuvant")
    form<-reformulate(c(g,covs),response="outcome");dd<-dd[complete.cases(dd[,all.vars(form)]),,drop=FALSE]
    f<-glm(form,data=dd,family=binomial());ss<-coef(summary(f))[g,]
    stdprob<-function(fit,pop){sapply(0:1,function(value){nd<-pop;nd[[g]]<-value;mean(predict(fit,newdata=nd,type="response"))})}
    pt<-stdprob(f,dd);bs<-matrix(NA_real_,B,2)
    for(b in seq_len(B)) {
      bt<-dd[sample.int(nrow(dd),replace=TRUE),,drop=FALSE]
      ff<-tryCatch(suppressWarnings(glm(form,data=bt,family=binomial())),error=function(e)NULL)
      if(!is.null(ff)&&ff$converged&&all(is.finite(coef(ff))))bs[b,]<-stdprob(ff,bt)
    }
    diff<-bs[,2]-bs[,1]
    pathres[[g]]<-data.frame(gene=g,outcome=label,n=nrow(dd),events=sum(dd$outcome),comparison="Variant detected",reference="Not detected",
      OR=exp(ss[1]),lower95=exp(ss[1]-1.96*ss[2]),upper95=exp(ss[1]+1.96*ss[2]),p=ss[4],
      probability_difference=pt[2]-pt[1],difference_lower95=unname(quantile(diff,.025,na.rm=TRUE)),difference_upper95=unname(quantile(diff,.975,na.rm=TRUE)),successful=sum(complete.cases(bs)))
    pathprobs[[g]]<-data.frame(gene=g,outcome=label,status=c("Not detected","Variant detected"),probability=pt,
      lower95=apply(bs,2,quantile,.025,na.rm=TRUE),upper95=apply(bs,2,quantile,.975,na.rm=TRUE))
  }
  pa<-do.call(rbind,pathres);pa$p_holm<-p.adjust(pa$p,"holm");pa$definition<-"OR = outcome odds in variant-detected group / not-detected group; not a probability ratio"
  pp<-do.call(rbind,pathprobs);tab(pa,"Supplementary_Table8_Adjusted_pathology.tsv");tab(pp,"Supplementary_Table8_Pathology_probabilities.tsv")
  ann<-pa;ann$text<-sprintf("Detected / not detected\nOR %.2f [%.2f, %.2f]\nHolm p %s",pa$OR,pa$lower95,pa$upper95,format.pval(pa$p_holm,digits=2,eps=.001))
  pplot<-ggplot2::ggplot(pp,ggplot2::aes(status,probability,color=status))+
    ggplot2::geom_point(size=3)+ggplot2::geom_errorbar(ggplot2::aes(ymin=lower95,ymax=upper95),width=.15)+
    ggplot2::facet_wrap(~gene+outcome,nrow=1)+ggplot2::scale_color_manual(values=c("#747474",pal[2]))+
    ggplot2::geom_text(data=ann,ggplot2::aes(x=1.5,y=.96,label=text),inherit.aes=FALSE,color="black",size=3.7,vjust=1)+
    ggplot2::scale_y_continuous(limits=c(0,1),labels=scales::percent)+
    ggplot2::labs(x=NULL,y="Adjusted outcome probability",title="Adjusted gene-pathology associations",
      caption="Age spline, sex, T/M category and neoadjuvant treatment adjusted. Pointwise bootstrap CIs.\nNodal involvement: N1/2 vs N0. Poor differentiation: PD vs WD/MD. OR > 1: higher outcome odds when a variant is detected.")+
    ggplot2::theme(legend.position="none")
  savefig(pplot,"Figure2D_Adjusted_pathology",12,6)

  # Prespecified context checks. Never infer treatment benefit from these HRs.
  inter_formula<-update(formulas$Clinical_KRAS,.~.+KRAS_group:Neoadjuvant)
  fi<-safe_cox(inter_formula,z);fr<-validation$fits$Clinical_KRAS
  lr<-2*(fi$loglik[2]-fr$loglik[2]);df<-length(coef(fi))-length(coef(fr))
  tab(data.frame(interaction="KRAS subtype by neoadjuvant treatment",LR=lr,df=df,p=pchisq(lr,df,lower.tail=FALSE),n=nrow(z),
    interpretation="Prognostic heterogeneity among selected surgical patients; not a treatment-effect estimate"),"Supplementary_Table11_Interaction.tsv")
  subpairs<-subdiagnostics<-list()
  subsets<-list(All=z,M0=z[z$M_stage=="0",],No_neoadjuvant=z[z$Neoadjuvant=="n",],Neoadjuvant=z[z$Neoadjuvant=="y",])
  for(nm in names(subsets)) {
    zz<-droplevels(subsets[[nm]]);covs<-extended
    covs<-covs[vapply(covs,function(a)all(vapply(all.vars(reformulate(a)),function(v)length(unique(zz[[v]]))>1,logical(1))),logical(1))]
    f<-safe_cox(reformulate(c("KRAS_group",covs),response="survival::Surv(OS_m,survive)"),zz)
    if(is.null(f))stop("Sensitivity model failed: ",nm)
    pw<-pairwise_cox(f,zz,"KRAS_group");pw$cohort<-nm;pw$n<-nrow(zz);pw$events<-sum(zz$survive)
    pw$interpretation<-ratio_direction(pw$HR,pw$group1,pw$group2,"all-cause death",pw$p_holm<.05)
    subpairs[[nm]]<-pw
    ph<-survival::cox.zph(f)$table;subdiagnostics[[nm]]<-data.frame(cohort=nm,term=rownames(ph),p=ph[,"p"])
  }
  # Measured pretreatment CA19-9 and ASA are available, never adjuvant yes/no.
  ca<-d;ca$ASA<-factor(ca$ASA)
  if(any(ca$CA19_9<0,na.rm=TRUE))stop("Negative CA19-9 requires review")
  caformula<-update(formulas$Clinical_KRAS,.~.+log1p(CA19_9)+ASA)
  ca<-droplevels(ca[complete.cases(ca[,all.vars(caformula)]),,drop=FALSE])
  fc<-safe_cox(caformula,ca)
  if(!is.null(fc)) {
    pw<-pairwise_cox(fc,ca,"KRAS_group");pw$cohort<-"CA19_9_ASA_adjusted";pw$n<-nrow(ca);pw$events<-sum(ca$survive)
    pw$interpretation<-ratio_direction(pw$HR,pw$group1,pw$group2,"all-cause death",pw$p_holm<.05)
    subpairs[["CA19_9_ASA_adjusted"]]<-pw
  }
  sub<-do.call(rbind,subpairs);tab(sub,"Supplementary_Table11_Sensitivity_pairs.tsv")
  tab(do.call(rbind,subdiagnostics),"Supplementary_Table11_Sensitivity_PH.tsv")
  sub<-subset(sub,group1=="G12D" & group2 %in% c("G12V","G12R"))
  sfp<-ggplot2::ggplot(sub,ggplot2::aes(HR,cohort,color=group2))+
    ggplot2::geom_vline(xintercept=1,linetype=2,color="grey50")+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=lower95,xmax=upper95),orientation="y",width=.2,position=ggplot2::position_dodge(.4))+
    ggplot2::geom_point(position=ggplot2::position_dodge(.4),size=2)+ggplot2::scale_x_log10()+
    ggplot2::labs(x="Death HR: G12D / reference subtype",y=NULL,color="Reference subtype",title="KRAS subtype associations across clinical contexts",
      caption=paste0("HR > 1: higher death hazard in G12D. Pointwise CIs. Global KRAS-by-treatment interaction p ",format.pval(pchisq(lr,df,lower.tail=FALSE),digits=3),".\nAll 10 pairs per model and Holm p values are available. M0 is a sensitivity cohort, not the main cohort."))
  savefig(sfp,"Supplementary_Figure6_Clinical_context",12,7)

  # Time-varying exposure effects for endpoints with previous PH signals.
  tvrows<-list()
  for(endpoint in c("Recorded recurrence","Reported TMB")) {
    zz<-d
    if(endpoint=="Recorded recurrence") {
      zz$time<-zz$RFS_m;zz$event<-zz$Recur
      for(k in 2:5)zz[[paste0("K",k)]]<-as.integer(zz$KRAS_group==levels(zz$KRAS_group)[k])
      ex<-paste0("K",2:5);labels<-levels(zz$KRAS_group)[2:5];refs<-rep("Not detected",4)
    } else {zz$time<-zz$OS_m;zz$event<-zz$survive;ex<-"TMB_SD";labels<-"TMB higher by 1 SD";refs<-"Lower TMB"}
    zz$Neo_yes<-as.integer(zz$Neoadjuvant=="y")
    covs<-setdiff(extended,"Neoadjuvant");covs<-c(covs,"Neo_yes")
    # Permit the previously flagged treatment effect to vary as well.
    tv<-c(ex,"Neo_yes")
    if(endpoint=="Reported TMB") {
      for(k in levels(zz$T_stage)[-1]) {nm<-paste0("Tcat",k);zz[[nm]]<-as.integer(zz$T_stage==k);tv<-c(tv,nm);covs<-c(covs,nm)}
      covs<-setdiff(covs,"T_stage")
    }
    zz<-droplevels(zz[complete.cases(zz[,unique(c("time","event",ex,all.vars(reformulate(covs))))]),,drop=FALSE])
    covs<-covs[vapply(covs,function(a)all(vapply(all.vars(reformulate(a)),function(v)length(unique(zz[[v]]))>1,logical(1))),logical(1))]
    form<-reformulate(c(ex,covs,paste0("tt(",tv,")")),response="survival::Surv(time,event)")
    ft<-survival::coxph(form,data=zz,tt=function(x,t,...)x*log(pmax(t,1)/12),ties="efron",x=TRUE)
    for(k in seq_along(ex)) for(t in c(12,24,36)) {
      namesb<-c(ex[k],paste0("tt(",ex[k],")"));a<-c(1,log(t/12));beta<-sum(coef(ft)[namesb]*a)
      se<-sqrt(as.numeric(t(a)%*%vcov(ft)[namesb,namesb]%*%a))
      tvrows[[length(tvrows)+1L]]<-data.frame(endpoint=endpoint,month=t,comparison=labels[k],reference=refs[k],n=nrow(zz),events=sum(zz$event),
        HR=exp(beta),lower95=exp(beta-1.96*se),upper95=exp(beta+1.96*se),p=2*pnorm(-abs(beta/se)),
        method="Cox with log-time exposure interaction centered at 12 months; pointwise Wald CI")
    }
  }
  tvout<-do.call(rbind,tvrows);tvout$p_holm<-ave(tvout$p,tvout$endpoint,FUN=function(x)p.adjust(x,"holm"))
  tvout$interpretation<-ratio_direction(tvout$HR,tvout$comparison,tvout$reference,ifelse(tvout$endpoint=="Reported TMB","all-cause death","recorded recurrence"),tvout$p_holm<.05)
  tab(tvout,"Supplementary_Table12_Time_varying_HR.tsv")
  tvp<-ggplot2::ggplot(tvout,ggplot2::aes(month,HR,color=comparison,group=comparison))+
    ggplot2::geom_hline(yintercept=1,linetype=2,color="grey50")+ggplot2::geom_line()+ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lower95,ymax=upper95),width=.8)+ggplot2::facet_wrap(~endpoint,scales="free_y")+
    ggplot2::scale_y_log10()+ggplot2::labs(x="Months from surgery",y="Event hazard ratio",color="Comparison group",title="Time-varying associations",
      caption="Recurrence: each subtype / not detected. TMB: one SD higher / lower value, with death as the event.\nHR > 1: higher event hazard in the comparison group at that time. Pointwise CIs; Holm p values in the result table.")
  savefig(tvp,"Supplementary_Figure7_Time_varying_effects",13,7)
  # Missing TMB is not assigned a zero value or reconstructed from variant count.
  mis<-list()
  for(status in c("Observed","Missing")) {
    zz<-d[if(status=="Observed")!is.na(d$TMB) else is.na(d$TMB),]
    rkm<-survival::survfit(survival::Surv(OS_m,1-survive)~1,data=zz)
    mis[[status]]<-data.frame(TMB_status=status,n=nrow(zz),deaths=sum(zz$survive),median_age=median(zz$Age),
      neoadjuvant_n=sum(zz$Neoadjuvant=="y"),M1_n=sum(zz$M_stage=="1"),median_reverse_KM_followup=unname(summary(rkm)$table["median"]))
  }
  tab(do.call(rbind,mis),"Supplementary_Table13_TMB_missingness.tsv")
  saveRDS(list(models=validation$fits,standardized=std,performance=validation$performance),file.path(out,"private","clinical_extension_models.rds"))
  message("Clinical extensions completed")
}
