# Fixed palettes: a variable/category keeps its colour across all oncoplots.
# Missing values use a separate grey colour, never the negative-category colour.
annotation_categorical_palettes <- function() list(
  KRAS_group=c(`Not detected`="#747474",G12D="#246A8A",G12V="#C17438",G12R="#438569",Other="#8865A3",`Other rare/multiple`="#D05A6E"),
  Sex=c(F="#E64B8A",M="#3B6FB6"),
  Age_group=c(`Age <=55`="#F3C623",`Age 56-69`="#00A6A6",`Age >=70`="#7A3E9D"),
  Differentiation=c(WD="#91CF60",MD="#FDD05A",PD="#B2182B"),
  Neoadjuvant=c(No="#F2D8A7",Yes="#24478F"),
  T_stage=c(T0="#FBE5D6",T1="#F4B183",T2="#CE5C17",T3="#6B2D0C",T4="#361304"),
  N_stage=c(N0="#DAC9EB",N1="#9B59B6",N2="#4B176C"),
  LN_positive=c(Negative="#B3DDE8",Positive="#006A7C"),
  M_stage=c(M0="#C5E1A5",M1="#7B3294"),
  AJCC_stage=c(IA="#A6CEE3",IB="#1F78B4",IIA="#B2DF8A",IIB="#33A02C",III="#FDBF6F",IV="#E31A1C"),
  Stage_Group=c(I="#F6E8C3",II="#DFC27D",III="#BF812D",IV="#543005"),
  Preop_platinum_exposure=c(No="#BFD3E6",Unknown="#E7BA52",Yes="#045A8D"),
  R_status=c(R0="#A8DDB5",R1="#E34A33"),
  MSI=c(`Probable stable`="#E7B800",`Probably stable`="#6846A5",Stable="#2B8CBE",Unstable="#CC1460"))

annotation_continuous_palette <- function(column) {
  if(column=="Size") column <- "Tumor_size"
  switch(column, Age=c("#DEEBF7","#08519C"), BMI=c("#E5F5E0","#237A36"),
    CA19_9=c("#FFF2CC","#B35806"), CEA=c("#EEE0F2","#762A83"),
    Tumor_size=c("#D9F0EF","#007F78"), TMB=c("#FCE0E8","#A5003B"),
    c("#E7E1EF","#54278F"))
}

annotation_discrete_colors <- function(column, levels) {
  aliases <- c(NAC="Neoadjuvant",T="T_stage",N="N_stage",M="M_stage",N_status="LN_positive",
    Stage="AJCC_stage",KRAS_clinical_group="KRAS_group",KRAS_MAF_group="KRAS_group")
  if(column %in% names(aliases)) column <- aliases[[column]]
  lookup <- levels
  if(column=="Neoadjuvant") lookup <- ifelse(lookup=="y","Yes",ifelse(lookup=="n","No",lookup))
  if(column=="LN_positive") lookup <- ifelse(lookup=="1","Positive",ifelse(lookup=="0","Negative",lookup))
  if(column %in% c("T_stage","N_stage","M_stage")) lookup <- ifelse(grepl("^[0-9]",lookup),paste0(substr(column,1,1),lookup),lookup)
  fixed <- annotation_categorical_palettes()[[column]]
  if(!is.null(fixed) && all(lookup %in% names(fixed))) return(stats::setNames(unname(fixed[lookup]),levels))
  # Stable fallback for legacy or new categories. It does not depend on patient
  # order or on which other annotation tracks happen to precede this one.
  hue <- sum(utf8ToInt(column)*seq_along(utf8ToInt(column))) %% 360
  lev <- sort(unique(levels))
  cols <- stats::setNames(grDevices::hcl(h=(hue+seq(0,300,length.out=length(lev)))%%360,c=65,l=60,fixup=TRUE),lev)
  cols[levels]
}
