# Clinical tracks for manuscript oncoplots. These are explicit display
# contracts: adding a variant-count bar must never silently drop a track.
source("R/annotation_palettes.R")
oncoplot_clinical_tracks <- function(profile = c("main", "repair")) {
  profile <- match.arg(profile)
  if (profile == "main") c("KRAS_group", "Sex", "Age_group", "Differentiation",
    "Neoadjuvant", "T_stage", "N_stage", "LN_positive", "M_stage", "BMI",
    "CA19_9", "CEA", "AJCC_stage", "Stage_Group", "Tumor_size", "TMB") else
    c("KRAS_group", "Age", "Sex", "Neoadjuvant", "Preop_platinum_exposure",
      "T_stage", "N_stage", "M_stage", "Differentiation", "R_status", "MSI", "TMB")
}

oncoplot_clinical_data <- function(clinical, sample_ids, profile = c("main", "repair")) {
  profile <- match.arg(profile)
  tracks <- oncoplot_clinical_tracks(profile)
  missing <- setdiff(c("Tumor_Sample_Barcode", tracks), names(clinical))
  if (length(missing)) stop("Required oncoplot clinical columns missing: ", paste(missing, collapse = ", "))
  if (anyDuplicated(clinical$Tumor_Sample_Barcode) || anyDuplicated(sample_ids)) stop("Oncoplot sample IDs must be unique.")
  index <- match(sample_ids, clinical$Tumor_Sample_Barcode)
  if (anyNA(index)) stop("Every oncoplot column must have a matched clinical record.")
  at <- as.data.frame(clinical[index, tracks, drop = FALSE])
  labels <- c(KRAS_group="KRAS subtype", Age="Age (years)", Age_group="Age group",
    Sex="Sex", Differentiation="Differentiation", Neoadjuvant="Neoadjuvant",
    Preop_platinum_exposure="Preop_platinum", T_stage="T category", N_stage="N category",
    LN_positive="LN positive", M_stage="M category", BMI="BMI", CA19_9="CA19-9",
    CEA="CEA", AJCC_stage="AJCC stage", Stage_Group="Stage group", Tumor_size="Tumor size",
    R_status="Margin", MSI="MSI", TMB="Reported TMB")
  at$Neoadjuvant <- ifelse(is.na(at$Neoadjuvant), NA_character_, ifelse(at$Neoadjuvant == "y", "Yes", "No"))
  if ("LN_positive" %in% tracks) at$LN_positive <- ifelse(is.na(at$LN_positive), NA_character_, ifelse(at$LN_positive == "1", "Positive", "Negative"))
  for (nm in c("T_stage", "N_stage", "M_stage")) at[[nm]] <- ifelse(is.na(at[[nm]]), NA_character_, paste0(substr(nm,1,1), at[[nm]]))
  colors <- list(); legend <- list(); audit <- list()
  fixed <- annotation_categorical_palettes()
  for (nm in tracks) {
    x <- at[[nm]]; observed <- sum(!is.na(x)); ceiling <- NA_real_; above <- 0L
    if (is.numeric(x)) {
      finite <- x[is.finite(x)]
      limits <- if (length(finite)) range(finite) else c(0,1)
      if (nm %in% c("CA19_9", "CEA") && length(finite)) {
        ceiling <- unname(stats::quantile(finite,.8)); limits[2] <- ceiling
        above <- sum(finite > ceiling)
      }
      if (limits[1] == limits[2]) limits <- limits + c(-1,1) * max(1,abs(limits[1])*.01)
      ramp <- annotation_continuous_palette(nm)
      colors[[nm]] <- circlize::colorRamp2(limits, ramp)
      breaks <- unique(c(limits[1],mean(limits),limits[2]))
      texts <- format(signif(breaks,3), trim=TRUE)
      if (!is.na(ceiling)) texts[length(texts)] <- paste0(">=", format(signif(ceiling,3),trim=TRUE))
      legend[[nm]] <- list(title=labels[[nm]],at=breaks,labels=texts)
    } else {
      levs <- sort(unique(as.character(stats::na.omit(x))))
      if (is.factor(x)) levs <- intersect(levels(x), levs)
      if (!length(levs)) levs <- "Not observed"
      # Clinical colors cannot depend on the different sample order of each plot.
      if(nm %in% names(fixed) && all(levs %in% names(fixed[[nm]]))) levs <- names(fixed[[nm]])[names(fixed[[nm]]) %in% levs]
      colors[[nm]] <- annotation_discrete_colors(nm,levs)
      legend[[nm]] <- list(title=labels[[nm]],at=levs,labels=levs)
    }
    audit[[nm]] <- data.frame(position=match(nm,tracks),column=nm,label=labels[[nm]],
      type=if(is.numeric(x)) "Continuous" else "Categorical", observed=observed,missing=sum(is.na(x)),
      categories=if(is.numeric(x)) NA_character_ else paste(names(colors[[nm]]),collapse="; "),
      color_ceiling=ceiling,values_above_color_ceiling=above,
      color_mapping=if(is.numeric(x)) paste(ramp,collapse=" -> ") else paste(paste(names(colors[[nm]]),colors[[nm]],sep="="),collapse="; "),
      missing_color="#BDBDBD")
  }
  list(data=at,colors=colors,legend=legend,labels=unname(labels[tracks]),audit=do.call(rbind,audit))
}
