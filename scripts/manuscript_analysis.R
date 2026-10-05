#!/usr/bin/env Rscript
# Publication analyses in reading order. Run from the repository root:
# Rscript --vanilla scripts/manuscript_analysis.R --config=config/local.R
# Figures 1-5: cohort, pathology, gene-level OS, KRAS OS and added prognostic value.
# Supplements: repair-gene profiles, recurrence, reported TMB, diagnostics.
# No clinical/MAF concordance or source-comparison outputs are produced here.
source("R/config.R")
source("R/data.R")
source("R/genomics.R")
source("R/repair_gene_sets.R")
source("R/pairwise_tests.R")
source("R/plot_helpers.R")
source("R/oncoplot_counts.R")
source("R/oncoplot_annotations.R")
source("R/cohort_flow.R")
source("R/forest_display.R")
source("R/effect_direction.R")
source("R/manuscript_story_extensions.R")
source("R/recurrence_signature_extensions.R")
config <- load_config_from_command_line()
if (!identical(config$clinical_schema, "curated_v3")) stop("The manuscript analysis requires curated_v3.")
if (!identical(config$updated_cohort, "all_unique")) stop("Manuscript cohort must include all stages, excluding duplicate-ID rows only.")
if (is.null(config$kras_source)) stop("Choose kras_source = 'clinical', 'maf', or 'workbook_maf' explicitly in the local config before publication analysis.")
assert_packages(c("cowplot", "ragg", "jsonlite", "png", "broom"))
out <- file.path(config$output_dir, config$manuscript_subdir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), showWarnings = FALSE)
dir.create(file.path(out, "tables"), showWarnings = FALSE)
dir.create(file.path(out, "private"), showWarnings = FALSE)
options(device = function(...) ragg::agg_png(file.path(out, "private", "layout_device.png"), width = 1600, height = 1200, res = 120))
cowplot::set_null_device(function(width, height) ragg::agg_png(tempfile(fileext=".png"), width=width, height=height, units="in", res=96))
input <- load_analysis_data(config)
set.seed(if (is.null(config$random_seed)) 20260928L else config$random_seed)
raw_variants <- as.data.frame(input$maf@data)
ids <- as.character(input$clinical$Tumor_Sample_Barcode)
vc <- oncoplot_variant_counts(raw_variants, ids, unique(raw_variants$Hugo_Symbol))
v <- vc$gene_events
d <- add_source_genomics(input$clinical, v)
d$Variant_count <- vc$patients$variant_count[match(d$Tumor_Sample_Barcode, vc$patients$sample_id)]
d$KRAS_group <- switch(config$kras_source,
  clinical = {
    x <- d$KRAS_clinical_group
    levels(x) <- c("Not detected", "G12D", "G12V", "G12R", "Other rare/multiple")
    x
  },
  maf = {
    x <- d$KRAS_MAF_group
    levels(x) <- c("Not detected", "G12D", "G12V", "G12R", "Other rare/multiple")
    x
  },
  # For the curated workbook source, preserve the exact analysis labels in
  # Main!KRAS_subtype_MAF rather than exposing internal normalization labels.
  workbook_maf = factor(as.character(d$KRAS_subtype_MAF),
    levels = c("Not detected", "G12D", "G12V", "G12R", "Other")),
  stop("Unsupported kras_source: ", config$kras_source))
genes <- c("KRAS", "TP53", "SMAD4", "CDKN2A")
for (g in genes) d[[g]] <- as.integer(d$Tumor_Sample_Barcode %in% v$Tumor_Sample_Barcode[v$Hugo_Symbol == g])
# One canonical KRAS definition also supplies binary status and driver count.
d$KRAS <- ifelse(is.na(d$KRAS_group), NA_integer_, as.integer(d$KRAS_group != "Not detected"))
d$driver_count <- rowSums(d[, genes])
d$driver_group <- factor(ifelse(d$driver_count >= 3, "3+", as.character(d$driver_count)), levels = c("1", "0", "2", "3+"))
d$Any_driver <- as.integer(d$driver_count > 0)
d$KRAS_TP53 <- factor(ifelse(d$KRAS + d$TP53 == 0, "Neither", ifelse(d$KRAS + d$TP53 == 1, "One gene", "Both genes")), levels = c("Neither", "One gene", "Both genes"))
d$G12D_vs_other <- factor(ifelse(is.na(d$KRAS) | d$KRAS == 0, NA, ifelse(d$KRAS_group == "G12D", "G12D", "Non-G12D mutant")), levels = c("Non-G12D mutant", "G12D"))
tp <- v[v$Hugo_Symbol == "TP53", ] |> group_by(Tumor_Sample_Barcode) |>
  summarise(TP53_truncation = if (n() != 1) NA_character_ else
    if (first(Variant_Classification) %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins")) "Truncating" else
    if (first(Variant_Classification) %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins")) "Nontruncating" else NA_character_, .groups = "drop")
d <- left_join(d, tp, by = "Tumor_Sample_Barcode")
d$TP53_truncation <- factor(d$TP53_truncation, levels = c("Nontruncating", "Truncating"))
d <- repair_group_status(d, v)
d$Age10 <- d$Age / 10
d$TMB_SD <- d$TMB / sd(d$TMB, na.rm = TRUE)
for (nm in c("Sex", "T_stage", "N_stage", "Neoadjuvant", "M_stage", "LVI", "PNI", "R_status")) d[[nm]] <- factor(d[[nm]])
d$Differentiation <- factor(d$Differentiation, levels = c("WD", "MD", "PD"))
stopifnot(!anyDuplicated(d$patient_id), !any(d$patient_duplicate))
write_tsv(d, file.path(out, "private", "analysis_dataset.tsv"), na = "NA")
write_tsv(vc$audit, file.path(out, "private", "variant_count_audit.tsv"))
tab <- function(x, name) {
  # Detailed result files support figures and reproducibility but are not
  # manuscript Supplementary Tables. The publication package labels them Input.
  name <- sub("^Supplementary_Table", "Input_", name)
  readr::write_tsv(as.data.frame(x), file.path(out, "tables", name), na = "NA")
}
fp <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))
pal <- c("#747474", "#246A8A", "#C17438", "#438569", "#8865A3")
theme_set(theme_classic(base_size = 12, base_family = "Pretendard") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 9, hjust = 0), plot.margin = margin(8, 10, 8, 8)))
savefig <- function(p, name, w = 12, h = 7) {
  ggplot2::ggsave(file.path(out, "figures", paste0(name, ".png")), p,
    width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white", limitsize = FALSE)
  # macOS Quartz embeds the installed font without requiring XQuartz/Cairo.
  pdf_device <- if (Sys.info()[["sysname"]] == "Darwin") function(filename, width, height, ...) grDevices::quartz(type="pdf", file=filename, width=width, height=height, family="Pretendard") else grDevices::cairo_pdf
  ggplot2::ggsave(file.path(out, "figures", paste0(name, ".pdf")), p,
    width = w, height = h, device = pdf_device, bg = "white", limitsize = FALSE)
}
imageplot <- function(path) cowplot::ggdraw() + cowplot::draw_grob(grid::rasterGrob(png::readPNG(path), interpolate = FALSE))

# 1. Cohort accounting and descriptive tables. No source comparison is exported.
flow <- data.frame(item = c("Source records", "Duplicate-ID records excluded", "Analysis patients", "M0", "M1", "OS deaths", "Recurrence evaluable", "Recorded recurrences", "TMB observed", "TMB missing"),
  n = c(nrow(input$clinical_all), sum(input$clinical_all$patient_duplicate), nrow(d), sum(d$M_stage == "0"), sum(d$M_stage == "1"),
    sum(d$survive), sum(complete.cases(d[, c("RFS_m", "Recur")])), sum(d$Recur, na.rm = TRUE), sum(!is.na(d$TMB)), sum(is.na(d$TMB))))
flow <- rbind(flow, data.frame(item = "Duplicate patient IDs excluded",
  n = length(unique(input$clinical_all$patient_id[input$clinical_all$patient_duplicate]))))
tab(flow, "Input_Cohort_flow.tsv")
fup <- survival::survfit(survival::Surv(OS_m, 1 - survive) ~ 1, data = d)
tab(data.frame(method = "Reverse Kaplan-Meier", median_months = unname(summary(fup)$table["median"])), "Input_Followup.tsv")
descvars <- c("Age", "Sex", "M_stage", "Neoadjuvant", "Preop_platinum_exposure", "T_stage", "N_stage", "Differentiation", "LVI", "PNI", "R_status", "Adjuvant_treatment", "TMB")
describe <- function(z, fields) bind_rows(lapply(fields, function(nm) {
  x <- z[[nm]]; n <- sum(!is.na(x)); miss <- sum(is.na(x))
  if (is.numeric(x)) {
    q <- if (n) quantile(x, c(.25,.5,.75), na.rm = TRUE) else rep(NA, 3)
    return(data.frame(variable = nm, category = "Median [IQR]", value = sprintf("%.1f [%.1f, %.1f]", q[2],q[1],q[3]), observed = n, missing = miss))
  }
  tt <- table(factor(as.character(x), levels = sort(unique(as.character(d[[nm]])))))
  data.frame(variable = nm, category = names(tt), value = sprintf("%d (%.1f%%)", as.integer(tt), 100 * as.integer(tt) / nrow(z)), observed = n, missing = miss)
}))
tab(describe(d, descvars), "Table1_Cohort_characteristics.tsv")
ktab <- bind_rows(lapply(levels(d$KRAS_group), function(g) {
  z <- d[!is.na(d$KRAS_group) & d$KRAS_group == g, ]; x <- describe(z, descvars[1:11]); x$group <- g; x$group_n <- nrow(z); x
}))
tab(ktab, "Table2_KRAS_characteristics.tsv")
missing <- data.frame(variable = c(descvars, "OS_m", "RFS_m", "Recur", "KRAS_group"))
missing$missing_n <- vapply(missing$variable, function(nm) sum(is.na(d[[nm]])), integer(1))
missing$total_n <- nrow(d)
tab(missing, "Supplementary_Table1_Missingness.tsv")
tab(data.frame(variable = names(d), definition = figure_label(names(d))), "Supplementary_Table1_Variable_labels.tsv")

# 2. Selected-source KRAS annotation and MAF gene oncoplots with true event counts.
draw_onco <- function(gene_list, filename, title, subset_ids = ids, clinical_profile = "main") {
  vv <- v; vv$Hugo_Symbol <- repair_paper_symbol(vv$Hugo_Symbol)
  # Every gene row, including KRAS, shows actual MAF variant classes. The
  # selected-source subtype annotation is a separate clinical/genomic variable.
  # Do not remove the KRAS gene row or invent its variant class from the subtype.
  cc <- oncoplot_variant_counts(vv, subset_ids, gene_list)
  gene_list <- cc$genes$gene[order(-cc$genes$patients, cc$genes$gene)]
  cc <- oncoplot_variant_counts(vv, subset_ids, gene_list)
  classes <- c(Missense = "#399B72", Frameshift = "#D67F37", Nonsense = "#C84659", Splice = "#8068A3", `In-frame` = "#D4AF37", Other = "#617C8B")
  vv$class <- case_when(vv$Variant_Classification == "Missense_Mutation" ~ "Missense",
    vv$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "Frameshift",
    vv$Variant_Classification == "Nonsense_Mutation" ~ "Nonsense", vv$Variant_Classification == "Splice_Site" ~ "Splice",
    vv$Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins") ~ "In-frame", TRUE ~ "Other")
  mm <- matrix("", length(gene_list), length(subset_ids), dimnames = list(gene_list, subset_ids))
  zz <- vv[vv$Hugo_Symbol %in% gene_list & vv$Tumor_Sample_Barcode %in% subset_ids, ]
  for (i in seq_len(nrow(zz))) {
    a <- zz$Hugo_Symbol[i]; b <- zz$Tumor_Sample_Barcode[i]
    mm[a,b] <- paste(unique(c(strsplit(mm[a,b], ";", fixed = TRUE)[[1]], zz$class[i]))[nzchar(unique(c(strsplit(mm[a,b], ";", fixed = TRUE)[[1]], zz$class[i])))], collapse = ";")
  }
  order_cols <- order(-colSums(mm != ""), -cc$patients$variant_count, subset_ids)
  mm <- mm[, order_cols, drop = FALSE]
  cc <- oncoplot_variant_counts(vv, colnames(mm), rownames(mm))
  clinical_annotations <- oncoplot_clinical_data(d, colnames(mm), clinical_profile)
  ca <- oncoplot_count_annotations(cc, 11)
  legend_params <- lapply(clinical_annotations$legend, function(x) c(x,list(title_gp=grid::gpar(fontsize=12,fontface="bold"),labels_gp=grid::gpar(fontsize=11))))
  ha <- c(ca$top, ComplexHeatmap::HeatmapAnnotation(df = clinical_annotations$data,
    col = clinical_annotations$colors, annotation_label = clinical_annotations$labels,
    annotation_legend_param = legend_params, na_col = "#BDBDBD",
    simple_anno_size=grid::unit(4.2,"mm"), gap=grid::unit(.7,"mm"), annotation_name_gp = grid::gpar(fontsize = 12)))
  af <- function(x,y,w,h,z) {
    grid::grid.rect(x,y,w,h, gp = grid::gpar(fill = "#F0F0F0", col = NA))
    active <- names(z)[z]
    if (length(active)) for (j in seq_along(active)) grid::grid.rect(x, y-h*.45+(j-.5)*h*.9/length(active),w*.95,h*.9/length(active), gp=grid::gpar(fill=classes[active[j]], col=NA))
  }
  ht <- ComplexHeatmap::oncoPrint(mm, alter_fun = af, alter_fun_is_vectorized = FALSE, col = classes,
    top_annotation = ha, right_annotation = ca$right, left_annotation = ca$left, show_pct = FALSE,
    column_order = seq_len(ncol(mm)), row_order = seq_len(nrow(mm)),
    remove_empty_columns = FALSE, remove_empty_rows = FALSE, show_column_names = FALSE,
    column_title = paste0(title, " (N = ", ncol(mm), ")"), column_title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
    row_names_gp = grid::gpar(fontsize = 12), heatmap_legend_param = list(title = "Variant class",
      at=names(classes)[names(classes) %in% unique(unlist(strsplit(mm[nzchar(mm)],";",fixed=TRUE)))],
      labels=names(classes)[names(classes) %in% unique(unlist(strsplit(mm[nzchar(mm)],";",fixed=TRUE)))],
      title_gp=grid::gpar(fontsize=12,fontface="bold"),labels_gp=grid::gpar(fontsize=11)))
  h <- if (length(gene_list) > 8) 9.2 else 7.8
  render <- function() {
    grid::grid.newpage(); grid::pushViewport(grid::viewport(gp=grid::gpar(fontfamily="Pretendard")))
    ComplexHeatmap::draw(ht, merge_legend=TRUE, padding=grid::unit(c(20,4,5,4), "mm"),
      annotation_legend_list=list(ComplexHeatmap::Legend(title="Clinical missingness",labels="Missing",legend_gp=grid::gpar(fill="#BDBDBD"))))
    kras_annotation_source <- switch(config$kras_source,
      clinical = "clinical KRAS_subtype",
      maf = "KRAS subtype re-derived from the MAF",
      workbook_maf = "workbook KRAS_subtype_MAF")
    grid::grid.text(paste0("Gene rows: MAF variant calls. KRAS subtype annotation: ", kras_annotation_source, ".\nTop: all retained unique variants per patient. Right: variants per gene. Left: patients n (% of displayed N).\nGrey mutation cells: no retained call; grey clinical tracks: missing. Gene-level assay coverage remains unverified.\n",
      if(clinical_profile=="main") "CA19-9 / CEA colour saturates at the 80th percentile (legend >=). Reported TMB and tumour-size source units are unverified." else "Reported TMB uses the reported source scale; its unit and assay method remain unverified."),
      x=.01,y=.018,just=c("left","bottom"),gp=grid::gpar(fontsize=9))
  }
  ragg::agg_png(file.path(out,"figures",paste0(filename,".png")), width=17,height=h,units="in",res=300); render(); dev.off()
  if (Sys.info()[["sysname"]] == "Darwin") grDevices::quartz(type="pdf",file=file.path(out,"figures",paste0(filename,".pdf")),width=17,height=h,family="Pretendard") else grDevices::cairo_pdf(file.path(out,"figures",paste0(filename,".pdf")),width=17,height=h,family="Pretendard")
  render(); dev.off()
  tab(cc$genes, paste0("Input_", filename, "_gene_counts.tsv"))
  tab(clinical_annotations$audit, paste0("Input_", filename, "_clinical_tracks.tsv"))
  # Patient identifiers are kept outside publication tables.
  write_tsv(cc$patients, file.path(out,"private",paste0(filename,"_patient_counts.tsv")))
}
top <- v |> distinct(Tumor_Sample_Barcode,Hugo_Symbol) |> count(Hugo_Symbol, sort=TRUE) |> slice_head(n=config$top_n) |> pull(Hugo_Symbol)
top <- unique(c("KRAS",top)) # Always retain the principal PDAC driver in the main overview.
draw_onco(top, "Figure1B_Oncoplot", "PDAC gene variants and clinical characteristics")
kfreq <- d |> count(KRAS_group, name="n") |> mutate(percent=100*n/nrow(d))
tab(kfreq,"Input_KRAS_distribution.tsv")
kdist <- ggplot(kfreq,aes(KRAS_group,n,fill=KRAS_group))+geom_col(width=.7)+geom_text(aes(label=sprintf("%d (%.1f%%)",n,percent)),vjust=-.4,size=3.8)+scale_fill_manual(values=pal)+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(x=NULL,y="Patients",title="KRAS subtype distribution")+theme(legend.position="none",axis.text.x=element_text(size=9))
savefig(kdist,"Figure1C_KRAS_distribution",8,4.5)
flow_spec <- cohort_flow_spec(flow)
jsonlite::write_json(flow_spec,file.path(out,"private","cohort_flow_diagram.json"),auto_unbox=TRUE,pretty=TRUE)
flowp <- plot_cohort_flow(flow_spec)
savefig(flowp,"Figure1A_Cohort",10,5.5)
fig1 <- cowplot::plot_grid(flowp,imageplot(file.path(out,"figures","Figure1B_Oncoplot.png")),kdist,
  ncol=1,rel_heights=c(.28,.52,.20))
savefig(fig1,"Figure1_Cohort_and_genomics",17,21)

# 3. Pathology association family: all 4 genes x 6 outcomes, BH across 24 tests.
cmh_list <- counts_list <- pair_list <- list()
for(g in genes) for(nm in c("T_stage","N_stage","Differentiation","LVI","PNI","R_status")) {
  z <- d[complete.cases(d[,c(g,nm,"Neoadjuvant")]),]
  a <- table(droplevels(factor(z[[nm]])),factor(z[[g]],levels=0:1),droplevels(z$Neoadjuvant))
  fit <- tryCatch(mantelhaen.test(a,correct=FALSE),error=identity)
  key <- paste(g,nm)
  cmh_list[[key]] <- data.frame(gene=g,characteristic=nm,n=nrow(z),p=if(inherits(fit,"error")) NA_real_ else fit$p.value)
  pw <- pairwise_categorical(z[[g]],z[[nm]],z$Neoadjuvant)
  if(nrow(pw)) pair_list[[key]] <- cbind(data.frame(gene=g,characteristic=nm),pw)
  counts_list[[key]] <- z |> transmute(gene=g,characteristic=nm,category=as.character(.data[[nm]]),status=ifelse(.data[[g]]==1,"Variant detected","Not detected")) |> count(gene,characteristic,category,status,name="n") |> group_by(gene,characteristic,status) |> mutate(denominator=sum(n),percent=100*n/denominator) |> ungroup()
}
cmh <- bind_rows(cmh_list); cmh$q_BH <- p.adjust(cmh$p,"BH")
cp <- bind_rows(pair_list); ct <- bind_rows(counts_list)
tab(cmh,"Supplementary_Table2_Pathology_CMH.tsv"); tab(cp,"Supplementary_Table2_Pathology_pairs.tsv"); tab(ct,"Input_Pathology_counts.tsv")
heat <- ggplot(cmh,aes(gene,figure_label(characteristic),fill=-log10(pmax(q_BH,1e-12))))+geom_tile(color="white")+geom_text(aes(label=paste0("p ",fp(p),"\nq ",fp(q_BH))),size=3.5)+scale_fill_gradient(low="#F5F7F8",high="#7BAAC1")+labs(x=NULL,y=NULL,fill="-log10 q",title="Gene-pathology associations",caption="Generalized CMH stratified by neoadjuvant treatment. BH adjustment across 24 tests.")
detail <- function(g,nm) {
  a <- ct[ct$gene==g & ct$characteristic==nm,]; pw <- cp[cp$gene==g & cp$characteristic==nm,]
  ggplot(a,aes(status,percent,fill=category))+geom_col(width=.65)+geom_text(aes(label=paste0(n,"\n",sprintf("%.1f%%",percent))),position=position_stack(vjust=.5),size=3.7)+scale_fill_manual(values=pal[2:5])+labs(x=NULL,y="Within-group patients (%)",fill=figure_label(nm),title=paste(g,figure_label(nm)),caption=pairwise_caption(pw,65))
}
savefig(heat,"Figure2A_Pathology_matrix",10,6)
savefig(detail("TP53","Differentiation"),"Figure2B_TP53_grade",7,6)
savefig(detail("CDKN2A","N_stage"),"Figure2C_CDKN2A_nodes",7,6)
savefig(cowplot::plot_grid(heat,cowplot::plot_grid(detail("TP53","Differentiation"),detail("CDKN2A","N_stage"),nrow=1),ncol=1,rel_heights=c(.48,.52)),"Figure2_Genes_and_pathology",14,11)

# 4. Survival. Full coefficient families stay in supplementary tables even when
# only selected terms are drawn. No p-value-driven covariate selection.
basic <- c("splines::ns(Age, df = 3)","Sex","T_stage","Neoadjuvant","M_stage")
extended <- c(basic,"N_stage","strata(Differentiation)","LVI","PNI","R_status")
cr <- pr <- dr <- contrasts <- list(); fits <- list()
fit_model <- function(exposure,adjust,family,endpoint="OS") {
  time <- if(endpoint=="OS") "OS_m" else "RFS_m"; event <- if(endpoint=="OS") "survive" else "Recur"
  needed <- unique(c(exposure,all.vars(reformulate(adjust)),time,event))
  z <- droplevels(d[complete.cases(d[,needed]),])
  z$time <- z[[time]]; z$event <- z[[event]]
  adjust <- setdiff(adjust,exposure)
  adjust <- adjust[vapply(adjust,function(a) all(vapply(all.vars(reformulate(a)),function(v)length(unique(z[[v]]))>1,logical(1))),logical(1))]
  form <- reformulate(c(exposure,adjust),response="survival::Surv(time,event)")
  id <- paste(endpoint,family,exposure,sep=": "); warnings <- character()
  f <- tryCatch(withCallingHandlers(survival::coxph(form,data=z,x=TRUE,model=TRUE,ties="efron"),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=identity)
  if(inherits(f,"error")) {dr[[id]] <<- data.frame(model=id,exposure,family,endpoint,n=nrow(z),events=sum(z$event),status=conditionMessage(f)); return(NULL)}
  s <- summary(f); ph <- tryCatch(survival::cox.zph(f),error=function(e) NULL)
  status <- if(length(warnings)||any(!is.finite(s$conf.int))||anyNA(coef(f))) "Unstable" else "Estimated"
  cr[[id]] <<- data.frame(model=id,exposure,family,endpoint,term=rownames(s$coefficients),n=nrow(z),events=sum(z$event),HR=s$conf.int[,1],lower95=s$conf.int[,3],upper95=s$conf.int[,4],p=s$coefficients[,5],is_exposure=seq_len(nrow(s$coefficients)) %in% f$assign[[exposure]],status,row.names=NULL)
  if(!is.null(ph)) pr[[id]] <<- data.frame(model=id,exposure,family,endpoint,term=rownames(ph$table),p=ph$table[,"p"])
  dr[[id]] <<- data.frame(model=id,exposure,family,endpoint,n=nrow(z),events=sum(z$event),parameters=length(coef(f)),events_per_parameter=sum(z$event)/length(coef(f)),status,warnings=paste(warnings,collapse="; "),formula=paste(deparse(form),collapse=" "))
  pw <- pairwise_cox(f,z,exposure)
  if(nrow(pw)) contrasts[[id]] <<- cbind(data.frame(model=id,exposure,family,endpoint),pw)
  fits[[id]] <<- f
  invisible(f)
}
exposures <- c(genes,"KRAS_group","driver_group","Any_driver","TP53_truncation","KRAS_TP53","G12D_vs_other")
for(ex in exposures) {fit_model(ex,character(),"Univariable genomic");fit_model(ex,basic,"Adjusted genomic");fit_model(ex,extended,"Extended genomic")}
fit_model("KRAS_group",basic,"Adjusted recurrence","Recorded recurrence")
fit_model("TMB_SD",basic,"Reported TMB exploratory")
cox <- bind_rows(cr) |> group_by(endpoint,family) |> mutate(q_BH=replace(rep(NA_real_,n()),which(is_exposure),p.adjust(p[is_exposure],"BH"))) |> ungroup()
cox <- coefficient_direction(cox,d)
ph <- bind_rows(pr); diag <- bind_rows(dr); cpair <- bind_rows(contrasts)
cpair$comparison_group <- cpair$group1; cpair$reference_group <- cpair$group2
cpair$event_definition <- ifelse(cpair$endpoint=="OS","all-cause death","recorded recurrence")
cpair$ratio_definition <- "Hazard in comparison_group / hazard in reference_group"
cpair$interpretation <- ratio_direction(cpair$HR,cpair$group1,cpair$group2,cpair$event_definition,cpair$p_holm<.05)
tab(cox,"Supplementary_Table3_Cox_coefficients.tsv");tab(diag,"Supplementary_Table3_Model_diagnostics.tsv");tab(ph,"Supplementary_Table3_PH_tests.tsv");tab(cpair,"Supplementary_Table4_Cox_pairs.tsv")
kmstats <- kmrisk <- kmpairs <- kmmed <- list()
kmplot <- function(variable,endpoint="OS",title=NULL) {
  tm <- if(endpoint=="OS") "OS_m" else "RFS_m"; ev <- if(endpoint=="OS") "survive" else "Recur"
  z <- d[complete.cases(d[,c(variable,tm,ev)]),];z$group<-droplevels(factor(z[[variable]]));z$time<-z[[tm]];z$event<-z[[ev]]
  if(variable %in% genes) levels(z$group)<-c("Not detected","Variant detected")
  if(variable=="driver_group") z$group<-factor(z$group,levels=c("0","1","2","3+"))
  f<-survival::survfit(survival::Surv(time,event)~group,data=z,conf.type="log-log")
  lr<-survival::survdiff(survival::Surv(time,event)~group,data=z);p<-pchisq(lr$chisq,length(lr$n)-1,lower.tail=FALSE)
  pw<-pairwise_logrank(z$time,z$event,z$group);id<-paste(endpoint,variable)
  kmstats[[id]]<<-data.frame(variable,endpoint,n=nrow(z),events=sum(z$event),p)
  kmpairs[[id]]<<-cbind(data.frame(variable,endpoint),pw)
  st<-summary(f,censored=TRUE)
  curve<-data.frame(time=st$time,prob=st$surv,lower=st$lower,upper=st$upper,censor=st$n.censor,group=sub("^group=","",as.character(st$strata))) |> bind_rows(data.frame(time=0,prob=1,lower=1,upper=1,censor=0,group=levels(z$group))) |> arrange(group,time)
  med<-as.data.frame(summary(f)$table);kmmed[[id]]<<-data.frame(variable,endpoint,group=sub("^group=","",rownames(med)),n=med$records,events=med$events,median=med$median,lower95=med[["0.95LCL"]],upper95=med[["0.95UCL"]])
  times<-seq(0,floor(max(z$time)/12)*12,12)
  risk<-expand_grid(group=levels(z$group),time=times) |> rowwise() |> mutate(n_risk=sum(z$group==group & z$time>=time)) |> ungroup()
  kmrisk[[id]]<<-mutate(risk,variable=variable,endpoint=endpoint)
  cols<-setNames(pal[seq_len(nlevels(z$group))],levels(z$group))
  plot<-ggplot(curve,aes(time,prob,color=group))+geom_step(linewidth=.8)+geom_point(data=curve[curve$censor>0,],shape=3,size=.7)+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1),labels=scales::percent)+scale_x_continuous(breaks=times)+labs(title=title,subtitle=sprintf("N=%d, events=%d, global log-rank p %s",nrow(z),sum(z$event),fp(p)),x=NULL,y=if(endpoint=="OS") "Overall survival" else "Freedom from recorded recurrence",color=NULL,caption=pairwise_caption(pw,100))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  risk$group<-factor(risk$group,levels=rev(levels(z$group)))
  rp<-ggplot(risk,aes(time,group,label=n_risk,color=group))+geom_text(size=3.1)+scale_color_manual(values=cols)+scale_x_continuous(breaks=times)+labs(x="Follow-up (months)",y=NULL)+theme_classic(base_family="Pretendard",base_size=10)+theme(legend.position="none",axis.line.y=element_blank(),axis.ticks.y=element_blank())
  cowplot::plot_grid(plot,rp,ncol=1,rel_heights=c(.77,.23),align="v",axis="lr")
}
forest <- function(x,title,labels=NULL) {
  forest_result_plot(x,title,labels)
}
core<-cox |> filter(family=="Extended genomic",is_exposure)
fg<-forest(filter(core,exposure %in% genes),"Major genes and OS",paste0(filter(core,exposure %in% genes)$exposure,": detected\nReference: not detected"))
dc<-kmplot("driver_group",title="Number of mutated major genes")
df<-filter(core,exposure=="driver_group")
fc<-forest(df,"Adjusted driver-gene count",paste0(sub("driver_group","",df$term)," mutated genes\nReference: 1 gene"))
savefig(fg,"Figure3A_Gene_OS",11,4)
savefig(dc,"Figure3B_Driver_count_OS",11,8)
savefig(fc,"Figure3C_Driver_count_adjusted",11,3.5)
savefig(cowplot::plot_grid(fg,dc,fc,ncol=1,rel_heights=c(.26,.49,.25)),"Figure3_Gene_level_OS",13,14)
kk<-kmplot("KRAS_group",title="Overall survival by KRAS subtype")
kp<-cpair |> filter(exposure=="KRAS_group",family=="Extended genomic",endpoint=="OS")
kf<-forest_result_plot(kp,"Adjusted KRAS subtype comparisons",paste0(kp$group1,"\nReference: ",kp$group2),p_column="p_holm",p_heading="Holm p")+labs(caption="HR = death hazard in the first group / reference group. HR > 1: higher hazard in the first group.\nHolm p across 10 pairs; pointwise 95% CIs. Age spline, sex, T/N/M, grade, treatment, LVI/PNI and margin adjusted.")
savefig(kk,"Figure4A_KRAS_OS",12,8)
savefig(kf,"Figure4B_KRAS_pairwise",14,7.5)
savefig(cowplot::plot_grid(kk,kf,ncol=1,rel_heights=c(.58,.42)),"Figure4_KRAS_subtype_OS",14,13)

# 5. Supplements. Repair gene lists are descriptive, not functional diagnoses.
sets<-repair_gene_sets()
for(s in names(sets)) draw_onco(sets[[s]],paste0("Supplementary_Figure",if(s=="HRD")1 else 2,"_",s),paste(s,if(s=="HRD")"related gene variants" else "gene variants"),clinical_profile="repair")
tab(bind_rows(lapply(names(sets),function(s)data.frame(gene_list=s,gene=sets[[s]],definition="At least one retained variant; not a functional deficiency diagnosis"))),"Supplementary_Table5_Repair_gene_lists.tsv")
repair_profile<-bind_rows(lapply(c("HRD_variant_status","MMR_variant_status"),function(nm) bind_rows(lapply(unique(d[[nm]]),function(g){a<-describe(d[d[[nm]]==g,],descvars[1:11]);a$gene_set<-nm;a$group<-as.character(g);a}))))
tab(repair_profile,"Supplementary_Table5_Repair_clinical_profiles.tsv")
rk<-kmplot("KRAS_group","Recorded recurrence","KRAS subtype and recorded recurrence")
rp<-d[d$Recur==1 & !is.na(d$Recurrence_pattern),] |> count(Recurrence_pattern,name="n") |> mutate(denominator=sum(n),percent=100*n/denominator)
tab(rp,"Supplementary_Table6_Recurrence_patterns.tsv")
prplot<-ggplot(rp,aes(Recurrence_pattern,percent,fill=Recurrence_pattern))+geom_col()+geom_text(aes(label=sprintf("%d (%.1f%%)",n,percent)),vjust=-.3)+scale_fill_manual(values=pal[2:5])+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(x=NULL,y="Among recorded recurrences (%)",title="Pattern among patients with recorded recurrence")+theme(legend.position="none")
savefig(cowplot::plot_grid(rk,prplot,ncol=1,rel_heights=c(.7,.3)),"Supplementary_Figure3_Recurrence",13,11)
# Reported TMB has unknown assay details and units: retain source values, never
# substitute count, invent a mut/Mb unit, or apply an external high/low cutoff.
tmbsummary<-data.frame(total=nrow(d),observed=sum(!is.na(d$TMB)),missing=sum(is.na(d$TMB)),median=median(d$TMB,na.rm=TRUE),q1=unname(quantile(d$TMB,.25,na.rm=TRUE)),q3=unname(quantile(d$TMB,.75,na.rm=TRUE)),sd=sd(d$TMB,na.rm=TRUE),unit="Unverified reported scale")
tab(tmbsummary,"Supplementary_Table7_Reported_TMB.tsv")
ta<-tpair<-list()
for(nm in c("KRAS_group","Stage_Group","Neoadjuvant","Differentiation","HRD_variant_status","MMR_variant_status")) {
  z<-d[complete.cases(d[,c("TMB",nm)]),]; p<-tryCatch(kruskal.test(z$TMB,factor(z[[nm]]))$p.value,error=function(e)NA_real_)
  ta[[nm]]<-data.frame(variable=nm,n=nrow(z),p)
  pw<-pairwise_numeric(z$TMB,z[[nm]]); if(nrow(pw)) tpair[[nm]]<-cbind(data.frame(variable=nm),pw)
}
ta<-bind_rows(ta);ta$q_BH<-p.adjust(ta$p,"BH");tpair<-bind_rows(tpair)
tab(ta,"Supplementary_Table7_TMB_associations.tsv");tab(tpair,"Supplementary_Table7_TMB_pairs.tsv")
histp<-ggplot(d[!is.na(d$TMB),],aes(TMB))+geom_histogram(binwidth=2,fill=pal[2],color="white")+labs(x="Reported TMB (source scale)",y="Patients",title="Reported TMB distribution",subtitle=sprintf("Observed %d, missing %d",tmbsummary$observed,tmbsummary$missing))
boxp<-ggplot(d[!is.na(d$TMB),],aes(KRAS_group,TMB,fill=KRAS_group))+geom_boxplot(outlier.alpha=.25)+scale_fill_manual(values=pal)+labs(x=NULL,y="Reported TMB (source scale)",title="Reported TMB by KRAS subtype",caption=pairwise_caption(tpair[tpair$variable=="KRAS_group",],85))+theme(legend.position="none",axis.text.x=element_text(size=9))
savefig(cowplot::plot_grid(histp,boxp,nrow=1),"Supplementary_Figure4_Reported_TMB",14,6)
sx<-cox |> filter(family %in% c("Adjusted genomic","Extended genomic"),exposure %in% genes,is_exposure)
savefig(forest(sx,"Core and extended adjustment",paste(sx$exposure,ifelse(sx$family=="Adjusted genomic","Core","Extended"))) +
  labs(caption="All comparisons: variant detected / not detected. HR > 1: higher death hazard in the detected group.\nPointwise 95% CIs; BH q values within each adjustment family. Extended models stratify by differentiation."),
  "Supplementary_Figure5_Model_sensitivity",14,7)
tab(bind_rows(kmstats),"Supplementary_Table4_Logrank.tsv");tab(bind_rows(kmpairs),"Supplementary_Table4_Logrank_pairs.tsv");tab(bind_rows(kmmed),"Supplementary_Table4_KM_estimates.tsv");tab(bind_rows(kmrisk),"Supplementary_Table4_Number_at_risk.tsv")
if(isTRUE(config$clinical_extensions)) {
  source("R/clinical_extensions.R")
  run_clinical_extensions(d,config,out,savefig,tab,pal,extended)
}
story_results <- run_manuscript_story_extensions(
  d = d, out = out, savefig = savefig, tab = tab, pal = pal, extended = extended,
  heat = heat, fg = fg, kk = kk, kf = kf, core = core, cpair = cpair
)
extension_results <- run_recurrence_signature_extensions(
  d = d, variants = v, out = out, savefig = savefig, tab = tab, pal = pal,
  config = config
)
provenance<-list(kras_source=config$kras_source,kras_source_column=if(config$kras_source=="workbook_maf") "KRAS_subtype_MAF" else if(config$kras_source=="clinical") "KRAS_subtype" else "MAF-derived",kras_source_final=TRUE,pretreatment_measurements=config$pretreatment_measurements,
  timing_basis="Investigator statement: diagnosis before neoadjuvant treatment; surgical pathology remains postoperative information",
  survival_origin="Surgery, per source dictionary",clinical_extensions=config$clinical_extensions,validation_bootstraps=config$validation_bootstraps,
  cohort_n=nrow(d),excluded_rows=sum(input$clinical_all$patient_duplicate),m0=sum(d$M_stage=="0"),m1=sum(d$M_stage=="1"),os_events=sum(d$survive),tmb_observed=sum(!is.na(d$TMB)),tmb_missing=sum(is.na(d$TMB)),tmb_unit="Unverified reported scale",main_figures=5,main_tables=2,source_comparisons_included=FALSE,signature_analysis=config$signature_analysis,signature_status=extension_results$signature_status,signature_interpretation="Exploratory targeted-panel SBS analysis; no treatment-induced claim",variant_count_definition="Unique sample / genome build / chromosome / position / ref / alt events after the existing MAF retention filter",cox_adjustment=extended,cox_bh_family="All exposure coefficients across the 10 genomic models within each adjustment family (16 contrasts in this cohort)",input_md5=as.list(tools::md5sum(c(config$clinical_file,config$maf_file))))
jsonlite::write_json(provenance,file.path(out,"private","provenance.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(c("# Publication outputs","",
  paste0("Clinical input: ", basename(config$clinical_file), " (sheet ", config$clinical_sheet, ")."),
  paste0("Genomic input: ", basename(config$maf_file), "."),
  "Main figures: 1-5. Main tables: 1-2. Selected supporting analyses are labelled Supplementary Figures.",
  "The original three-group recurrence analysis is retained; distant-only and five-level expanded recurrence analyses are added.",
  "Detailed calculation files in tables/ are labelled Input_*.tsv; no Supplementary Table series is generated.",
  "Patient-level data and source hashes stay in private/; never distribute this folder.",
  "Reported TMB assay details and units remain unverified. No count substitution or external threshold.",
  "Signature analysis is exploratory because the MAF contains selected coding variants, per-sample SNV counts are low, and no panel BED/callable territory is available.",
  "All study samples are considered preoperative-treatment-naive; SBS31/35 results are associations, not treatment-induced effects.",
  "No clinical/MAF source comparison is included in this publication package."),file.path(out,"README.md"))
message("Manuscript analyses complete: ",normalizePath(out))
