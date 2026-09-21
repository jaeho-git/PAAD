ensure_analysis_data <- function(config, data = NULL) {
  if (is.null(data)) load_analysis_data(config) else data
}

workflow_oncoplots <- function(config, data = NULL) {
  data <- ensure_analysis_data(config, data)
  annotations <- annotation_variables(data)
  message("[1/2] Main PDAC oncoplot")
  draw_custom_oncoplot(data$maf, config$top_n, "Top 20 Mutated Genes in PDAC", output_path(config, "1_Oncoplot_Main_PDAC.png"), annotations, config$plot_dpi)
  valid <- data$clinical |> dplyr::filter(!is.na(T), !is.na(N)) |> dplyr::pull(Tumor_Sample_Barcode)
  if (length(valid)) {
    message("[2/2] T/N-filtered PDAC oncoplot")
    subset <- maftools::subsetMaf(data$maf, tsb = valid)
    draw_custom_oncoplot(subset, config$top_n, "PDAC Oncoplot (Filtered by T & N)", output_path(config, "1_2_Oncoplot_TN_Filtered_PDAC.png"), annotations, config$plot_dpi)
  } else message("Skipping T/N-filtered oncoplot: no samples have both T and N.")
  invisible(TRUE)
}

workflow_clinical_summary <- function(config, data = NULL) {
  data <- ensure_analysis_data(config, data)
  write_clinical_summary(data$clinical, annotation_variables(data), output_path(config, "2_Clinical_Table_PDAC.xlsx"))
}

workflow_comparison_plots <- function(config, data = NULL) {
  data <- ensure_analysis_data(config, data)
  # Reset at the stochastic workflow boundary so an individual run and run_all
  # produce the same simulated Fisher results when a seed is explicitly set.
  reset_workflow_seed(config)
  variables <- intersect(c("Differentiation", "NAC", "T", "N", "N_status", "Stage", "Stage_Group"), names(data$clinical))
  run_comparison_analysis(data$maf, data$clinical, config$output_dir, "PDAC", variables, config$top_n, config$plot_dpi)
}

workflow_heatmaps <- function(config, data = NULL) {
  data <- ensure_analysis_data(config, data)
  annotations <- annotation_variables(data)
  draw_clustered_heatmap(data$maf, config$top_n, "PDAC Clustered Heatmap (Top 20)", output_path(config, "5_Heatmap_Top20_PDAC.png"), annotations, config$plot_dpi)
  draw_clustered_heatmap(data$maf, min(5L, config$top_n), "PDAC Clustered Heatmap (Top 5)", output_path(config, "6_Heatmap_Top5_PDAC.png"), annotations, config$plot_dpi)
  draw_tmb_heatmap(data$maf, "TMB Heatmap (PDAC)", output_path(config, "7_TMB_Heatmap_PDAC.png"), annotations, config$plot_dpi)
  invisible(TRUE)
}

workflow_somatic_interactions <- function(config, data = NULL) {
  data <- ensure_analysis_data(config, data)
  draw_somatic_interactions(data$maf, output_path(config, "9_Somatic_Interactions_PDAC.png"), output_path(config, "9_Somatic_Interactions_Results_PDAC.xlsx"), config$plot_dpi)
}

workflow_stage_counts <- function(config, data = NULL) {
  data <- ensure_analysis_data(config, data)
  draw_stage_counts(data$clinical, output_path(config, "10_Stage_Counts_PDAC.png"), config$plot_dpi)
}

workflow_driver_kras_survival <- function(config, data = NULL) {
  data <- ensure_analysis_data(config, data)
  reset_workflow_seed(config)
  run_driver_kras_outputs(data$maf, data$clinical, config$output_dir, config$plot_dpi)
}

workflow_all <- function(config) {
  data <- load_analysis_data(config)
  workflows <- list(
    oncoplots = workflow_oncoplots,
    clinical_summary = workflow_clinical_summary,
    comparison_plots = workflow_comparison_plots,
    heatmaps = workflow_heatmaps,
    somatic_interactions = workflow_somatic_interactions,
    stage_counts = workflow_stage_counts,
    driver_kras_survival = workflow_driver_kras_survival
  )
  for (name in names(workflows)) {
    message("\n=== Running ", name, " ===")
    tryCatch(workflows[[name]](config, data), error = function(error) {
      stop("Workflow '", name, "' failed: ", conditionMessage(error), call. = FALSE)
    })
  }
  message("\nAll configured PDAC workflows completed successfully.")
  invisible(TRUE)
}
