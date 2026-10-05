# Copy to config/local.R; local.R is ignored by Git.
# The new workbook contains Main and Data_dictionary sheets.
paad_config <- list(
  clinical_file = "data/raw/260927_v3_PDAC_ANALYSIS_with_MAF_annotations.xlsx",
  previous_clinical_file = "data/raw/PDAC_clinical_info.xlsx",
  target_patients_file = "data/raw/result_527명.txt",
  maf_file = "data/raw/PDAC_oncopanel.maf",
  output_dir = file.path(PROJECT_ROOT, "outputs", "updated_v3_run"),
  clinical_schema = "curated_v3",
  analysis_profile = "manuscript", # Publication figures only; no source-comparison figures.
  kras_source = "workbook_maf", # Uses Main!KRAS_subtype_MAF; do not mix sources.
  manuscript_subdir = "manuscript",
  # Mutational-signature analysis is optional because it requires a separate,
  # pre-installed SigProfiler environment and a local GRCh37 reference.
  signature_analysis = FALSE,
  signature_python = NULL,
  signature_reference_volume = NULL,
  signature_min_snv = 10L,
  signature_bootstraps = 100L,
  signature_positive_min_count = 5L,
  signature_positive_min_proportion = 0.20,
  signature_min_cosine = 0.90,
  signature_min_stability = 0.80,
  clinical_extensions = FALSE,
  validation_bootstraps = 500L,
  pretreatment_measurements = TRUE, # Investigator-confirmed clinical/molecular timing.
  clinical_sheet = "Main",
  updated_cohort = "all_unique",
  target_encoding = "CP949",
  compatibility_mode = "v19",
  sex_column_policy = "legacy", # Ignored for curated_v3: use the workbook's F/M.
  random_seed = 20260928L,
  top_n = 20L,
  plot_dpi = 1000L
)
