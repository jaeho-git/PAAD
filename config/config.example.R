# Copy this file to config/local.R and edit only the local copy.
# Actual patient data and generated outputs must not be committed.

paad_config <- list(
  clinical_file = "data/raw/PDAC_clinical_info.xlsx",
  target_patients_file = "data/raw/result_527명.txt",
  maf_file = "data/raw/PDAC_oncopanel.maf",
  output_dir = file.path(PROJECT_ROOT, "outputs", "local_run"),

  target_encoding = "CP949",
  clinical_sheet = 1,

  # Preserve the analysis definitions and filenames of the v19 integrated script.
  compatibility_mode = "v19",

  # v19 requests `Sex`, but the observed clinical workbook uses `성별코드` and
  # the target-patient file uses `sex`. The default keeps v19 behavior and does
  # not silently substitute either field. Explicit alternatives are:
  #   "target_patients"  - map target-patient `sex` to `Sex`
  #   "clinical_korean"  - map clinical-workbook `성별코드` to `Sex`
  sex_column_policy = "legacy",

  # v19 uses simulated Fisher tests without setting a seed. Keep NULL for exact
  # compatibility, or set an integer locally for reproducible simulated p-values.
  random_seed = NULL,

  top_n = 20L,
  # 1000 preserves v19 production figures. A lower value can be used locally
  # for a faster smoke test without changing the analysis data or filenames.
  plot_dpi = 1000L
)
