# Synthetic-data smoke-test configuration. The files under data/example contain
# no real patient records. Generated artifacts stay under outputs/.

paad_config <- list(
  clinical_file = "data/example/PDAC_clinical_info.example.xlsx",
  target_patients_file = "data/example/result_527명.example.txt",
  maf_file = "data/example/PDAC_oncopanel.example.maf",
  output_dir = file.path(PROJECT_ROOT, "outputs", "synthetic_run"),
  target_encoding = "CP949",
  clinical_sheet = 1,
  compatibility_mode = "v19",
  sex_column_policy = "target_patients",
  random_seed = 20250921L,
  top_n = 20L,
  plot_dpi = 150L
)
