scripts <- c(
  "scripts/02_prepare_2024_panel.R",
  "scripts/03_build_station_networks.R",
  "scripts/04_summarise_data_and_networks.R",
  "scripts/05_plot_data_diagnostics.R",
  "scripts/06_corbit_diagnostics.R",
  "scripts/07_select_gnar_network_and_order.R",
  "scripts/08_fit_naive_baseline.R",
  "scripts/09_fit_ar1_benchmark.R",
  "scripts/10_fit_unrestricted_var.R",
  "scripts/11_fit_ridge_var.R",
  "scripts/12_fit_bayesian_gaussian_gnar.R",
  "scripts/13_fit_hierarchical_bayesian_gnar.R",
  "scripts/14_plot_hierarchical_intervals.R",
  "scripts/15_simulation_validation.R",
  "scripts/16_missing_data_robustness.R",
  "scripts/17_calibration_by_decile.R",
  "scripts/18_student_t_sensitivity.R",
  "scripts/19_peak_weighted_gnar.R",
  "scripts/20_external_year_validation_2025.R"
)

required_processed <- c(
  "data/processed/defra_aurn/aurn_pm25_daily_panel_2022.csv",
  "data/processed/defra_aurn/aurn_pm25_daily_panel_2023.csv",
  "data/processed/defra_aurn/aurn_pm25_daily_panel_2024.csv",
  "data/processed/defra_aurn/aurn_pm25_daily_panel_2025.csv",
  "data/processed/defra_aurn/aurn_pm25_selected_sites_fixed_2024_nodes.csv"
)

if (any(!file.exists(required_processed))) {
  message("Processed AURN files not found. Running scripts/01_extract_aurn_data.R")
  source("scripts/01_extract_aurn_data.R")
}

for (script in scripts) {
  message("\nRunning ", script)
  source(script)
}
