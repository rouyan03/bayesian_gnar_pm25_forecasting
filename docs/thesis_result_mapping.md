# Thesis result mapping

This file links the main thesis results to the cleaned code and output files.

| Thesis item | Script | Output |
|---|---|---|
| 2024 data summary | `04_summarise_data_and_networks.R` | `outputs/tables/main_2024_site_summary.csv` |
| network summary | `04_summarise_data_and_networks.R` | `outputs/tables/network_summary_2024_nodes.csv` |
| station map | `05_plot_data_diagnostics.R` | `outputs/figures/poster_station_map_mean_only_2024.png` |
| example PM2.5 time series | `05_plot_data_diagnostics.R` | `outputs/figures/pm25_daily_2024_example_sites_6.png` |
| corbit diagnostics | `06_corbit_diagnostics.R` | `outputs/figures/corbit_pnacf_equal_neighbour_2024_selection.png` |
| GNAR network/order selection | `07_select_gnar_network_and_order.R` | `outputs/tables/ols_gnar_selection_validation_2024_corrected.csv` |
| selection heatmap | `07_select_gnar_network_and_order.R` | `outputs/figures/ols_gnar_selection_validation_heatmap_2024_thesis_blue_red.png` |
| naive baseline | `08_fit_naive_baseline.R` | `outputs/tables/baseline_forecasts_2024.csv` |
| AR(1) benchmark | `09_fit_ar1_benchmark.R` | `outputs/tables/station_ar_test_metrics_2024.csv` |
| unrestricted VAR | `10_fit_unrestricted_var.R` | `outputs/tables/unrestricted_var_2024.csv` |
| ridge VAR | `11_fit_ridge_var.R` | `outputs/tables/ols_var_vs_ridge_var_2024.csv` |
| OLS GNAR | `07_select_gnar_network_and_order.R` | `outputs/tables/ols_gnar_selected_test_2024_corrected.csv` |
| Bayesian Gaussian GNAR | `12_fit_bayesian_gaussian_gnar.R` | `outputs/tables/bayes_gaussian_gnar_multichain_2024.csv` |
| Bayesian Gaussian MCMC diagnostics | `12_fit_bayesian_gaussian_gnar.R` | `outputs/tables/bayes_gaussian_gnar_multichain_diagnostics_2024.csv` |
| hierarchical Bayesian GNAR | `13_fit_hierarchical_bayesian_gnar.R` | `outputs/tables/hierarchical_station_intercept_bayes_gnar_coefficients_2024.csv` |
| hierarchical predictive intervals | `14_plot_hierarchical_intervals.R` | `outputs/figures/poster_hierarchical_bayes_gnar_predictive_intervals_6stations_2024.png` |
| simulation validation | `15_simulation_validation.R` | `outputs/tables/bayesian_gnar_synthetic_summary.csv` |
| missing-data robustness | `16_missing_data_robustness.R` | `outputs/tables/hierarchical_missingness_summary_with_ar1_2024.csv` |
| missing-data robustness figure | `16_missing_data_robustness.R` | `outputs/figures/poster_missingness_robustness_compact_with_ar1_2024.png` |
| calibration and forecast error by decile | `17_calibration_by_decile.R` | `outputs/tables/hierarchical_peak_diagnostics_by_decile_2024.csv` |
| calibration by decile figure | `17_calibration_by_decile.R` | `outputs/figures/poster_hierarchical_peak_calibration_bias_2024.png` |
| Student-t sensitivity | `18_student_t_sensitivity.R` | `outputs/tables/bayes_student_t_gnar_nu_validation_3_to_50_2024.csv` |
| peak-weighted GNAR | `19_peak_weighted_gnar.R` | `outputs/tables/peak_student_t_weighted_comparison_2024.csv` |
| peak forecast-error figure | `19_peak_weighted_gnar.R` | `outputs/figures/peak_models_decile_forecast_error_comparison_2024.png` |
| external-year validation | `20_external_year_validation_2025.R` | `outputs/tables/external_year_pm25_model_comparison_2025.csv` |
