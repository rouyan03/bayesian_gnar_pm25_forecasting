root <- normalizePath(getwd(), mustWork = TRUE)
table_dir <- file.path(root, "outputs", "tables")
figure_dir <- file.path(root, "outputs", "figures")

near <- function(x, value, tolerance = 1e-6) {
  is.finite(x) && abs(x - value) < tolerance
}

check_file <- function(path) {
  data.frame(
    check = path,
    status = if (file.exists(file.path(root, path))) "ok" else "missing",
    stringsAsFactors = FALSE
  )
}

checks <- do.call(rbind, lapply(c(
  "outputs/tables/baseline_forecasts_2024.csv",
  "outputs/tables/station_ar_test_metrics_2024.csv",
  "outputs/tables/unrestricted_var_2024.csv",
  "outputs/tables/var_benchmark_2024.csv",
  "outputs/tables/ols_gnar_selected_test_2024_corrected.csv",
  "outputs/tables/bayes_gaussian_gnar_multichain_2024.csv",
  "outputs/tables/hierarchical_station_intercept_bayes_gnar_2024.csv",
  "outputs/tables/hierarchical_missingness_summary_with_ar1_2024.csv",
  "outputs/tables/hierarchical_peak_diagnostics_by_decile_2024.csv",
  "outputs/tables/bayes_student_t_gnar_nu_validation_3_to_50_2024.csv",
  "outputs/tables/peak_student_t_weighted_comparison_2024.csv",
  "outputs/tables/external_year_pm25_model_comparison_2025.csv",
  "outputs/tables/bayesian_gnar_synthetic_summary.csv",
  "outputs/figures/corbit_pnacf_equal_neighbour_2024_selection.png",
  "outputs/figures/ols_gnar_selection_validation_heatmap_2024_thesis_blue_red.png",
  "outputs/figures/poster_hierarchical_bayes_gnar_predictive_intervals_6stations_2024.png",
  "outputs/figures/poster_missingness_robustness_compact_with_ar1_2024.png",
  "outputs/figures/poster_hierarchical_peak_calibration_bias_2024.png",
  "outputs/figures/peak_models_decile_forecast_error_comparison_2024.png"
), check_file))

add_check <- function(name, ok) {
  checks <<- rbind(
    checks,
    data.frame(check = name, status = if (ok) "ok" else "failed")
  )
}

baseline <- read.csv(file.path(table_dir, "baseline_forecasts_2024.csv"))
add_check("2024 naive RMSPE", near(baseline$rmse[1], 4.303624))

ar1 <- read.csv(file.path(table_dir, "station_ar_test_metrics_2024.csv"))
add_check("2024 AR(1) RMSPE", near(ar1$rmse[1], 3.970447))

ols_var <- read.csv(file.path(table_dir, "unrestricted_var_2024.csv"))
add_check("2024 unrestricted VAR RMSPE", near(ols_var$test_rmse[1], 4.158548))

ridge_var <- read.csv(file.path(table_dir, "var_benchmark_2024.csv"))
add_check("2024 ridge VAR RMSPE", near(ridge_var$test_rmse[1], 3.821319))
add_check("2024 ridge VAR lambda", ridge_var$selected_lambda[1] == 1000)

ols_gnar <- read.csv(file.path(table_dir, "ols_gnar_selected_test_2024_corrected.csv"))
add_check("2024 selected GNAR order", ols_gnar$model_order[1] == "GNAR(2,[1,1])")
add_check("2024 selected GNAR network", ols_gnar$network[1] == "equal_neighbour")
add_check("2024 OLS GNAR RMSPE", near(ols_gnar$test_rmse[1], 3.833179))

bayes <- read.csv(file.path(table_dir, "bayes_gaussian_gnar_multichain_2024.csv"))
add_check("2024 Bayesian GNAR coverage", near(bayes$coverage_95[1], 0.9279279))

hier <- read.csv(file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_2024.csv"))
add_check("2024 hierarchical GNAR RMSPE", near(hier$rmse[1], 3.787855))
add_check("2024 hierarchical GNAR coverage", near(hier$coverage_95[1], 0.927027))

missingness <- read.csv(file.path(table_dir, "hierarchical_missingness_summary_with_ar1_2024.csv"))
add_check("missingness table rows", nrow(missingness) == 28)
add_check(
  "hierarchical GNAR 40 percent missingness RMSPE",
  near(missingness$rmse[
    missingness$model == "Hierarchical Bayesian GNAR" &
      missingness$missing_percent == 40
  ], 3.814103)
)

nu <- read.csv(file.path(table_dir, "bayes_student_t_gnar_nu_validation_3_to_50_2024.csv"))
add_check("Student-t nu grid", nrow(nu) == 48 && identical(sort(nu$nu), 3:50))
add_check("Student-t validation selected nu", nu$nu[which.min(nu$validation_rmspe)] == 46)

peak <- read.csv(file.path(table_dir, "peak_student_t_weighted_comparison_2024.csv"))
add_check("peak table rows", nrow(peak) == 3)
add_check(
  "peak-weighted peak RMSPE",
  near(peak$peak_rmspe[peak$model == "Peak-weighted GNAR"], 7.523920)
)
add_check(
  "Student-t fitted nu",
  peak$selected_nu[peak$model == "Student-t Bayesian GNAR"] == 3
)

external <- read.csv(file.path(table_dir, "external_year_pm25_model_comparison_2025.csv"))
add_check("external-year table rows", nrow(external) == 7)
add_check(
  "2025 hierarchical GNAR RMSPE",
  near(external$rmse[external$model == "Hierarchical Bayesian GNAR"], 4.498542)
)
add_check("2025 observed station-days", unique(external$observed_test_cells) == 10778)

simulation <- read.csv(file.path(table_dir, "bayesian_gnar_synthetic_summary.csv"))
add_check("simulation parameter coverage", simulation$parameter_coverage_95[1] == 1)
add_check("simulation forecast coverage", near(simulation$coverage_95[1], 0.9511111))

print(checks, row.names = FALSE)

if (any(checks$status != "ok")) {
  stop("Some checks failed.")
}

message("All checked thesis outputs are present and match the expected values.")
