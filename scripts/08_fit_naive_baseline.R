project_root <- normalizePath(getwd(), mustWork = TRUE)

analysis_dir <- file.path(project_root, "data", "processed", "analysis")
table_dir <- file.path(project_root, "outputs", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "R", "gnar_bayes_utils.R"))

panel <- read.csv(
  file.path(analysis_dir, "main_2024_pm25_daily_panel.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
panel$date <- as.Date(panel$date)

site_ids <- setdiff(names(panel), "date")
X <- as.matrix(panel[, site_ids, drop = FALSE])
storage.mode(X) <- "numeric"

split <- analysis_split_2024(panel$date)
test_dates <- panel$date[split$test_origins + 1]
observed <- X[split$test_origins + 1, , drop = FALSE]

#naive forecast (tomorrow = today's observed PM2.5)
forecast <- X[split$test_origins, , drop = FALSE]
metrics <- forecast_metrics(observed, forecast)

summary_table <- data.frame(
  model = "Naive previous-day forecast",
  parameters = 0,
  test_start = as.character(min(test_dates)),
  test_end = as.character(max(test_dates)),
  test_days = length(test_dates),
  test_cells = length(observed),
  rmse = metrics$rmse,
  mae = metrics$mae,
  stringsAsFactors = FALSE
)

forecast_table <- data.frame(
  date = rep(test_dates, each = length(site_ids)),
  site_id = rep(site_ids, times = length(test_dates)),
  observed = as.vector(t(observed)),
  forecast = as.vector(t(forecast))
)

write.csv(summary_table, file.path(table_dir, "baseline_forecasts_2024.csv"), row.names = FALSE)
write.csv(forecast_table, file.path(table_dir, "baseline_forecast_values_2024.csv"), row.names = FALSE)
