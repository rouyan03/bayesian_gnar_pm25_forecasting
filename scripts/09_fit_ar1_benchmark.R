project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(project_root, "data", "processed", "analysis")
table_dir <- file.path(project_root, "outputs", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "R", "gnar_bayes_utils.R"))

panel_file <- file.path(analysis_dir, "main_2024_pm25_daily_panel.csv")
if (!file.exists(panel_file)) {
  stop("Missing main 2024 panel. Run final_code/scripts/02_prepare_2024_panel.R.")
}

panel <- read.csv(panel_file, check.names = FALSE, stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
site_ids <- setdiff(names(panel), "date")

X <- as.matrix(panel[, site_ids, drop = FALSE])
storage.mode(X) <- "numeric"
if (any(!is.finite(X))) stop("The AR(1) benchmark requires a complete 2024 panel.")
#fixed 2024 train-val-test split
split <- analysis_split_2024(panel$date)
#fit AR(1) model to one station
fit_ar1 <- function(series, end_index) {
  response <- series[2:end_index]
  lag1 <- series[1:(end_index - 1)]
  design <- cbind(intercept = 1, lag1 = lag1)
  coefficients <- lm.fit(design, response)$coefficients
  coefficients[is.na(coefficients)] <- 0
  coefficients
}
#store test forecasts and fitted coefficients
forecast <- matrix(
  NA_real_,
  nrow = length(split$test_origins),
  ncol = length(site_ids),
  dimnames = list(NULL, site_ids)
)
coefficients <- matrix(
  NA_real_,
  nrow = length(site_ids),
  ncol = 2,
  dimnames = list(site_ids, c("intercept", "lag1"))
)

for (site in site_ids) {
  coefficients[site, ] <- fit_ar1(X[, site], split$fit_end)
  forecast[, site] <- coefficients[site, "intercept"] +
    coefficients[site, "lag1"] * X[split$test_origins, site]
}

observed <- X[split$test_origins + 1, , drop = FALSE]
errors <- observed - forecast
metrics <- forecast_metrics(observed, forecast)
#save test performance
summary_table <- data.frame(
  model = "AR(1)",
  parameters = 2 * length(site_ids),
  test_start = as.character(panel$date[min(split$test_origins) + 1]),
  test_end = as.character(panel$date[max(split$test_origins) + 1]),
  rmse = metrics$rmse,
  mae = metrics$mae,
  stringsAsFactors = FALSE
)
#save fitted coefficients
coefficient_table <- data.frame(
  site_id = site_ids,
  intercept = coefficients[, "intercept"],
  lag1 = coefficients[, "lag1"],
  stringsAsFactors = FALSE
)

forecast_table <- data.frame(
  date = rep(panel$date[split$test_origins + 1], each = length(site_ids)),
  site_id = rep(site_ids, times = length(split$test_origins)),
  observed = as.vector(t(observed)),
  forecast = as.vector(t(forecast)),
  error = as.vector(t(errors))
)

write.csv(summary_table, file.path(table_dir, "station_ar_test_metrics_2024.csv"), row.names = FALSE)
write.csv(coefficient_table, file.path(table_dir, "station_ar_coefficients_2024.csv"), row.names = FALSE)
write.csv(forecast_table, file.path(table_dir, "station_ar_forecasts_2024.csv"), row.names = FALSE)
