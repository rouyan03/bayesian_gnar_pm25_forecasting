project_root <- normalizePath(getwd(), mustWork = TRUE)
ver_root <- project_root
analysis_dir <- file.path(ver_root, "data", "processed", "analysis")
table_dir <- file.path(ver_root, "outputs", "tables")
figure_dir <- file.path(ver_root, "outputs", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(project_root, "R", "gnar_bayes_utils.R"))

panel <- read.csv(
  file.path(analysis_dir, "main_2024_pm25_daily_panel.csv"),
  check.names = FALSE
)
panel$date <- as.Date(panel$date)
site_ids <- setdiff(names(panel), "date")
X <- as.matrix(panel[, site_ids, drop = FALSE])
storage.mode(X) <- "numeric"
split <- analysis_split_2024(panel$date)
#build VAR response and lagged predictor matrices
var_design <- function(X, p) {
  Y <- X[(p + 1):nrow(X), , drop = FALSE]
  Z <- matrix(1, nrow(Y), 1)
  colnames(Z) <- "intercept"
  for (lag in seq_len(p)) {
    lagged <- X[(p + 1 - lag):(nrow(X) - lag), , drop = FALSE]
    colnames(lagged) <- paste0(colnames(X), "_lag", lag)
    Z <- cbind(Z, lagged)
  }
  list(Y = Y, Z = Z)
}
#fit ridge-regularised VAR
fit_ridge_var <- function(X, p, lambda) {
  design <- var_design(X, p)
  penalty <- diag(lambda, ncol(design$Z))
  penalty[1, 1] <- 0
  coefficients <- solve(
    crossprod(design$Z) + penalty,
    crossprod(design$Z, design$Y)
  )
  list(coefficients = coefficients, p = p, lambda = lambda)
}

forecast_var <- function(X, fit, origins) {
  forecast <- matrix(NA_real_, length(origins), ncol(X))
  for (index in seq_along(origins)) {
    origin <- origins[index]
    predictor <- 1
    for (lag in seq_len(fit$p)) {
      predictor <- c(predictor, X[origin - lag + 1, ])
    }
    forecast[index, ] <- as.numeric(predictor %*% fit$coefficients)
  }
  colnames(forecast) <- colnames(X)
  forecast
}
#compare VAR orders and ridge penalties on the validation period
grid <- expand.grid(
  p = c(1, 2),
  lambda = c(0.1, 1, 10, 100, 1000),
  KEEP.OUT.ATTRS = FALSE
)
validation <- do.call(rbind, lapply(seq_len(nrow(grid)), function(index) {
  spec <- grid[index, ]
  fit <- fit_ridge_var(
    X[seq_len(split$selection_train_end), , drop = FALSE],
    spec$p,
    spec$lambda
  )
  forecast <- forecast_var(X, fit, split$validation_origins)
  observed <- X[split$validation_origins + 1, , drop = FALSE]
  data.frame(
    model = "Ridge VAR",
    p = spec$p,
    lambda = spec$lambda,
    parameters = ncol(X) * (1 + ncol(X) * spec$p),
    validation_rmse = rmse(observed - forecast),
    validation_mae = mae(observed - forecast)
  )
}))
validation <- validation[order(validation$validation_rmse), ]
write.csv(
  validation,
  file.path(table_dir, "var_benchmark_validation_2024.csv"),
  row.names = FALSE
)
selected <- validation[1, ] #select the best specification and refit on all train+val data
fit <- fit_ridge_var(
  X[seq_len(split$fit_end), , drop = FALSE],
  selected$p,
  selected$lambda
)
forecast <- forecast_var(X, fit, split$test_origins) #eval the selected ridge VAR on test dataset
observed <- X[split$test_origins + 1, , drop = FALSE]
metrics <- forecast_metrics(observed, forecast)
summary <- data.frame( #save final test results
  model = "Ridge VAR",
  model_order = paste0("VAR(", selected$p, ")"),
  selected_lambda = selected$lambda,
  selection_rule = "minimum validation RMSE",
  parameters = selected$parameters,
  test_rmse = metrics$rmse,
  test_mae = metrics$mae,
  note = "VAR uses no network; ridge regularisation shrinks the high-dimensional VAR coefficients"
)
write.csv(summary, file.path(table_dir, "var_benchmark_2024.csv"), row.names = FALSE)
write.csv(
  data.frame(
    date = rep(panel$date[split$test_origins + 1], each = ncol(X)),
    site_id = rep(site_ids, times = length(split$test_origins)),
    observed = as.vector(t(observed)),
    forecast = as.vector(t(forecast))
  ),
  file.path(table_dir, "var_benchmark_forecasts_2024.csv"),
  row.names = FALSE
)
