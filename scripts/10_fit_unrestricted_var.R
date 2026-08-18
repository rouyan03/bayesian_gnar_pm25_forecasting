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
  check.names = FALSE,
  stringsAsFactors = FALSE
)
panel$date <- as.Date(panel$date)
site_ids <- setdiff(names(panel), "date")
X <- as.matrix(panel[, site_ids, drop = FALSE])
storage.mode(X) <- "numeric"
if (any(!is.finite(X))) stop("The unrestricted VAR diagnostic requires a complete panel.")

split <- analysis_split_2024(panel$date)
#build VAR response and lagged design matrices
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
#fit an unrestricted VAR by OLS
fit_ols_var <- function(X, p) {
  design <- var_design(X, p)
  decomposition <- qr(design$Z)
  fit <- lm.fit(design$Z, design$Y)
  coefficients <- fit$coefficients
  if (any(!is.finite(coefficients))) {
    stop("Unrestricted VAR coefficients are not uniquely estimable.")
  }

  singular_values <- svd(design$Z, nu = 0, nv = 0)$d
  list(
    coefficients = coefficients,
    p = p,
    time_rows = nrow(design$Z),
    predictors_per_equation = ncol(design$Z),
    total_coefficients = ncol(design$Y) * ncol(design$Z),
    design_rank = decomposition$rank,
    full_column_rank = decomposition$rank == ncol(design$Z),
    condition_number = max(singular_values) / min(singular_values),
    residual_variance = colSums(fit$residuals^2) /
      max(nrow(design$Z) - ncol(design$Z), 1)
  )
}
#one-day-ahead VAR forecasts
forecast_var <- function(X, fit, origins) {
  output <- matrix(
    NA_real_,
    length(origins),
    ncol(X),
    dimnames = list(NULL, colnames(X))
  )
  for (index in seq_along(origins)) {
    origin <- origins[index]
    predictor <- 1
    for (lag in seq_len(fit$p)) {
      predictor <- c(predictor, X[origin - lag + 1, ])
    }
    output[index, ] <- as.numeric(predictor %*% fit$coefficients)
  }
  output
}
#compare VAR(1) and VAR(2) on validation period
orders <- 1:2
validation_rows <- lapply(orders, function(p) {
  fit <- fit_ols_var(
    X[seq_len(split$selection_train_end), , drop = FALSE],
    p
  )
  forecast <- forecast_var(X, fit, split$validation_origins)
  observed <- X[split$validation_origins + 1, , drop = FALSE]
  data.frame(
    model = "Unrestricted OLS VAR",
    p = p,
    time_rows = fit$time_rows,
    predictors_per_equation = fit$predictors_per_equation,
    total_coefficients = fit$total_coefficients,
    design_rank = fit$design_rank,
    full_column_rank = fit$full_column_rank,
    condition_number = fit$condition_number,
    validation_rmse = rmse(observed - forecast),
    validation_mae = mae(observed - forecast),
    stringsAsFactors = FALSE
  )
})
validation <- do.call(rbind, validation_rows)
validation <- validation[order(validation$validation_rmse, validation$p), ]
validation$selected <- seq_len(nrow(validation)) == 1

write.csv(
  validation,
  file.path(table_dir, "unrestricted_var_validation_2024.csv"),
  row.names = FALSE
)
#select the best VAR order and refit using all train+val data
selected_p <- validation$p[1]
final_fit <- fit_ols_var(
  X[seq_len(split$fit_end), , drop = FALSE],
  selected_p
)
#evaluate the selected VAR on test period
forecast <- forecast_var(X, final_fit, split$test_origins)
observed <- X[split$test_origins + 1, , drop = FALSE]
metrics <- forecast_metrics(observed, forecast)

coefficient_norms <- apply(final_fit$coefficients[-1, , drop = FALSE], 2, function(x) {
  sqrt(sum(x^2))
})
maximum_absolute_coefficient <- apply(
  abs(final_fit$coefficients[-1, , drop = FALSE]),
  2,
  max
)
#save results
summary_table <- data.frame(
  model = "Unrestricted OLS VAR",
  model_order = paste0("VAR(", selected_p, ")"),
  selection_rule = "minimum chronological validation RMSE",
  time_rows = final_fit$time_rows,
  predictors_per_equation = final_fit$predictors_per_equation,
  total_coefficients = final_fit$total_coefficients,
  design_rank = final_fit$design_rank,
  full_column_rank = final_fit$full_column_rank,
  condition_number = final_fit$condition_number,
  median_coefficient_l2_norm = median(coefficient_norms),
  maximum_absolute_coefficient = max(maximum_absolute_coefficient),
  test_rmse = metrics$rmse,
  test_mae = metrics$mae,
  note = paste(
    "OLS VAR is estimable because each station equation has fewer",
    "predictors than time rows, but it remains parameter-rich and",
    "requires complete multivariate lag vectors."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(table_dir, "unrestricted_var_2024.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    date = rep(panel$date[split$test_origins + 1], each = ncol(X)),
    site_id = rep(site_ids, times = length(split$test_origins)),
    observed = as.vector(t(observed)),
    forecast = as.vector(t(forecast))
  ),
  file.path(table_dir, "unrestricted_var_forecasts_2024.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    station = site_ids,
    residual_variance = final_fit$residual_variance,
    coefficient_l2_norm = coefficient_norms,
    maximum_absolute_coefficient = maximum_absolute_coefficient,
    stringsAsFactors = FALSE
  ),
  file.path(table_dir, "unrestricted_var_stability_2024.csv"),
  row.names = FALSE
)

comparison_file <- file.path(table_dir, "var_benchmark_2024.csv")
if (file.exists(comparison_file)) {
  ridge <- read.csv(comparison_file, stringsAsFactors = FALSE)
  comparison <- data.frame(
    model = c("Unrestricted OLS VAR", "Ridge VAR"),
    model_order = c(summary_table$model_order, ridge$model_order[1]),
    regularisation = c("None", paste0("Ridge lambda=", ridge$selected_lambda[1])),
    total_coefficients = c(
      summary_table$total_coefficients,
      ridge$parameters[1]
    ),
    test_rmse = c(summary_table$test_rmse, ridge$test_rmse[1]),
    test_mae = c(summary_table$test_mae, ridge$test_mae[1]),
    stringsAsFactors = FALSE
  )
  comparison <- comparison[order(comparison$test_rmse), ]
  write.csv(
    comparison,
    file.path(table_dir, "ols_var_vs_ridge_var_2024.csv"),
    row.names = FALSE
  )
}
