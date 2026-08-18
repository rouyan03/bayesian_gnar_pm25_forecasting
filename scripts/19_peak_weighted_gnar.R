project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(project_root, "data", "processed", "analysis")
network_dir <- file.path(project_root, "data", "processed", "networks")
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "outputs", "figures")
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
if (any(!is.finite(X))) stop("Peak-focused extensions require the complete 2024 panel.")

split <- analysis_split_2024(panel$date)
A <- read_matrix_csv(
  file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv")
)
A <- row_normalise(A[site_ids, site_ids, drop = FALSE])
p <- 2
stages <- c(1, 1)
W <- stage_matrix(A, 1)

#design matrix for equal-neighbour GNAR(2,[1,1])
origin_design <- function(origin) {
  lag1 <- as.numeric(X[origin, ])
  lag2 <- as.numeric(X[origin - 1, ])
  cbind(
    intercept = 1,
    alpha_lag1 = lag1,
    beta_lag1_stage1 = as.numeric(W %*% lag1),
    alpha_lag2 = lag2,
    beta_lag2_stage1 = as.numeric(W %*% lag2)
  )
}

forecast_coefficients <- function(coefficients, origins) {
  output <- matrix(
    NA_real_,
    length(origins),
    ncol(X),
    dimnames = list(NULL, site_ids)
  )
  for (index in seq_along(origins)) {
    output[index, ] <- as.numeric(origin_design(origins[index]) %*% coefficients)
  }
  output
}

#overall scores and highest decile scores for a point forecast
score_point_forecast <- function(observed, forecast, label, period) {
  observed_vector <- as.vector(observed)
  forecast_vector <- as.vector(forecast)
  cutoff <- as.numeric(quantile(observed_vector, 0.90, names = FALSE))
  peak <- observed_vector >= cutoff
  data.frame(
    model = label,
    period = period,
    rmse = rmse(observed_vector - forecast_vector),
    mae = mae(observed_vector - forecast_vector),
    mean_bias_forecast_minus_observed = mean(forecast_vector - observed_vector),
    peak_threshold = cutoff,
    peak_count = sum(peak),
    peak_rmse = rmse(observed_vector[peak] - forecast_vector[peak]),
    peak_mae = mae(observed_vector[peak] - forecast_vector[peak]),
    peak_bias_forecast_minus_observed =
      mean(forecast_vector[peak] - observed_vector[peak]),
    stringsAsFactors = FALSE
  )
}

#peak-weighted least squares gives high observations more influence in fitting
fit_peak_weighted <- function(R, y, quantile_level, peak_weight) {
  threshold <- as.numeric(quantile(y, quantile_level, names = FALSE))
  weights <- ifelse(y >= threshold, peak_weight, 1)
  coefficients <- lm.wfit(R, y, w = weights)$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  setNames(coefficients, colnames(R))
}

selection_design <- build_gnar_design(
  X[seq_len(split$selection_train_end), , drop = FALSE],
  A,
  p,
  stages
)
final_design <- build_gnar_design(
  X[seq_len(split$fit_end), , drop = FALSE],
  A,
  p,
  stages
)
observed_validation <- X[split$validation_origins + 1, , drop = FALSE]
observed_test <- X[split$test_origins + 1, , drop = FALSE]

baseline_selection <- setNames(
  lm.fit(selection_design$R, selection_design$y)$coefficients,
  colnames(selection_design$R)
)
baseline_final <- setNames(
  lm.fit(final_design$R, final_design$y)$coefficients,
  colnames(final_design$R)
)

validation_rows <- list()
counter <- 1

add_validation <- function(label, forecast, family, setting) {
  score <- score_point_forecast(
    observed_validation,
    forecast,
    label,
    "validation"
  )
  score$family <- family
  score$setting <- setting
  validation_rows[[counter]] <<- score
  counter <<- counter + 1
}

baseline_validation <- forecast_coefficients(
  baseline_selection,
  split$validation_origins
)
add_validation("Gaussian GNAR baseline", baseline_validation, "baseline", "OLS mean")

#try a small grid of peak thresholds and fitting weights (hyperparam tuning)
for (threshold_quantile in c(0.80, 0.90, 0.95)) {
  for (peak_weight in c(2, 4, 8, 12)) {
    coefficients <- fit_peak_weighted(
      selection_design$R,
      selection_design$y,
      threshold_quantile,
      peak_weight
    )
    label <- paste0(
      "Peak-weighted GNAR q=", threshold_quantile,
      " weight=", peak_weight
    )
    add_validation(
      label,
      forecast_coefficients(coefficients, split$validation_origins),
      "peak_weighted",
      paste0("q=", threshold_quantile, ";weight=", peak_weight)
    )
  }
}

validation <- do.call(rbind, validation_rows)
baseline_validation_rmse <- validation$rmse[validation$family == "baseline"]

#wanna keep average RMSPE close to baseline while reducing peak bias (5%)
validation$rmse_guard <- validation$rmse <= baseline_validation_rmse * 1.05
validation$selection_score <-
  abs(validation$peak_bias_forecast_minus_observed) +
  2 * pmax(0, validation$rmse / baseline_validation_rmse - 1)
validation <- validation[
  order(!validation$rmse_guard, validation$selection_score, validation$rmse),
]

write.csv(
  validation,
  file.path(table_dir, "peak_focused_gnar_validation_2024.csv"),
  row.names = FALSE
)

selected_peak_weighted <- validation[
  validation$family == "peak_weighted" & validation$rmse_guard,
  ,
  drop = FALSE
][1, , drop = FALSE]

#refit the selected peak-weighted model on the full 2024 training period (weight=4, q=0.95)
test_rows <- list()
test_forecasts <- list()
test_counter <- 1

add_test <- function(label, forecast, family, setting) {
  score <- score_point_forecast(observed_test, forecast, label, "test")
  score$family <- family
  score$setting <- setting
  test_rows[[test_counter]] <<- score
  test_forecasts[[label]] <<- forecast
  test_counter <<- test_counter + 1
}

baseline_test <- forecast_coefficients(baseline_final, split$test_origins)
add_test("Gaussian GNAR baseline", baseline_test, "baseline", "OLS mean")

#apply it to the test period
parts <- strsplit(selected_peak_weighted$setting, ";", fixed = TRUE)[[1]]
threshold_quantile <- as.numeric(sub("q=", "", parts[1]))
peak_weight <- as.numeric(sub("weight=", "", parts[2]))
peak_weighted_coefficients <- fit_peak_weighted(
  final_design$R,
  final_design$y,
  threshold_quantile,
  peak_weight
)
peak_weighted_forecast <- forecast_coefficients(
  peak_weighted_coefficients,
  split$test_origins
)
add_test(
  selected_peak_weighted$model,
  peak_weighted_forecast,
  "peak_weighted",
  selected_peak_weighted$setting
)

test <- do.call(rbind, test_rows)
baseline_test_bias <- abs(
  test$peak_bias_forecast_minus_observed[test$family == "baseline"]
)
baseline_test_rmse <- test$rmse[test$family == "baseline"]
test$peak_bias_reduction_fraction <-
  1 - abs(test$peak_bias_forecast_minus_observed) / baseline_test_bias
test$success_peak_bias <- test$peak_bias_reduction_fraction >= 0.50
test$success_rmse_guard <- test$rmse <= baseline_test_rmse * 1.05
test$success <- test$success_peak_bias & test$success_rmse_guard

write.csv(
  test,
  file.path(table_dir, "peak_focused_gnar_test_2024.csv"),
  row.names = FALSE
)

forecast_output <- do.call(rbind, lapply(names(test_forecasts), function(label) {
  forecast <- test_forecasts[[label]]
  data.frame(
    model = label,
    date = rep(panel$date[split$test_origins + 1], each = length(site_ids)),
    site_id = rep(site_ids, times = length(split$test_origins)),
    observed = as.vector(t(observed_test)),
    forecast = as.vector(t(forecast)),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  forecast_output,
  file.path(table_dir, "peak_focused_gnar_forecasts_2024.csv"),
  row.names = FALSE
)

#table values reported in the upper-tail subsection
peak_summary <- function(model, observed, forecast, lower = NULL, upper = NULL, selected_nu = NA_real_) {
  observed_vector <- as.vector(observed)
  forecast_vector <- as.vector(forecast)
  peak_threshold <- as.numeric(quantile(observed_vector, 0.90, names = FALSE))
  peak <- observed_vector >= peak_threshold

  data.frame(
    model = model,
    rmspe = rmse(observed_vector - forecast_vector),
    mae = mae(observed_vector - forecast_vector),
    peak_bias = mean(forecast_vector[peak] - observed_vector[peak]),
    peak_rmspe = rmse(observed_vector[peak] - forecast_vector[peak]),
    peak_coverage = if (is.null(lower) || is.null(upper)) {
      NA_real_
    } else {
      lower_vector <- as.vector(lower)
      upper_vector <- as.vector(upper)
      mean(observed_vector[peak] >= lower_vector[peak] & observed_vector[peak] <= upper_vector[peak])
    },
    peak_threshold = peak_threshold,
    peak_count = sum(peak),
    selected_nu = selected_nu,
    stringsAsFactors = FALSE
  )
}

#mean forecast error by observed PM2.5 decile
decile_errors <- function(model, observed, forecast) {
  observed_vector <- as.vector(observed)
  forecast_vector <- as.vector(forecast)
  decile <- cut(
    observed_vector,
    breaks = quantile(observed_vector, seq(0, 1, 0.1), names = FALSE),
    include.lowest = TRUE,
    labels = FALSE
  )
  aggregate(
    forecast_vector - observed_vector,
    by = list(model = rep(model, length(decile)), decile = decile),
    FUN = mean
  ) |>
    setNames(c("model", "decile", "forecast_error"))
}

hierarchical_file <- file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_forecasts_2024.csv")
student_t_file <- file.path(table_dir, "bayes_student_t_gnar_forecasts_2024.csv")
student_t_summary_file <- file.path(table_dir, "bayes_student_t_gnar_2024.csv")
if (!file.exists(hierarchical_file) || !file.exists(student_t_file)) {
  stop("Run scripts 13 and 18 before script 19.")
}

hierarchical <- read.csv(hierarchical_file, stringsAsFactors = FALSE)
student_t <- read.csv(student_t_file, stringsAsFactors = FALSE)
student_t_summary <- read.csv(student_t_summary_file, stringsAsFactors = FALSE)
peak_weighted_row <- test[test$family == "peak_weighted", ][1, ]
peak_weighted_forecast <- test_forecasts[[peak_weighted_row$model]]

peak_comparison <- rbind(
  peak_summary(
    "Hierarchical Bayesian GNAR",
    observed_test,
    matrix(hierarchical$forecast_mean, ncol = length(site_ids), byrow = TRUE),
    matrix(hierarchical$lower_95, ncol = length(site_ids), byrow = TRUE),
    matrix(hierarchical$upper_95, ncol = length(site_ids), byrow = TRUE)
  ),
  peak_summary(
    "Student-t Bayesian GNAR",
    observed_test,
    matrix(student_t$forecast_mean, ncol = length(site_ids), byrow = TRUE),
    matrix(student_t$lower_95, ncol = length(site_ids), byrow = TRUE),
    matrix(student_t$upper_95, ncol = length(site_ids), byrow = TRUE),
    selected_nu = student_t_summary$selected_nu[1]
  ),
  peak_summary(
    "Peak-weighted GNAR",
    observed_test,
    peak_weighted_forecast
  )
)
write.csv(
  peak_comparison,
  file.path(table_dir, "peak_student_t_weighted_comparison_2024.csv"),
  row.names = FALSE
)

decile_table <- rbind(
  decile_errors(
    "Hierarchical Bayesian GNAR",
    observed_test,
    matrix(hierarchical$forecast_mean, ncol = length(site_ids), byrow = TRUE)
  ),
  decile_errors(
    "Student-t Bayesian GNAR",
    observed_test,
    matrix(student_t$forecast_mean, ncol = length(site_ids), byrow = TRUE)
  ),
  decile_errors("Peak-weighted GNAR", observed_test, peak_weighted_forecast)
)
write.csv(
  decile_table,
  file.path(table_dir, "peak_models_decile_forecast_error_comparison_2024.csv"),
  row.names = FALSE
)

#plot figure to compare upper-tail bias across the three models.
if (requireNamespace("ggplot2", quietly = TRUE)) {
    decile_table$model <- factor(
      decile_table$model,
      levels = c(
        "Hierarchical Bayesian GNAR",
        "Peak-weighted GNAR",
        "Student-t Bayesian GNAR"
      )
    )

    peak_plot <- ggplot2::ggplot(
      decile_table,
      ggplot2::aes(
        decile,
        forecast_error,
        colour = model,
        linetype = model,
        shape = model
      )
    ) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.35) +
      ggplot2::geom_line(linewidth = 0.85) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = 1:10) +
      ggplot2::scale_colour_manual(
        values = c(
          "Hierarchical Bayesian GNAR" = "blue",
          "Peak-weighted GNAR" = "green3",
          "Student-t Bayesian GNAR" = "red"
        )
      ) +
      ggplot2::scale_linetype_manual(
        values = c(
          "Hierarchical Bayesian GNAR" = "solid",
          "Peak-weighted GNAR" = "dotdash",
          "Student-t Bayesian GNAR" = "dashed"
        )
      ) +
      ggplot2::scale_shape_manual(
        values = c(
          "Hierarchical Bayesian GNAR" = 16,
          "Peak-weighted GNAR" = 15,
          "Student-t Bayesian GNAR" = 17
        )
      ) +
      ggplot2::labs(
        x = "Observed PM2.5 decile",
        y = "Forecast error",
        colour = NULL,
        linetype = NULL,
        shape = NULL
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")

    ggplot2::ggsave(
      file.path(figure_dir, "peak_models_decile_forecast_error_comparison_2024.png"),
      peak_plot,
      width = 5.2,
      height = 3.1,
      dpi = 320
    )
}
