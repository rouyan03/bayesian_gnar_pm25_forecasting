args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)

project_root <- if (length(script_arg) > 0) {
  script_file <- normalizePath(sub("^--file=", "", script_arg[1]))
  if (basename(dirname(script_file)) == "scripts") {
    dirname(dirname(script_file))
  } else {
    dirname(script_file)
  }
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

analysis_dir <- file.path(project_root, "data", "processed", "analysis")
network_dir <- file.path(project_root, "data", "processed", "networks")
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "outputs", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "R", "gnar_bayes_utils.R"))

#complete 2024 panel
panel <- read.csv(
  file.path(analysis_dir, "main_2024_pm25_daily_panel.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
panel$date <- as.Date(panel$date)
site_ids <- setdiff(names(panel), "date")
X_truth <- as.matrix(panel[, site_ids, drop = FALSE])
storage.mode(X_truth) <- "numeric"
if (any(!is.finite(X_truth))) {
  stop("This script expects the complete 2024 panel.")
}

split <- analysis_split_2024(panel$date)
observed_test <- X_truth[split$test_origins + 1, , drop = FALSE]

#selected GNAR specification (gnar(2,(1,1)))
p <- 2
stages <- c(1, 1)
A <- read_matrix_csv(
  file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv")
)
A <- row_normalise(A[site_ids, site_ids, drop = FALSE])

n_rep <- as.integer(Sys.getenv("MISSINGNESS_REPS", "10"))
missing_rates <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)
refit_missingness <- Sys.getenv("REFIT_MISSINGNESS", "0") == "1"
main_results_file <- file.path(table_dir, "missingness_main_model_replicates_2024.csv")
hierarchical_results_file <- file.path(table_dir, "missingness_hierarchical_replicates_2024.csv")
ar1_results_file <- file.path(table_dir, "station_ar1_extended_missingness_results_2024.csv")

#fill training gaps by linear interpolation
fill_training <- function(X) {
  output <- X
  for (node in seq_len(ncol(X))) {
    observed <- which(is.finite(X[, node]))
    if (length(observed) < 2) stop("Insufficient observations for ", colnames(X)[node])
    output[, node] <- approx(
      observed,
      X[observed, node],
      seq_len(nrow(X)),
      rule = 2
    )$y
  }
  output
}

#same training observations missing each time
make_missing_index <- function(missing_rate, replication) {
  set.seed(300000 + round(100 * missing_rate) * 1000 + replication)
  missing_index <- matrix(FALSE, nrow(X_truth), ncol(X_truth))
  if (missing_rate > 0) {
    missing_index[seq_len(split$fit_end), ] <- matrix(
      runif(split$fit_end * ncol(missing_index)) < missing_rate,
      nrow = split$fit_end,
      ncol = ncol(missing_index)
    )
  }
  missing_index[seq_len(p), ] <- FALSE
  missing_index[(split$fit_end + 1):nrow(missing_index), ] <- FALSE
  missing_index
}

available_neighbour <- function(values, weights) {
  available <- is.finite(values) & weights > 0
  if (!any(available)) return(NA_real_)
  sum(weights[available] * values[available]) / sum(weights[available])
}

build_missing_gnar_design <- function(X, station_intercepts = FALSE) {
  stage_matrices <- lapply(seq_len(max(stages)), function(stage) stage_matrix(A, stage))
  rows <- list()
  response <- numeric()
  station_id <- character()

  for (time in (p + 1):nrow(X)) {
    for (node in seq_len(ncol(X))) {
      if (!is.finite(X[time, node])) next
      own_lags <- X[time - seq_len(p), node]
      if (any(!is.finite(own_lags))) next

      row <- if (station_intercepts) numeric() else c(intercept = 1)
      valid <- TRUE
      for (lag in seq_len(p)) {
        row[paste0("alpha_lag", lag)] <- own_lags[lag]
        for (stage in seq_len(stages[lag])) {
          neighbour <- available_neighbour(
            X[time - lag, ],
            stage_matrices[[stage]][node, ]
          )
          if (!is.finite(neighbour)) {
            valid <- FALSE
            break
          }
          row[paste0("beta_lag", lag, "_stage", stage)] <- neighbour
        }
        if (!valid) break
      }

      if (valid) {
        rows[[length(rows) + 1]] <- row
        response <- c(response, X[time, node])
        station_id <- c(station_id, colnames(X)[node])
      }
    }
  }

  if (length(rows) == 0) stop("No valid GNAR rows remain.")
  R <- do.call(rbind, rows)
  if (station_intercepts) {
    station <- factor(station_id, levels = colnames(X))
    Z <- model.matrix(~station - 1)
    colnames(Z) <- paste0("station_", colnames(X))
    R <- cbind(R, Z)
  }
  storage.mode(R) <- "numeric"
  list(y = response, R = R)
}

#one step GNAR forecast
predict_missing_gnar <- function(history, theta) {
  stage_matrices <- lapply(seq_len(max(stages)), function(stage) stage_matrix(A, stage))
  forecast <- rep(if ("intercept" %in% names(theta)) theta["intercept"] else 0, ncol(history))
  fallback <- mean(history[nrow(history), ], na.rm = TRUE)

  for (lag in seq_len(p)) {
    values <- history[nrow(history) - lag + 1, ]
    own <- values
    own[!is.finite(own)] <- fallback
    forecast <- forecast + theta[paste0("alpha_lag", lag)] * own

    for (stage in seq_len(stages[lag])) {
      neighbour <- vapply(seq_len(ncol(history)), function(node) {
        available_neighbour(values, stage_matrices[[stage]][node, ])
      }, numeric(1))
      neighbour[!is.finite(neighbour)] <- fallback
      forecast <- forecast + theta[paste0("beta_lag", lag, "_stage", stage)] * neighbour
    }
  }
  forecast
}

forecast_gnar_missing <- function(X_missing, theta, draws = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  history <- X_missing
  history[seq_len(split$fit_end), ] <- fill_training(
    X_missing[seq_len(split$fit_end), , drop = FALSE]
  )
  point <- lower <- upper <- matrix(
    NA_real_,
    length(split$test_origins),
    ncol(X_missing),
    dimnames = list(NULL, colnames(X_missing))
  )

  for (index in seq_along(split$test_origins)) {
    origin <- split$test_origins[index]
    if (is.null(draws)) {
      point[index, ] <- predict_missing_gnar(history[seq_len(origin), , drop = FALSE], theta)
    } else {
      predictors <- matrix(NA_real_, nrow(draws$theta), ncol(X_missing))
      for (draw in seq_len(nrow(draws$theta))) {
        predictors[draw, ] <- predict_missing_gnar(
          history[seq_len(origin), , drop = FALSE],
          draws$theta[draw, ]
        )
      }
      predictive <- predictors + matrix(
        rnorm(length(predictors)) * rep(sqrt(draws$sigma2), ncol(X_missing)),
        nrow = nrow(predictors)
      )
      point[index, ] <- colMeans(predictors)
      lower[index, ] <- apply(predictive, 2, quantile, 0.025)
      upper[index, ] <- apply(predictive, 2, quantile, 0.975)
    }

    next_time <- origin + 1
    available <- is.finite(X_missing[next_time, ])
    history[next_time, ] <- point[index, ]
    history[next_time, available] <- X_missing[next_time, available]
  }
  list(mean = point, lower = lower, upper = upper)
}

#unrestricted VAR(1) 
fit_var1 <- function(X) {
  Y <- X[2:nrow(X), , drop = FALSE]
  Z <- cbind(intercept = 1, X[1:(nrow(X) - 1), , drop = FALSE])
  fit <- lm.fit(Z, Y)
  coefficients <- fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  list(coefficients = coefficients, p = 1)
}

forecast_var_missing <- function(X_missing, fit) {
  history <- X_missing
  history[seq_len(split$fit_end), ] <- fill_training(
    X_missing[seq_len(split$fit_end), , drop = FALSE]
  )
  output <- matrix(
    NA_real_,
    length(split$test_origins),
    ncol(X_missing),
    dimnames = list(NULL, colnames(X_missing))
  )
  for (index in seq_along(split$test_origins)) {
    origin <- split$test_origins[index]
    output[index, ] <- as.numeric(c(1, history[origin, ]) %*% fit$coefficients)
    available <- is.finite(X_missing[origin + 1, ])
    history[origin + 1, ] <- output[index, ]
    history[origin + 1, available] <- X_missing[origin + 1, available]
  }
  output
}

#hierarchical Bayesian GNAR 
fit_hierarchical_missing_chain <- function(y, R, n_iter = 1800, burn = 500, seed = 1) {
  set.seed(seed)
  station_columns <- paste0("station_", site_ids)
  station_position <- match(station_columns, colnames(R))
  common_position <- setdiff(seq_len(ncol(R)), station_position)
  d <- length(site_ids)
  n <- length(y)

  theta <- lm.fit(R, y)$coefficients
  theta[!is.finite(theta)] <- 0
  names(theta) <- colnames(R)
  sigma2 <- var(y - as.vector(R %*% theta))
  mu0 <- mean(theta[station_position])
  tau2_mu <- var(theta[station_position])
  if (!is.finite(tau2_mu) || tau2_mu <= 0) tau2_mu <- 1

  kept <- n_iter - burn
  theta_draws <- matrix(NA_real_, kept, ncol(R), dimnames = list(NULL, colnames(R)))
  sigma2_draws <- tau2_draws <- mu0_draws <- numeric(kept)
  XtX <- crossprod(R)
  Xty <- crossprod(R, y)
  a0 <- b0 <- 0.01
  mu0_prior_variance <- 100^2

  for (iteration in seq_len(n_iter)) {
    prior_precision <- numeric(ncol(R))
    prior_mean <- numeric(ncol(R))
    prior_precision[common_position] <- 1
    prior_precision[station_position] <- 1 / tau2_mu
    prior_mean[station_position] <- mu0

    precision <- XtX / sigma2 + diag(prior_precision, ncol(R))
    rhs <- as.vector(Xty) / sigma2 + prior_precision * prior_mean
    theta <- rmvnorm_precision(as.vector(solve(precision, rhs)), precision)
    names(theta) <- colnames(R)

    station_intercepts <- theta[station_position]
    mu0_variance <- 1 / (d / tau2_mu + 1 / mu0_prior_variance)
    mu0_mean <- mu0_variance * sum(station_intercepts) / tau2_mu
    mu0 <- rnorm(1, mu0_mean, sqrt(mu0_variance))
    tau2_mu <- rinvgamma(
      1,
      shape = a0 + d / 2,
      rate = b0 + sum((station_intercepts - mu0)^2) / 2
    )
    sigma2 <- rinvgamma(
      1,
      shape = a0 + n / 2,
      rate = b0 + sum((y - as.vector(R %*% theta))^2) / 2
    )

    if (iteration > burn) {
      index <- iteration - burn
      theta_draws[index, ] <- theta
      sigma2_draws[index] <- sigma2
      tau2_draws[index] <- tau2_mu
      mu0_draws[index] <- mu0
    }
  }
  list(theta = theta_draws, sigma2 = sigma2_draws, tau2_mu = tau2_draws, mu0 = mu0_draws)
}

forecast_hierarchical_missing <- function(X_missing, draws) {
  history <- X_missing
  history[seq_len(split$fit_end), ] <- fill_training(
    X_missing[seq_len(split$fit_end), , drop = FALSE]
  )
  keep <- unique(round(seq(1, nrow(draws$theta), length.out = min(500, nrow(draws$theta)))))
  point <- matrix(NA_real_, length(split$test_origins), ncol(X_missing), dimnames = list(NULL, site_ids))

  for (index in seq_along(split$test_origins)) {
    origin <- split$test_origins[index]
    predictors <- matrix(NA_real_, length(keep), ncol(X_missing))
    for (draw_index in seq_along(keep)) {
      theta <- draws$theta[keep[draw_index], ]
      lag1 <- history[origin, ]
      lag2 <- history[origin - 1, ]
      fallback <- mean(lag1, na.rm = TRUE)
      lag1[!is.finite(lag1)] <- fallback
      lag2[!is.finite(lag2)] <- fallback
      predictors[draw_index, ] <- as.numeric(theta[paste0("station_", site_ids)]) +
        theta["alpha_lag1"] * lag1 +
        theta["beta_lag1_stage1"] * as.numeric(A %*% lag1) +
        theta["alpha_lag2"] * lag2 +
        theta["beta_lag2_stage1"] * as.numeric(A %*% lag2)
    }
    point[index, ] <- colMeans(predictors)
    available <- is.finite(X_missing[origin + 1, ])
    history[origin + 1, ] <- point[index, ]
    history[origin + 1, available] <- X_missing[origin + 1, available]
  }
  point
}

#refit the standard GNAR and VAR models for each missingness level
run_main_missingness_models <- function() {
  rows <- list()
  row_id <- 1
  bayes_seed <- 1

  for (missing_rate in missing_rates) {
    for (replication in seq_len(n_rep)) {
      message("Main missingness: ", 100 * missing_rate, "%, replication ", replication, "/", n_rep)
      missing_index <- make_missing_index(missing_rate, replication)
      X_missing <- X_truth
      X_missing[missing_index] <- NA_real_

      design <- build_missing_gnar_design(X_missing[seq_len(split$fit_end), , drop = FALSE])
      ols_fit <- lm.fit(design$R, design$y)
      theta_ols <- ols_fit$coefficients
      theta_ols[is.na(theta_ols)] <- 0
      names(theta_ols) <- colnames(design$R)

      bayes <- fit_gaussian_chain(
        design$y,
        design$R,
        coefficient_tau2 = 1,
        intercept_tau2 = 100^2,
        n_iter = 1800,
        burn = 500,
        seed = 305000 + bayes_seed
      )
      keep <- unique(round(seq(1, nrow(bayes$theta), length.out = min(500, nrow(bayes$theta)))))
      bayes_forecast <- forecast_gnar_missing(
        X_missing,
        colMeans(bayes$theta),
        list(theta = bayes$theta[keep, , drop = FALSE], sigma2 = bayes$sigma2[keep]),
        seed = 306000 + round(100 * missing_rate) * 100 + replication
      )

      training <- fill_training(X_missing[seq_len(split$fit_end), , drop = FALSE])
      var_forecast <- forecast_var_missing(X_missing, fit_var1(training))
      ols_gnar_forecast <- forecast_gnar_missing(X_missing, theta_ols)$mean

      forecasts <- list(
        `Unrestricted VAR` = var_forecast,
        `OLS GNAR` = ols_gnar_forecast,
        `Bayesian Gaussian GNAR` = bayes_forecast$mean
      )

      for (model in names(forecasts)) {
        rows[[row_id]] <- data.frame(
          missing_rate = missing_rate,
          missing_percent = 100 * missing_rate,
          replication = replication,
          model = model,
          realised_training_missing_fraction = mean(missing_index[seq_len(split$fit_end), , drop = FALSE]),
          training_rows = length(design$y),
          rmse = rmse(observed_test - forecasts[[model]]),
          mae = mae(observed_test - forecasts[[model]]),
          coverage_95 = if (model == "Bayesian Gaussian GNAR") {
            mean(observed_test >= bayes_forecast$lower & observed_test <= bayes_forecast$upper)
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1
      }
      bayes_seed <- bayes_seed + 4
    }
  }
  do.call(rbind, rows)
}

#refit the hierarchical model separately 
run_hierarchical_missingness_model <- function() {
  rows <- list()
  row_id <- 1

  for (missing_rate in missing_rates) {
    for (replication in seq_len(n_rep)) {
      message("Hierarchical missingness: ", 100 * missing_rate, "%, replication ", replication, "/", n_rep)
      missing_index <- make_missing_index(missing_rate, replication)
      X_missing <- X_truth
      X_missing[missing_index] <- NA_real_

      design <- build_missing_gnar_design(
        X_missing[seq_len(split$fit_end), , drop = FALSE],
        station_intercepts = TRUE
      )
      fit <- fit_hierarchical_missing_chain(
        design$y,
        design$R,
        n_iter = 1800,
        burn = 500,
        seed = 540000 + round(100 * missing_rate) * 100 + replication
      )
      forecast <- forecast_hierarchical_missing(X_missing, fit)

      rows[[row_id]] <- data.frame(
        missing_rate = missing_rate,
        missing_percent = 100 * missing_rate,
        replication = replication,
        model = "Hierarchical Bayesian GNAR",
        realised_training_missing_fraction = mean(missing_index[seq_len(split$fit_end), , drop = FALSE]),
        training_rows = length(design$y),
        rmse = rmse(observed_test - forecast),
        mae = mae(observed_test - forecast),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }
  }
  do.call(rbind, rows)
}

#separate AR(1) model for each station
fit_station_ar1 <- function(X_train) {
  coefficients <- matrix(
    NA_real_,
    nrow = ncol(X_train),
    ncol = 2,
    dimnames = list(colnames(X_train), c("intercept", "ar1"))
  )
  for (site in colnames(X_train)) {
    y <- X_train[2:nrow(X_train), site]
    lag1 <- X_train[1:(nrow(X_train) - 1), site]
    keep <- is.finite(y) & is.finite(lag1)
    fit <- lm.fit(cbind(intercept = 1, ar1 = lag1[keep]), y[keep])
    coefficients[site, ] <- ifelse(is.na(fit$coefficients), 0, fit$coefficients)
  }
  coefficients
}

forecast_ar1_missing <- function(X_missing, coefficients) {
  history <- X_missing
  history[seq_len(split$fit_end), ] <- fill_training(
    X_missing[seq_len(split$fit_end), , drop = FALSE]
  )
  output <- matrix(
    NA_real_,
    length(split$test_origins),
    ncol(X_missing),
    dimnames = list(NULL, colnames(X_missing))
  )
  for (index in seq_along(split$test_origins)) {
    origin <- split$test_origins[index]
    output[index, ] <- coefficients[, "intercept"] + coefficients[, "ar1"] * history[origin, ]
    available <- is.finite(X_missing[origin + 1, ])
    history[origin + 1, ] <- output[index, ]
    history[origin + 1, available] <- X_missing[origin + 1, available]
  }
  output
}

run_ar1_missingness_model <- function() {
  rows <- list()
  row_id <- 1

  for (missing_rate in missing_rates) {
    for (replication in seq_len(n_rep)) {
      message("AR(1) missingness: ", 100 * missing_rate, "%, replication ", replication, "/", n_rep)
      missing_index <- make_missing_index(missing_rate, replication)
      X_missing <- X_truth
      X_missing[missing_index] <- NA_real_

      fit <- fit_station_ar1(X_missing[seq_len(split$fit_end), , drop = FALSE])
      forecast <- forecast_ar1_missing(X_missing, fit)

      rows[[row_id]] <- data.frame(
        missing_rate = missing_rate,
        missing_percent = 100 * missing_rate,
        replication = replication,
        model = "AR(1)",
        realised_training_missing_fraction = mean(missing_index[seq_len(split$fit_end), , drop = FALSE]),
        training_rows = sum(complete.cases(cbind(
          as.vector(X_missing[2:split$fit_end, , drop = FALSE]),
          as.vector(X_missing[1:(split$fit_end - 1), , drop = FALSE])
        ))),
        rmse = rmse(observed_test - forecast),
        mae = mae(observed_test - forecast),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }
  }
  do.call(rbind, rows)
}

#average over repeated random missingness runs
summarise_missingness <- function(results) {
  mean_summary <- aggregate(
    cbind(rmse, mae, realised_training_missing_fraction) ~
      missing_rate + missing_percent + model,
    data = results,
    FUN = mean
  )
  rmse_uncertainty <- aggregate(
    rmse ~ missing_rate + missing_percent + model,
    data = results,
    FUN = function(values) {
      c(
        mcse = sd(values) / sqrt(length(values)),
        lower = mean(values) - 1.96 * sd(values) / sqrt(length(values)),
        upper = mean(values) + 1.96 * sd(values) / sqrt(length(values))
      )
    }
  )
  parts <- if (is.matrix(rmse_uncertainty$rmse)) {
    as.data.frame(rmse_uncertainty$rmse)
  } else {
    as.data.frame(do.call(rbind, rmse_uncertainty$rmse))
  }
  out <- cbind(
    rmse_uncertainty[c("missing_rate", "missing_percent", "model")],
    parts
  )
  names(out)[4:6] <- c("rmse_mcse", "rmse_lower_95_mc", "rmse_upper_95_mc")
  merge(
    mean_summary,
    out,
    by = c("missing_rate", "missing_percent", "model"),
    all.x = TRUE
  )
}

#first run creates these files. Later runs read them from saved files
if (!refit_missingness && file.exists(main_results_file)) {
  main_results <- read.csv(main_results_file, stringsAsFactors = FALSE)
} else {
  main_results <- run_main_missingness_models()
  write.csv(main_results, main_results_file, row.names = FALSE)
}

if (!refit_missingness && file.exists(hierarchical_results_file)) {
  hierarchical_results <- read.csv(hierarchical_results_file, stringsAsFactors = FALSE)
} else {
  hierarchical_results <- run_hierarchical_missingness_model()
  write.csv(hierarchical_results, hierarchical_results_file, row.names = FALSE)
}

missingness_columns <- c(
  "missing_rate",
  "missing_percent",
  "replication",
  "model",
  "realised_training_missing_fraction",
  "rmse",
  "mae"
)
base_summary <- summarise_missingness(rbind(
  subset(
    main_results,
    model %in% c("Bayesian Gaussian GNAR", "Unrestricted VAR")
  )[, missingness_columns],
  hierarchical_results[, missingness_columns]
))
base_summary <- base_summary[order(base_summary$missing_percent, base_summary$model), ]
write.csv(main_results, file.path(table_dir, "extended_missingness_results_2024.csv"), row.names = FALSE)
write.csv(base_summary, file.path(table_dir, "three_model_missingness_summary_2024.csv"), row.names = FALSE)

if (!refit_missingness && file.exists(ar1_results_file)) {
  ar_results <- read.csv(ar1_results_file, stringsAsFactors = FALSE)
} else {
  ar_results <- run_ar1_missingness_model()
  write.csv(ar_results, ar1_results_file, row.names = FALSE)
}

ar_summary <- summarise_missingness(ar_results)
keep_columns <- c(
  "missing_rate",
  "missing_percent",
  "model",
  "rmse",
  "mae",
  "realised_training_missing_fraction",
  "rmse_mcse",
  "rmse_lower_95_mc",
  "rmse_upper_95_mc"
)

#final table 
comparison <- rbind(
  base_summary[
    base_summary$model %in% c(
      "Hierarchical Bayesian GNAR",
      "Bayesian Gaussian GNAR",
      "Unrestricted VAR"
    ),
    keep_columns
  ],
  ar_summary[, keep_columns]
)
comparison <- comparison[order(comparison$missing_percent, comparison$model), ]
write.csv(
  comparison,
  file.path(table_dir, "hierarchical_missingness_summary_with_ar1_2024.csv"),
  row.names = FALSE
)

comparison$model <- factor(
  comparison$model,
  levels = c(
    "Hierarchical Bayesian GNAR",
    "Bayesian Gaussian GNAR",
    "AR(1)",
    "Unrestricted VAR"
  )
)
comparison_40 <- subset(comparison, missing_percent <= 40)

#plot figure
plot_40 <- ggplot2::ggplot(
  comparison_40,
  ggplot2::aes(
    missing_percent,
    rmse,
    colour = model,
    linetype = model,
    shape = model
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = rmse_lower_95_mc,
      ymax = rmse_upper_95_mc
    ),
    width = 1.3,
    linewidth = 0.48,
    alpha = 0.7
  ) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::scale_colour_manual(
    values = c(
      "Hierarchical Bayesian GNAR" = "#1D5D7C",
      "Bayesian Gaussian GNAR" = "#168B67",
      "AR(1)" = "#5B5B5B",
      "Unrestricted VAR" = "#A33C31"
    )
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "Hierarchical Bayesian GNAR" = "solid",
      "Bayesian Gaussian GNAR" = "dashed",
      "AR(1)" = "longdash",
      "Unrestricted VAR" = "dotdash"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Hierarchical Bayesian GNAR" = 17,
      "Bayesian Gaussian GNAR" = 16,
      "AR(1)" = 18,
      "Unrestricted VAR" = 15
    )
  ) +
  ggplot2::scale_x_continuous(
    breaks = 100 * missing_rates,
    labels = paste0(100 * missing_rates, "%")
  ) +
  ggplot2::labs(
    x = "Training observations removed (%)",
    y = "Test RMSPE",
    colour = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(1.35, "cm"),
    legend.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 15),
    axis.text = ggplot2::element_text(size = 12, colour = "#333333")
  )

ggplot2::ggsave(
  file.path(figure_dir, "poster_missingness_robustness_compact_with_ar1_2024.png"),
  plot_40,
  width = 9.2,
  height = 5.6,
  dpi = 320
)
