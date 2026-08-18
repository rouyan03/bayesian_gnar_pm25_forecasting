project_root <- normalizePath(getwd(), mustWork = TRUE)
processed_dir <- file.path(project_root, "data", "processed", "defra_aurn")
network_dir <- file.path(project_root, "data", "processed", "networks")
table_dir <- file.path(project_root, "outputs", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "R", "gnar_bayes_utils.R"))

#read one yearly daily PM2.5 panel
read_panel <- function(year) {
  panel <- read.csv(
    file.path(processed_dir, paste0("aurn_pm25_daily_panel_", year, ".csv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  panel$date <- as.Date(panel$date)
  panel
}

#fill missing training values within each station series
interpolate_training <- function(X) {
  output <- X
  for (column in seq_len(ncol(X))) {
    observed <- which(is.finite(X[, column]))
    if (length(observed) < 2) {
      stop("Insufficient training observations for ", colnames(X)[column])
    }
    output[, column] <- approx(
      x = observed,
      y = X[observed, column],
      xout = seq_len(nrow(X)),
      method = "linear",
      rule = 2
    )$y
  }
  output
}

#missing history is filled using only past values for 2025 (test data)
causal_fill_test <- function(X_train_filled, X_test_observed) {
  output <- rbind(X_train_filled, X_test_observed)
  first_test <- nrow(X_train_filled) + 1
  for (row in first_test:nrow(output)) {
    missing <- !is.finite(output[row, ])
    output[row, missing] <- output[row - 1, missing]
  }
  output
}

#separate AR(1) for each station
fit_ar1 <- function(X_train) {
  coefficients <- matrix(
    NA_real_,
    nrow = ncol(X_train),
    ncol = 2,
    dimnames = list(colnames(X_train), c("intercept", "lag1"))
  )
  for (site in colnames(X_train)) {
    y <- X_train[2:nrow(X_train), site]
    lag1 <- X_train[1:(nrow(X_train) - 1), site]
    fit <- lm.fit(cbind(intercept = 1, lag1 = lag1), y)
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    coefficients[site, ] <- coef
  }
  coefficients
}

forecast_ar1 <- function(X, coefficients, origins) {
  forecast <- matrix(
    NA_real_,
    nrow = length(origins),
    ncol = ncol(X),
    dimnames = list(NULL, colnames(X))
  )
  for (index in seq_along(origins)) {
    forecast[index, ] <-
      coefficients[, "intercept"] + coefficients[, "lag1"] * X[origins[index], ]
  }
  forecast
}

#regression design for VAR
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

fit_ols_var <- function(X, p) {
  design <- var_design(X, p)
  fit <- lm.fit(design$Z, design$Y)
  coefficients <- fit$coefficients
  if (any(!is.finite(coefficients))) {
    stop("Unrestricted OLS VAR coefficients are not uniquely estimable.")
  }
  list(coefficients = coefficients, p = p)
}

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
  forecast <- matrix(
    NA_real_,
    nrow = length(origins),
    ncol = ncol(X),
    dimnames = list(NULL, colnames(X))
  )
  for (index in seq_along(origins)) {
    origin <- origins[index]
    predictor <- 1
    for (lag in seq_len(fit$p)) {
      predictor <- c(predictor, X[origin - lag + 1, ])
    }
    forecast[index, ] <- as.numeric(predictor %*% fit$coefficients)
  }
  forecast
}

#station labels for the hierarchical 
make_station_index <- function(X, p) {
  factor(
    rep(colnames(X), times = nrow(X) - p),
    levels = colnames(X)
  )
}

#hierarchical Bayesian GNAR fitted with station-specific intercepts
fit_hierarchical_station_intercept_chain <- function(
    y,
    R_common,
    station,
    site_ids,
    n_iter = 3500,
    burn = 1000,
    seed = 1
) {
  set.seed(seed)
  Z <- stats::model.matrix(~station - 1)
  colnames(Z) <- paste0("station_", site_ids)
  R <- cbind(R_common, Z)
  k_common <- ncol(R_common)
  d <- length(site_ids)
  n <- length(y)

  fit <- lm.fit(R, y)
  theta <- fit$coefficients
  theta[is.na(theta)] <- 0
  sigma2 <- var(y - as.vector(R %*% theta))
  mu0 <- mean(theta[(k_common + 1):(k_common + d)])
  tau2_mu <- var(theta[(k_common + 1):(k_common + d)])
  if (!is.finite(tau2_mu) || tau2_mu <= 0) tau2_mu <- 1

  draws <- matrix(
    NA_real_,
    n_iter - burn,
    ncol(R),
    dimnames = list(NULL, colnames(R))
  )
  sigma2_draws <- tau2_draws <- mu0_draws <- numeric(n_iter - burn)
  XtX <- crossprod(R)
  Xty <- crossprod(R, y)
  a0 <- b0 <- 0.01
  mu0_prior_var <- 100^2

  for (iter in seq_len(n_iter)) {
    prior_precision <- c(rep(1, k_common), rep(1 / tau2_mu, d))
    prior_mean <- c(rep(0, k_common), rep(mu0, d))
    precision <- XtX / sigma2 + diag(prior_precision, ncol(R))
    rhs <- Xty / sigma2 + prior_precision * prior_mean
    posterior_mean <- solve(precision, rhs)
    theta <- rmvnorm_precision(as.vector(posterior_mean), precision)
    names(theta) <- colnames(R)

    station_intercepts <- theta[(k_common + 1):(k_common + d)]
    mu0_var <- 1 / (d / tau2_mu + 1 / mu0_prior_var)
    mu0_mean <- mu0_var * sum(station_intercepts) / tau2_mu
    mu0 <- rnorm(1, mu0_mean, sqrt(mu0_var))
    tau2_mu <- rinvgamma(
      1,
      shape = a0 + d / 2,
      rate = b0 + sum((station_intercepts - mu0)^2) / 2
    )

    residuals <- y - as.vector(R %*% theta)
    sigma2 <- rinvgamma(
      1,
      shape = a0 + n / 2,
      rate = b0 + sum(residuals^2) / 2
    )

    if (iter > burn) {
      kept <- iter - burn
      draws[kept, ] <- theta
      sigma2_draws[kept] <- sigma2
      tau2_draws[kept] <- tau2_mu
      mu0_draws[kept] <- mu0
    }
  }
  list(theta = draws, sigma2 = sigma2_draws, tau2_mu = tau2_draws, mu0 = mu0_draws)
}

#posterior predictive forecasts for the hierarchical GNAR
posterior_forecast_hierarchical <- function(
    X,
    A,
    draws,
    p,
    stages,
    origins,
    observed = NULL,
    seed = 1
) {
  set.seed(seed)
  n_draws <- nrow(draws$theta)
  common_names <- grep(
    "^station_",
    colnames(draws$theta),
    invert = TRUE,
    value = TRUE
  )
  station_names <- paste0("station_", colnames(X))
  forecast_mean <- lower <- upper <- matrix(
    NA_real_,
    nrow = length(origins),
    ncol = ncol(X),
    dimnames = list(NULL, colnames(X))
  )
  crps <- interval_score_95 <- pit <- matrix(
    NA_real_,
    nrow = length(origins),
    ncol = ncol(X),
    dimnames = list(NULL, colnames(X))
  )

  for (index in seq_along(origins)) {
    origin <- origins[index]
    predictors <- matrix(NA_real_, n_draws, ncol(X))
    for (draw in seq_len(n_draws)) {
      common_theta <- draws$theta[draw, common_names]
      predictors[draw, ] <- gnar_predictor(
        X[seq_len(origin), , drop = FALSE],
        A,
        common_theta,
        p,
        stages
      ) + draws$theta[draw, station_names]
    }
    predictive <- predictors +
      matrix(
        rnorm(length(predictors), sd = rep(sqrt(draws$sigma2), ncol(X))),
        nrow = n_draws
      )
    forecast_mean[index, ] <- colMeans(predictors)
    lower[index, ] <- apply(predictive, 2, quantile, 0.025)
    upper[index, ] <- apply(predictive, 2, quantile, 0.975)
    if (!is.null(observed)) {
      scores <- predictive_draw_scores(predictive, observed[index, ])
      crps[index, ] <- scores$crps
      interval_score_95[index, ] <- scores$interval_score_95
      pit[index, ] <- scores$pit
    }
  }
  list(
    mean = forecast_mean,
    lower = lower,
    upper = upper,
    crps = crps,
    interval_score_95 = interval_score_95,
    pit = pit
  )
}

#point forecast models share the same output columns
summarise_point <- function(model, parameters, forecast) {
  metrics <- forecast_metrics(X_test_observed, forecast)
  data.frame(
    model = model,
    parameters = parameters,
    rmse = metrics$rmse,
    mae = metrics$mae,
    coverage_95 = NA_real_,
    interval_width_95 = NA_real_,
    max_rhat = NA_real_,
    min_bulk_ess = NA_real_,
    min_tail_ess = NA_real_,
    stringsAsFactors = FALSE
  )
}

#intervals and MCMC diagnostics
summarise_probabilistic <- function(model, parameters, forecast, diagnostics) {
  metrics <- forecast_metrics(
    X_test_observed,
    forecast$mean,
    forecast$lower,
    forecast$upper
  )
  data.frame(
    model = model,
    parameters = parameters,
    rmse = metrics$rmse,
    mae = metrics$mae,
    coverage_95 = metrics$coverage_95,
    interval_width_95 = metrics$interval_width_95,
    max_rhat = max(diagnostics$rhat, na.rm = TRUE),
    min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#train on 2022-2024 and test on 2025 data
panels <- lapply(2022:2025, read_panel)
names(panels) <- as.character(2022:2025)
common_sites <- Reduce(intersect, lapply(panels, function(x) setdiff(names(x), "date")))
train_panel <- do.call(
  rbind,
  lapply(panels[c("2022", "2023", "2024")], function(x) x[, c("date", common_sites)])
)
test_panel <- panels[["2025"]][, c("date", common_sites)]

X_train_observed <- as.matrix(train_panel[, common_sites, drop = FALSE])
X_test_observed <- as.matrix(test_panel[, common_sites, drop = FALSE])
storage.mode(X_train_observed) <- "numeric"
storage.mode(X_test_observed) <- "numeric"
X_train <- interpolate_training(X_train_observed)
X_all <- causal_fill_test(X_train, X_test_observed)
train_n <- nrow(X_train)
test_origins <- train_n:(nrow(X_all) - 1)

#keep the same station network and GNAR order as the 2024 analysis
A <- read_matrix_csv(file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv"))
A <- row_normalise(A[common_sites, common_sites, drop = FALSE])
p <- 2
stages <- c(1, 1)

results <- list()

#previous day forecast
persistence <- X_all[test_origins, , drop = FALSE]
results[["Persistence"]] <- summarise_point("Persistence", 0, persistence)

#AR(1) 
ar1_coefficients <- fit_ar1(X_train)
ar1_forecast <- forecast_ar1(X_all, ar1_coefficients, test_origins)
results[["AR(1)"]] <- summarise_point("AR(1)", 2 * ncol(X_train), ar1_forecast)

#unrestricted and ridge VAR 
ols_var_fit <- fit_ols_var(X_train, p = 1)
ols_var_forecast <- forecast_var(X_all, ols_var_fit, test_origins)
results[["Unrestricted OLS VAR"]] <- summarise_point(
  "Unrestricted OLS VAR",
  ncol(X_train) * (1 + ncol(X_train)),
  ols_var_forecast
)

ridge_var_fit <- fit_ridge_var(X_train, p = 2, lambda = 1000)
ridge_var_forecast <- forecast_var(X_all, ridge_var_fit, test_origins)
results[["Ridge VAR"]] <- summarise_point(
  "Ridge VAR",
  ncol(X_train) * (1 + ncol(X_train) * 2),
  ridge_var_forecast
)

#OLS and Bayesian Gaussian GNAR 
gnar_design <- build_gnar_design(X_train, A, p, stages)
ols_gnar_fit <- lm.fit(gnar_design$R, gnar_design$y)
ols_gnar_theta <- ols_gnar_fit$coefficients
ols_gnar_theta[is.na(ols_gnar_theta)] <- 0
names(ols_gnar_theta) <- colnames(gnar_design$R)
ols_gnar_forecast <- forecast_point_path(X_all, A, ols_gnar_theta, p, stages, test_origins)
results[["OLS GNAR"]] <- summarise_point(
  "OLS GNAR",
  ncol(gnar_design$R),
  ols_gnar_forecast
)

gaussian_chains <- lapply(65001:65004, function(seed) {
  fit_gaussian_chain(
    gnar_design$y,
    gnar_design$R,
    coefficient_tau2 = 1,
    intercept_tau2 = 100^2,
    n_iter = 3500,
    burn = 1000,
    seed = seed
  )
})
gaussian_diagnostics <- chain_diagnostics(gaussian_chains)
gaussian_combined <- combine_chains(gaussian_chains)
set.seed(65010)
gaussian_keep <- sample(seq_len(nrow(gaussian_combined$theta)), 2500)
gaussian_forecast <- posterior_forecast(
  X_all,
  A,
  list(
    theta = gaussian_combined$theta[gaussian_keep, , drop = FALSE],
    sigma2 = gaussian_combined$sigma2[gaussian_keep]
  ),
  p,
  stages,
  test_origins,
  observed = X_test_observed,
  seed = 65011
)
results[["Bayesian Gaussian GNAR"]] <- summarise_probabilistic(
  "Bayesian Gaussian GNAR",
  ncol(gnar_design$R),
  gaussian_forecast,
  gaussian_diagnostics
)

#hierarchical Bayesian GNAR 
hier_design <- build_gnar_design(
  X_train,
  A,
  p,
  stages,
  include_intercept = FALSE
)
station <- make_station_index(X_train, p)
hier_chains <- lapply(65101:65104, function(seed) {
  fit_hierarchical_station_intercept_chain(
    hier_design$y,
    hier_design$R,
    station,
    common_sites,
    n_iter = 3500,
    burn = 1000,
    seed = seed
  )
})
hier_combined <- list(
  theta = do.call(rbind, lapply(hier_chains, `[[`, "theta")),
  sigma2 = unlist(lapply(hier_chains, `[[`, "sigma2"), use.names = FALSE),
  tau2_mu = unlist(lapply(hier_chains, `[[`, "tau2_mu"), use.names = FALSE),
  mu0 = unlist(lapply(hier_chains, `[[`, "mu0"), use.names = FALSE)
)
if (requireNamespace("posterior", quietly = TRUE)) {
  iterations <- nrow(hier_chains[[1]]$theta)
  parameter_names <- c(colnames(hier_chains[[1]]$theta), "sigma2", "tau2_mu", "mu0")
  draw_array <- array(
    NA_real_,
    c(iterations, length(hier_chains), length(parameter_names)),
    dimnames = list(NULL, paste0("chain", seq_along(hier_chains)), parameter_names)
  )
  for (chain in seq_along(hier_chains)) {
    draw_array[, chain, ] <- cbind(
      hier_chains[[chain]]$theta,
      sigma2 = hier_chains[[chain]]$sigma2,
      tau2_mu = hier_chains[[chain]]$tau2_mu,
      mu0 = hier_chains[[chain]]$mu0
    )
  }
  hier_diagnostics <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(draw_array),
    "rhat",
    "ess_bulk",
    "ess_tail"
  ))
} else {
  hier_diagnostics <- data.frame(rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_)
}
set.seed(65110)
hier_keep <- sample(seq_len(nrow(hier_combined$theta)), 2500)
hier_forecast <- posterior_forecast_hierarchical(
  X_all,
  A,
  list(
    theta = hier_combined$theta[hier_keep, , drop = FALSE],
    sigma2 = hier_combined$sigma2[hier_keep]
  ),
  p,
  stages,
  test_origins,
  observed = X_test_observed,
  seed = 65111
)
results[["Hierarchical Bayesian GNAR"]] <- summarise_probabilistic(
  "Hierarchical Bayesian GNAR",
  ncol(hier_design$R) + length(common_sites),
  hier_forecast,
  hier_diagnostics
)

summary_table <- do.call(rbind, results)
summary_table$train_years <- "2022-2024"
summary_table$test_year <- 2025
summary_table$stations <- length(common_sites)
summary_table$observed_test_cells <- sum(is.finite(X_test_observed))
summary_table$training_missing_method <- "Within-training linear interpolation"
summary_table$test_history_missing_method <- "Causal last observation carried forward"

model_order <- c(
  "Persistence",
  "AR(1)",
  "Unrestricted OLS VAR",
  "Ridge VAR",
  "OLS GNAR",
  "Bayesian Gaussian GNAR",
  "Hierarchical Bayesian GNAR"
)
summary_table <- summary_table[match(model_order, summary_table$model), ]

#main table used in the report
write.csv(
  summary_table,
  file.path(table_dir, "external_year_pm25_model_comparison_2025.csv"),
  row.names = FALSE
)

#extra
write.csv(
  data.frame(
    date = rep(test_panel$date, each = length(common_sites)),
    site_id = rep(common_sites, times = nrow(test_panel)),
    observed = as.vector(t(X_test_observed)),
    persistence = as.vector(t(persistence)),
    ar1 = as.vector(t(ar1_forecast)),
    unrestricted_ols_var = as.vector(t(ols_var_forecast)),
    ridge_var = as.vector(t(ridge_var_forecast)),
    ols_gnar = as.vector(t(ols_gnar_forecast)),
    bayesian_gaussian_gnar = as.vector(t(gaussian_forecast$mean)),
    hierarchical_bayesian_gnar = as.vector(t(hier_forecast$mean))
  ),
  file.path(table_dir, "external_year_pm25_model_forecasts_2025.csv"),
  row.names = FALSE
)
