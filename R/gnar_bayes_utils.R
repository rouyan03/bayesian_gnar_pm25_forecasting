read_matrix_csv <- function(file) {
  A <- as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
  storage.mode(A) <- "numeric"
  A
}

row_normalise <- function(A) {
  row_sums <- rowSums(A)
  W <- matrix(0, nrow(A), ncol(A), dimnames = dimnames(A))
  keep <- is.finite(row_sums) & row_sums > 0
  W[keep, ] <- A[keep, , drop = FALSE] / row_sums[keep]
  W
}

stage_matrix <- function(A, stage) {
  if (stage == 0) return(diag(nrow(A)))
  if (stage == 1) return(row_normalise(A))

  adjacency <- (A > 0) * 1
  diag(adjacency) <- 0
  reach_previous <- adjacency
  reached_earlier <- adjacency

  for (current_stage in 2:stage) {
    reach_current <- (reach_previous %*% adjacency) > 0
    reach_current[reached_earlier > 0] <- FALSE
    diag(reach_current) <- FALSE
    if (current_stage == stage) {
      dimnames(reach_current) <- dimnames(A)
      return(row_normalise(reach_current * 1))
    }
    reach_previous <- reach_current * 1
    reached_earlier <- (reached_earlier + reach_current) > 0
  }

  matrix(0, nrow(A), ncol(A), dimnames = dimnames(A))
}

build_gnar_design <- function(X, A, p, stages_by_lag, include_intercept = TRUE) {
  stopifnot(length(stages_by_lag) == p)
  if (nrow(X) <= p) stop("Time series length must exceed p.")

  d <- ncol(X)
  columns <- list()
  if (include_intercept) {
    columns$intercept <- rep(1, (nrow(X) - p) * d)
  }

  for (lag in seq_len(p)) {
    x_lag <- X[(p + 1 - lag):(nrow(X) - lag), , drop = FALSE]
    columns[[paste0("alpha_lag", lag)]] <- as.vector(t(x_lag))

    if (stages_by_lag[lag] > 0) {
      for (stage in seq_len(stages_by_lag[lag])) {
        W_stage <- stage_matrix(A, stage)
        if (all(W_stage == 0)) {
          return(list(
            y = numeric(0),
            R = matrix(numeric(0), 0, 0),
            status = paste0("not_applicable_empty_stage: beta_lag", lag, "_stage", stage)
          ))
        }
        neighbour_lag <- x_lag %*% t(W_stage)
        columns[[paste0("beta_lag", lag, "_stage", stage)]] <- as.vector(t(neighbour_lag))
      }
    }
  }

  y <- as.vector(t(X[(p + 1):nrow(X), , drop = FALSE]))
  R <- as.matrix(as.data.frame(columns, check.names = FALSE))
  storage.mode(R) <- "numeric"
  keep <- complete.cases(cbind(y, R))
  list(y = y[keep], R = R[keep, , drop = FALSE], status = "ok")
}

gnar_predictor <- function(history, A, theta, p, stages_by_lag) {
  d <- ncol(history)
  prediction <- rep(0, d)
  if ("intercept" %in% names(theta)) prediction <- prediction + theta["intercept"]

  for (lag in seq_len(p)) {
    x_lag <- as.numeric(history[nrow(history) - lag + 1, ])
    prediction <- prediction + theta[paste0("alpha_lag", lag)] * x_lag
    if (stages_by_lag[lag] > 0) {
      for (stage in seq_len(stages_by_lag[lag])) {
        W_stage <- stage_matrix(A, stage)
        prediction <- prediction +
          theta[paste0("beta_lag", lag, "_stage", stage)] *
          as.numeric(W_stage %*% x_lag)
      }
    }
  }

  prediction
}

forecast_point_path <- function(X, A, theta, p, stages_by_lag, origins) {
  forecasts <- matrix(
    NA_real_,
    nrow = length(origins),
    ncol = ncol(X),
    dimnames = list(NULL, colnames(X))
  )
  for (index in seq_along(origins)) {
    origin <- origins[index]
    forecasts[index, ] <- gnar_predictor(
      X[seq_len(origin), , drop = FALSE],
      A,
      theta,
      p,
      stages_by_lag
    )
  }
  forecasts
}

rmse <- function(error) sqrt(mean(error^2, na.rm = TRUE))
mae <- function(error) mean(abs(error), na.rm = TRUE)

forecast_metrics <- function(observed, forecast, lower = NULL, upper = NULL) {
  error <- observed - forecast
  out <- list(rmse = rmse(error), mae = mae(error))
  if (!is.null(lower) && !is.null(upper)) {
    out$coverage_95 <- mean(observed >= lower & observed <= upper, na.rm = TRUE)
    out$interval_width_95 <- mean(upper - lower, na.rm = TRUE)
  } else {
    out$coverage_95 <- NA_real_
    out$interval_width_95 <- NA_real_
  }
  out
}

sample_crps <- function(draws, observation) {
  draws <- sort(draws[is.finite(draws)])
  if (!is.finite(observation) || length(draws) == 0) return(NA_real_)
  m <- length(draws)
  mean(abs(draws - observation)) -
    sum((2 * seq_len(m) - m - 1) * draws) / (m^2)
}

interval_score <- function(observation, lower, upper, alpha = 0.05) {
  if (!is.finite(observation) || !is.finite(lower) || !is.finite(upper)) {
    return(NA_real_)
  }
  upper - lower +
    (2 / alpha) * (lower - observation) * as.numeric(observation < lower) +
    (2 / alpha) * (observation - upper) * as.numeric(observation > upper)
}

predictive_draw_scores <- function(predictive, observed, alpha = 0.05) {
  if (length(observed) != ncol(predictive)) {
    stop("observed must have one value per predictive-draw column.")
  }
  lower <- apply(predictive, 2, quantile, probs = alpha / 2)
  upper <- apply(predictive, 2, quantile, probs = 1 - alpha / 2)
  data.frame(
    crps = vapply(
      seq_along(observed),
      function(column) sample_crps(predictive[, column], observed[column]),
      numeric(1)
    ),
    interval_score_95 = vapply(
      seq_along(observed),
      function(column) interval_score(
        observed[column],
        lower[column],
        upper[column],
        alpha
      ),
      numeric(1)
    ),
    pit = vapply(
      seq_along(observed),
      function(column) {
        if (!is.finite(observed[column])) return(NA_real_)
        mean(predictive[, column] <= observed[column])
      },
      numeric(1)
    )
  )
}

rinvgamma <- function(n, shape, rate) {
  1 / rgamma(n, shape = shape, rate = rate)
}

rmvnorm_precision <- function(mean, precision) {
  R <- chol(precision)
  as.vector(mean + backsolve(R, rnorm(length(mean))))
}

fit_gaussian_chain <- function(
    y,
    R,
    coefficient_tau2 = 1,
    intercept_tau2 = 100^2,
    a0 = 0.01,
    b0 = 0.01,
    n_iter = 5000,
    burn = 1000,
    seed = 1
) {
  set.seed(seed)
  n <- length(y)
  k <- ncol(R)
  theta <- as.vector(lm.fit(R, y)$coefficients)
  theta[is.na(theta)] <- 0
  names(theta) <- colnames(R)
  sigma2 <- var(y - as.vector(R %*% theta))

  prior_variance <- rep(coefficient_tau2, k)
  names(prior_variance) <- colnames(R)
  if ("intercept" %in% names(prior_variance)) {
    prior_variance["intercept"] <- intercept_tau2
  }
  prior_precision <- diag(1 / prior_variance, k)

  theta_draws <- matrix(NA_real_, n_iter - burn, k, dimnames = list(NULL, colnames(R)))
  sigma2_draws <- numeric(n_iter - burn)
  XtX <- crossprod(R)
  Xty <- crossprod(R, y)

  for (iter in seq_len(n_iter)) {
    precision <- XtX / sigma2 + prior_precision
    covariance_mean <- solve(precision, Xty / sigma2)
    theta <- rmvnorm_precision(as.vector(covariance_mean), precision)
    names(theta) <- colnames(R)

    residuals <- y - as.vector(R %*% theta)
    sigma2 <- rinvgamma(
      1,
      shape = a0 + n / 2,
      rate = b0 + sum(residuals^2) / 2
    )

    if (iter > burn) {
      kept <- iter - burn
      theta_draws[kept, ] <- theta
      sigma2_draws[kept] <- sigma2
    }
  }

  list(theta = theta_draws, sigma2 = sigma2_draws)
}

combine_chains <- function(chains) {
  theta <- do.call(rbind, lapply(chains, `[[`, "theta"))
  sigma2 <- unlist(lapply(chains, `[[`, "sigma2"), use.names = FALSE)
  list(theta = theta, sigma2 = sigma2)
}

chain_diagnostics <- function(chains) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("The posterior package is required for multi-chain diagnostics.")
  }
  parameter_names <- c(colnames(chains[[1]]$theta), "sigma2")
  iterations <- nrow(chains[[1]]$theta)
  array_draws <- array(
    NA_real_,
    dim = c(iterations, length(chains), length(parameter_names)),
    dimnames = list(NULL, paste0("chain", seq_along(chains)), parameter_names)
  )
  for (chain in seq_along(chains)) {
    array_draws[, chain, ] <- cbind(chains[[chain]]$theta, sigma2 = chains[[chain]]$sigma2)
  }
  draws <- posterior::as_draws_array(array_draws)
  diagnostics <- posterior::summarise_draws(
    draws,
    "mean",
    "sd",
    ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
    "rhat",
    "ess_bulk",
    "ess_tail",
    "mcse_mean"
  )
  as.data.frame(diagnostics)
}

posterior_forecast <- function(
    X,
    A,
    draws,
    p,
    stages_by_lag,
    origins,
    observed = NULL,
    seed = 1
) {
  set.seed(seed)
  theta_draws <- draws$theta
  sigma2_draws <- draws$sigma2
  n_draws <- nrow(theta_draws)
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
      predictors[draw, ] <- gnar_predictor(
        X[seq_len(origin), , drop = FALSE],
        A,
        theta_draws[draw, ],
        p,
        stages_by_lag
      )
    }
    predictive <- predictors +
      matrix(rnorm(length(predictors), sd = rep(sqrt(sigma2_draws), ncol(X))), nrow = n_draws)
    forecast_mean[index, ] <- colMeans(predictors)
    lower[index, ] <- apply(predictive, 2, quantile, probs = 0.025)
    upper[index, ] <- apply(predictive, 2, quantile, probs = 0.975)
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

companion_spectral_radius <- function(theta, A, p, stages_by_lag) {
  d <- nrow(A)
  lag_matrices <- vector("list", p)
  for (lag in seq_len(p)) {
    lag_matrix <- theta[paste0("alpha_lag", lag)] * diag(d)
    if (stages_by_lag[lag] > 0) {
      for (stage in seq_len(stages_by_lag[lag])) {
        lag_matrix <- lag_matrix +
          theta[paste0("beta_lag", lag, "_stage", stage)] * stage_matrix(A, stage)
      }
    }
    lag_matrices[[lag]] <- lag_matrix
  }

  if (p == 1) return(max(Mod(eigen(lag_matrices[[1]], only.values = TRUE)$values)))

  companion <- matrix(0, d * p, d * p)
  companion[seq_len(d), ] <- do.call(cbind, lag_matrices)
  companion[(d + 1):(d * p), seq_len(d * (p - 1))] <- diag(d * (p - 1))
  max(Mod(eigen(companion, only.values = TRUE)$values))
}

stationarity_summary <- function(theta_draws, A, p, stages_by_lag, thin = 10) {
  indices <- seq(1, nrow(theta_draws), by = thin)
  radii <- vapply(
    indices,
    function(index) companion_spectral_radius(theta_draws[index, ], A, p, stages_by_lag),
    numeric(1)
  )
  data.frame(
    draws_checked = length(radii),
    stationary_probability = mean(radii < 1),
    spectral_radius_mean = mean(radii),
    spectral_radius_q025 = quantile(radii, 0.025),
    spectral_radius_q975 = quantile(radii, 0.975),
    spectral_radius_max = max(radii)
  )
}

analysis_split_2024 <- function(dates) {
  selection_train_end <- max(which(dates <= as.Date("2024-08-31")))
  fit_end <- max(which(dates <= as.Date("2024-10-18")))
  list(
    selection_train_end = selection_train_end,
    fit_end = fit_end,
    validation_origins = selection_train_end:(fit_end - 1),
    test_origins = fit_end:(length(dates) - 1)
  )
}

canonical_orders <- function() {
  list(
    list(name = "GNAR(1,[0])", p = 1, stages = c(0)),
    list(name = "GNAR(1,[1])", p = 1, stages = c(1)),
    list(name = "GNAR(1,[2])", p = 1, stages = c(2)),
    list(name = "GNAR(2,[0,0])", p = 2, stages = c(0, 0)),
    list(name = "GNAR(2,[1,0])", p = 2, stages = c(1, 0)),
    list(name = "GNAR(2,[1,1])", p = 2, stages = c(1, 1)),
    list(name = "GNAR(2,[2,0])", p = 2, stages = c(2, 0)),
    list(name = "GNAR(2,[2,1])", p = 2, stages = c(2, 1)),
    list(name = "GNAR(2,[2,2])", p = 2, stages = c(2, 2))
  )
}

load_canonical_networks <- function(network_dir, site_ids) {
  equal <- read_matrix_csv(
    file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv")
  )[site_ids, site_ids, drop = FALSE]

  distance <- read_matrix_csv(
    file.path(network_dir, "aurn_station_distance_km_2024_nodes.csv")
  )[site_ids, site_ids, drop = FALSE]
  inverse_raw <- 1 / distance
  diag(inverse_raw) <- 0
  inverse_raw[!is.finite(inverse_raw)] <- 0

  knn <- read_matrix_csv(
    file.path(network_dir, "aurn_network_knn4_2024_nodes.csv")
  )[site_ids, site_ids, drop = FALSE]
  knn_union <- ((knn > 0) | (t(knn) > 0)) * 1
  diag(knn_union) <- 0

  threshold <- read_matrix_csv(
    file.path(network_dir, "aurn_network_threshold150km_2024_nodes.csv")
  )[site_ids, site_ids, drop = FALSE]
  threshold <- (threshold > 0) * 1
  diag(threshold) <- 0

  list(
    equal_neighbour = row_normalise(equal),
    inverse_distance = row_normalise(inverse_raw),
    knn4_undirected_union = row_normalise(knn_union),
    threshold150km = row_normalise(threshold)
  )
}

fit_ols_gnar <- function(X, A, p, stages_by_lag) {
  design <- build_gnar_design(X, A, p, stages_by_lag)
  if (design$status != "ok") return(list(status = design$status))
  fit <- lm.fit(design$R, design$y)
  theta <- fit$coefficients
  theta[is.na(theta)] <- 0
  names(theta) <- colnames(design$R)
  residuals <- design$y - as.vector(design$R %*% theta)
  list(
    status = "ok",
    theta = theta,
    residuals = residuals,
    fitted = as.vector(design$R %*% theta),
    design = design
  )
}

fit_student_t_chain <- function(
    y,
    R,
    nu = 8,
    coefficient_tau2 = 1,
    intercept_tau2 = 100^2,
    a0 = 0.01,
    b0 = 0.01,
    n_iter = 5000,
    burn = 1000,
    seed = 1
) {
  if (nu <= 2) stop("nu must exceed 2 for finite innovation variance.")
  set.seed(seed)
  n <- length(y)
  k <- ncol(R)
  theta <- as.vector(lm.fit(R, y)$coefficients)
  theta[is.na(theta)] <- 0
  names(theta) <- colnames(R)
  sigma2 <- var(y - as.vector(R %*% theta))
  omega <- rep(1, n)

  prior_variance <- rep(coefficient_tau2, k)
  names(prior_variance) <- colnames(R)
  if ("intercept" %in% names(prior_variance)) {
    prior_variance["intercept"] <- intercept_tau2
  }
  prior_precision <- diag(1 / prior_variance, k)

  theta_draws <- matrix(
    NA_real_,
    n_iter - burn,
    k,
    dimnames = list(NULL, colnames(R))
  )
  sigma2_draws <- numeric(n_iter - burn)

  for (iter in seq_len(n_iter)) {
    weighted_R <- R * sqrt(omega)
    weighted_y <- y * sqrt(omega)
    precision <- crossprod(weighted_R) / sigma2 + prior_precision
    posterior_mean <- solve(precision, crossprod(weighted_R, weighted_y) / sigma2)
    theta <- rmvnorm_precision(as.vector(posterior_mean), precision)
    names(theta) <- colnames(R)

    residuals <- y - as.vector(R %*% theta)
    sigma2 <- rinvgamma(
      1,
      shape = a0 + n / 2,
      rate = b0 + sum(omega * residuals^2) / 2
    )
    omega <- rgamma(
      n,
      shape = (nu + 1) / 2,
      rate = (nu + residuals^2 / sigma2) / 2
    )

    if (iter > burn) {
      kept <- iter - burn
      theta_draws[kept, ] <- theta
      sigma2_draws[kept] <- sigma2
    }
  }

  list(theta = theta_draws, sigma2 = sigma2_draws, nu = nu)
}

posterior_forecast_student_t <- function(
    X,
    A,
    draws,
    p,
    stages_by_lag,
    origins,
    observed = NULL,
    seed = 1
) {
  set.seed(seed)
  n_draws <- nrow(draws$theta)
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
      predictors[draw, ] <- gnar_predictor(
        X[seq_len(origin), , drop = FALSE],
        A,
        draws$theta[draw, ],
        p,
        stages_by_lag
      )
    }
    innovation <- matrix(
      rt(length(predictors), df = draws$nu) *
        rep(sqrt(draws$sigma2), ncol(X)),
      nrow = n_draws
    )
    predictive <- predictors + innovation
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
