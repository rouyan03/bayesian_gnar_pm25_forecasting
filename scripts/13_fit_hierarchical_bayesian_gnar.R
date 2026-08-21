project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(project_root, "data", "processed", "analysis")
network_dir <- file.path(project_root, "data", "processed", "networks")
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
if (any(!is.finite(X))) stop("The hierarchical model requires a complete 2024 panel.")
#load selected equal-neigh network
A <- read_matrix_csv(file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv"))
A <- row_normalise(A[site_ids, site_ids, drop = FALSE])
#use the selected GNAR(2,[1,1]) 
split <- analysis_split_2024(panel$date)
p <- 2
stages <- c(1, 1)
#match each stacked response to its monitoring station
make_station_index <- function(X_fit, p) {
  factor(
    rep(colnames(X_fit), times = nrow(X_fit) - p),
    levels = colnames(X_fit)
  )
}
#fit one hierarchical Bayesian GNAR chain
fit_hierarchical_chain <- function(y, R_common, station, site_ids, n_iter, burn, seed) {
  set.seed(seed)
  Z <- model.matrix(~station - 1) #add one intercept column for each station
  colnames(Z) <- paste0("station_", site_ids)
  R <- cbind(R_common, Z)
  k_common <- ncol(R_common)
  d <- length(site_ids)
  n <- length(y)

  #use OLS values to initialise the chain
  theta <- lm.fit(R, y)$coefficients
  theta[is.na(theta)] <- 0
  sigma2 <- var(y - as.vector(R %*% theta))
  mu0 <- mean(theta[(k_common + 1):(k_common + d)])
  tau2_mu <- var(theta[(k_common + 1):(k_common + d)])
  if (!is.finite(tau2_mu) || tau2_mu <= 0) tau2_mu <- 1
  
  theta_draws <- matrix( #store post-burn-in draws
    NA_real_,
    n_iter - burn,
    ncol(R),
    dimnames = list(NULL, colnames(R))
  )
  sigma2_draws <- tau2_draws <- mu0_draws <- numeric(n_iter - burn)

  XtX <- crossprod(R)
  Xty <- crossprod(R, y)
  a0 <- b0 <- 0.01
  mu0_prior_variance <- 100^2

  for (iteration in seq_len(n_iter)) {  #gibbs sampler
    prior_precision <- c(rep(1, k_common), rep(1 / tau2_mu, d))
    prior_mean <- c(rep(0, k_common), rep(mu0, d))
    precision <- XtX / sigma2 + diag(prior_precision, ncol(R))
    rhs <- as.vector(Xty) / sigma2 + prior_precision * prior_mean
    theta <- rmvnorm_precision(as.vector(solve(precision, rhs)), precision)
    names(theta) <- colnames(R)
    #update the population mean of the station intercepts
    station_intercepts <- theta[(k_common + 1):(k_common + d)]
    mu0_variance <- 1 / (d / tau2_mu + 1 / mu0_prior_variance)
    mu0_mean <- mu0_variance * sum(station_intercepts) / tau2_mu
    mu0 <- rnorm(1, mu0_mean, sqrt(mu0_variance))
    tau2_mu <- rinvgamma(   #update the between-station intercept variance
      1,
      shape = a0 + d / 2,
      rate = b0 + sum((station_intercepts - mu0)^2) / 2
    )

    residual <- y - as.vector(R %*% theta)  #update the residual variance
    sigma2 <- rinvgamma(
      1,
      shape = a0 + n / 2,
      rate = b0 + sum(residual^2) / 2
    )
    #keep draws after burn-in
    if (iteration > burn) {
      kept <- iteration - burn
      theta_draws[kept, ] <- theta
      sigma2_draws[kept] <- sigma2
      tau2_draws[kept] <- tau2_mu
      mu0_draws[kept] <- mu0
    }
  }

  list(theta = theta_draws, sigma2 = sigma2_draws, tau2_mu = tau2_draws, mu0 = mu0_draws)
}
#generate hierarchical posterior predictive forecasts
forecast_hierarchical <- function(X, A, draws, p, stages, origins, observed, seed) {
  set.seed(seed)
  theta_draws <- draws$theta
  sigma2_draws <- draws$sigma2
  n_draws <- nrow(theta_draws)

  forecast_mean <- lower <- upper <- crps <- interval_score_95 <- pit <- matrix(
    NA_real_,
    nrow = length(origins),
    ncol = ncol(X),
    dimnames = list(NULL, colnames(X))
  )

  for (index in seq_along(origins)) {
    origin <- origins[index]
    predictors <- matrix(NA_real_, n_draws, ncol(X), dimnames = list(NULL, colnames(X)))
    #calculate the conditional mean for each posterior draw
    for (draw in seq_len(n_draws)) {
      theta <- theta_draws[draw, ]
      common_theta <- theta[!grepl("^station_", names(theta))]
      station_intercepts <- theta[paste0("station_", colnames(X))]
      predictors[draw, ] <- station_intercepts + gnar_predictor(
        X[seq_len(origin), , drop = FALSE],
        A,
        common_theta,
        p,
        stages
      )
    }

    predictive <- predictors +
      matrix(rnorm(length(predictors), sd = rep(sqrt(sigma2_draws), ncol(X))), nrow = n_draws)

    forecast_mean[index, ] <- colMeans(predictors)
    lower[index, ] <- apply(predictive, 2, quantile, 0.025)
    upper[index, ] <- apply(predictive, 2, quantile, 0.975)
    scores <- predictive_draw_scores(predictive, observed[index, ])
    crps[index, ] <- scores$crps
    interval_score_95[index, ] <- scores$interval_score_95
    pit[index, ] <- scores$pit
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
#build the common GNAR predictors without a global intercept
design <- build_gnar_design(
  X[seq_len(split$fit_end), , drop = FALSE],
  A,
  p,
  stages,
  include_intercept = FALSE
)
station <- make_station_index(X[seq_len(split$fit_end), , drop = FALSE], p)
#run 4 independent MCMC chains
chains <- lapply(42301:42304, function(seed) {
  fit_hierarchical_chain(
    design$y,
    design$R,
    station,
    site_ids,
    n_iter = 5000,
    burn = 1200,
    seed = seed
  )
})

combined <- list(
  theta = do.call(rbind, lapply(chains, `[[`, "theta")),
  sigma2 = unlist(lapply(chains, `[[`, "sigma2"), use.names = FALSE),
  tau2_mu = unlist(lapply(chains, `[[`, "tau2_mu"), use.names = FALSE),
  mu0 = unlist(lapply(chains, `[[`, "mu0"), use.names = FALSE)
)
#check stationarity
common_names <- grep("^station_", colnames(combined$theta), invert = TRUE, value = TRUE)
stationarity <- stationarity_summary(
  combined$theta[, common_names, drop = FALSE],
  A,
  p,
  stages,
  thin = 10
)

if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Package posterior is required for MCMC diagnostics.")
}

iterations <- nrow(chains[[1]]$theta)
parameter_names <- c(colnames(chains[[1]]$theta), "sigma2", "tau2_mu", "mu0")
draw_array <- array(
  NA_real_,
  c(iterations, length(chains), length(parameter_names)),
  dimnames = list(NULL, paste0("chain", seq_along(chains)), parameter_names)
)
for (chain in seq_along(chains)) {
  draw_array[, chain, ] <- cbind(
    chains[[chain]]$theta,
    sigma2 = chains[[chain]]$sigma2,
    tau2_mu = chains[[chain]]$tau2_mu,
    mu0 = chains[[chain]]$mu0
  )
}
diagnostics <- as.data.frame(posterior::summarise_draws(
  posterior::as_draws_array(draw_array),
  "rhat",
  "ess_bulk",
  "ess_tail"
))
#summarise posterior coefficient distributions
interval <- apply(combined$theta, 2, quantile, c(0.025, 0.975))
coefficient_table <- data.frame(
  parameter = colnames(combined$theta),
  posterior_mean = colMeans(combined$theta),
  lower_95 = interval[1, ],
  upper_95 = interval[2, ],
  stringsAsFactors = FALSE
)

set.seed(42399)
keep <- sample(seq_len(nrow(combined$theta)), min(3500, nrow(combined$theta)))
observed <- X[split$test_origins + 1, , drop = FALSE]
forecast <- forecast_hierarchical(
  X,
  A,
  list(theta = combined$theta[keep, , drop = FALSE], sigma2 = combined$sigma2[keep]),
  p,
  stages,
  split$test_origins,
  observed,
  seed = 42400
)
#calculate point and probabilistic forecast scores
metrics <- forecast_metrics(observed, forecast$mean, forecast$lower, forecast$upper)
summary_table <- data.frame(
  model = "Hierarchical Bayesian GNAR",
  network = "equal_neighbour",
  model_order = "GNAR(2,[1,1])",
  parameters = ncol(combined$theta),
  rmse = metrics$rmse,
  mae = metrics$mae,
  coverage_95 = metrics$coverage_95,
  interval_width_95 = metrics$interval_width_95,
  crps = mean(forecast$crps),
  interval_score_95 = mean(forecast$interval_score_95),
  max_rhat = max(diagnostics$rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
  stationary_probability = stationarity$stationary_probability,
  tau_mu_mean = mean(sqrt(combined$tau2_mu)),
  stringsAsFactors = FALSE
)

forecast_table <- data.frame(
  date = rep(panel$date[split$test_origins + 1], each = length(site_ids)),
  site_id = rep(site_ids, times = length(split$test_origins)),
  observed = as.vector(t(observed)),
  forecast_mean = as.vector(t(forecast$mean)),
  lower_95 = as.vector(t(forecast$lower)),
  upper_95 = as.vector(t(forecast$upper)),
  crps = as.vector(t(forecast$crps)),
  interval_score_95 = as.vector(t(forecast$interval_score_95)),
  pit = as.vector(t(forecast$pit)),
  stringsAsFactors = FALSE
)

theta_hat <- colMeans(combined$theta)
Z <- model.matrix(~station - 1)
colnames(Z) <- paste0("station_", site_ids)
R_hierarchical <- cbind(design$R, Z)
fitted_train <- as.vector(R_hierarchical %*% theta_hat)
residual_train <- design$y - fitted_train
residual_matrix <- matrix(
  residual_train,
  nrow = split$fit_end - p,
  ncol = length(site_ids),
  byrow = TRUE,
  dimnames = list(NULL, site_ids)
)
mean_residual_by_day <- rowMeans(residual_matrix)
acf_values <- as.numeric(stats::acf(mean_residual_by_day, plot = FALSE, lag.max = 14)$acf)[-1]
residual_table <- data.frame(
  date = rep(panel$date[(p + 1):split$fit_end], each = length(site_ids)),
  site_id = rep(site_ids, times = split$fit_end - p),
  observed = design$y,
  fitted = fitted_train,
  residual = residual_train,
  stringsAsFactors = FALSE
)
residual_summary <- data.frame(
  model = "Hierarchical Bayesian GNAR",
  n_residuals = length(residual_train),
  mean_residual = mean(residual_train),
  sd_residual = sd(residual_train),
  q025 = as.numeric(quantile(residual_train, 0.025)),
  q500 = as.numeric(quantile(residual_train, 0.500)),
  q975 = as.numeric(quantile(residual_train, 0.975)),
  mean_abs_residual = mean(abs(residual_train)),
  acf_lag1_daily_mean = acf_values[1],
  acf_lag2_daily_mean = acf_values[2],
  acf_lag7_daily_mean = acf_values[7],
  stringsAsFactors = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  figure_dir <- file.path(project_root, "outputs", "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  qq_data <- data.frame(
    theoretical = qnorm(ppoints(length(residual_train))),
    residual = sort(residual_train)
  )

  residual_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(theoretical, residual)) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = sd(residual_train),
      colour = "#A33C31",
      linetype = "dashed"
    ) +
    ggplot2::geom_point(size = 0.8, alpha = 0.55, colour = "#1D5D7C") +
    ggplot2::labs(x = "Normal quantile", y = "Residual") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "#333333")
    )

  ggplot2::ggsave(
    file.path(figure_dir, "hierarchical_bayes_gnar_residual_diagnostics_2024.png"),
    residual_plot,
    width = 5.8,
    height = 4.6,
    dpi = 320
  )
}

write.csv(summary_table, file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_2024.csv"), row.names = FALSE)
write.csv(coefficient_table, file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_coefficients_2024.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_diagnostics_2024.csv"), row.names = FALSE)
write.csv(forecast_table, file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_forecasts_2024.csv"), row.names = FALSE)
write.csv(residual_table, file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_residuals_2024.csv"), row.names = FALSE)
write.csv(residual_summary, file.path(table_dir, "hierarchical_station_intercept_bayes_gnar_residual_summary_2024.csv"), row.names = FALSE)
