project_root <- normalizePath(getwd(), mustWork = TRUE)
ver_root <- project_root
table_dir <- file.path(ver_root, "outputs", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "R", "gnar_bayes_utils.R"))
#simple ring network (each node has two neighbours)
ring_network <- function(nodes) {
  A <- matrix(0, nodes, nodes)
  for (node in seq_len(nodes)) {
    A[node, ifelse(node == 1, nodes, node - 1)] <- 1
    A[node, ifelse(node == nodes, 1, node + 1)] <- 1
  }
  dimnames(A) <- list(paste0("node", seq_len(nodes)), paste0("node", seq_len(nodes)))
  row_normalise(A)
}
#simulate data from a GNAR(2,[1,1]) model with known parameters
simulate_known_gnar <- function(Tn, A, theta, sigma, burn_in = 200, seed = 6101) {
  set.seed(seed)
  total <- Tn + burn_in
  X <- matrix(0, total, nrow(A), dimnames = list(NULL, rownames(A)))

  for (time in 3:total) {
    mean_value <- theta["intercept"] +
      theta["alpha_lag1"] * X[time - 1, ] +
      theta["beta_lag1_stage1"] * as.numeric(A %*% X[time - 1, ]) +
      theta["alpha_lag2"] * X[time - 2, ] +
      theta["beta_lag2_stage1"] * as.numeric(A %*% X[time - 2, ])
    X[time, ] <- mean_value + rnorm(nrow(A), 0, sigma)
  }

  X[(burn_in + 1):(burn_in + Tn), , drop = FALSE]
}

theta_true <- c(  #true simulation parameters
  intercept = 1.20,
  alpha_lag1 = 0.35,
  beta_lag1_stage1 = 0.30,
  alpha_lag2 = -0.10,
  beta_lag2_stage1 = 0.08
)
sigma_true <- 0.75

nodes <- 30
Tn <- 240
train_end <- 180
p <- 2
stages <- c(1, 1)
#simulate the synthetic network time series
A <- ring_network(nodes)
X <- simulate_known_gnar(Tn, A, theta_true, sigma_true)
#build GNAR design matrix using training data
design <- build_gnar_design(X[seq_len(train_end), , drop = FALSE], A, p, stages)
if (design$status != "ok") stop(design$status)

chain_seeds <- c(6111, 6112, 6113, 6114) #run 4 Bayesian Gaussian GNAR chains
chains <- lapply(
  chain_seeds,
  function(seed) {
    fit_gaussian_chain(
      design$y,
      design$R,
      coefficient_tau2 = 1,
      intercept_tau2 = 100^2,
      n_iter = 5000,
      burn = 1000,
      seed = seed
    )
  }
)

diagnostics <- chain_diagnostics(chains)
diagnostics$simulation <- "known_gnar_2_1_1"
write.csv(
  diagnostics,
  file.path(table_dir, "bayesian_gnar_synthetic_diagnostics.csv"),
  row.names = FALSE
)

combined <- combine_chains(chains) #combine chains and compare estimates with the known values
coefficient_interval <- apply(combined$theta, 2, quantile, probs = c(0.025, 0.975))
posterior_mean <- colMeans(combined$theta)
#fit OLS to the same simulated training data for comparison
ols_fit <- lm.fit(design$R, design$y)
ols_theta <- ols_fit$coefficients
ols_theta[is.na(ols_theta)] <- 0
names(ols_theta) <- colnames(design$R)
#check parameter recovery
parameter_table <- data.frame(
  parameter = names(theta_true),
  true_value = as.numeric(theta_true),
  ols_estimate = as.numeric(ols_theta[names(theta_true)]),
  posterior_mean = as.numeric(posterior_mean[names(theta_true)]),
  lower_95 = as.numeric(coefficient_interval[1, names(theta_true)]),
  upper_95 = as.numeric(coefficient_interval[2, names(theta_true)]),
  covered_95 = as.numeric(
    theta_true >= coefficient_interval[1, names(theta_true)] &
      theta_true <= coefficient_interval[2, names(theta_true)]
  ),
  stringsAsFactors = FALSE
)
parameter_table$absolute_error <- abs(parameter_table$posterior_mean - parameter_table$true_value)

write.csv(
  parameter_table,
  file.path(table_dir, "bayesian_gnar_synthetic_parameter_recovery.csv"),
  row.names = FALSE
)

sigma_interval <- quantile(sqrt(combined$sigma2), probs = c(0.025, 0.975))
origins <- train_end:(nrow(X) - 1)
observed <- X[origins + 1, , drop = FALSE]

set.seed(6199)
forecast_indices <- sample(seq_len(nrow(combined$theta)), min(3000, nrow(combined$theta)))
forecast <- posterior_forecast(
  X,
  A,
  list(
    theta = combined$theta[forecast_indices, , drop = FALSE],
    sigma2 = combined$sigma2[forecast_indices]
  ),
  p,
  stages,
  origins,
  observed = observed,
  seed = 6200
)
metrics <- forecast_metrics(observed, forecast$mean, forecast$lower, forecast$upper)

summary_table <- data.frame(
  simulation = "known GNAR(2,[1,1])",
  nodes = nodes,
  total_days = Tn,
  training_days = train_end,
  test_days = Tn - train_end,
  training_rows = length(design$y),
  true_sigma = sigma_true,
  posterior_sigma_mean = mean(sqrt(combined$sigma2)),
  posterior_sigma_lower_95 = as.numeric(sigma_interval[1]),
  posterior_sigma_upper_95 = as.numeric(sigma_interval[2]),
  parameter_mean_absolute_error = mean(parameter_table$absolute_error),
  parameter_coverage_95 = mean(parameter_table$covered_95),
  forecast_rmse = metrics$rmse,
  forecast_mae = metrics$mae,
  coverage_95 = metrics$coverage_95,
  interval_width_95 = metrics$interval_width_95,
  mean_crps = mean(forecast$crps, na.rm = TRUE),
  max_rhat = max(diagnostics$rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(table_dir, "bayesian_gnar_synthetic_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    time = rep(origins + 1, each = nodes),
    node = rep(colnames(X), times = length(origins)),
    observed = as.vector(t(observed)),
    forecast_mean = as.vector(t(forecast$mean)),
    lower_95 = as.vector(t(forecast$lower)),
    upper_95 = as.vector(t(forecast$upper))
  ),
  file.path(table_dir, "bayesian_gnar_synthetic_forecasts.csv"),
  row.names = FALSE
)
