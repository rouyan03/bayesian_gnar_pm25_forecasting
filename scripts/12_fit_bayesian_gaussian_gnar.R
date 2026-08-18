project_root <- normalizePath(getwd(), mustWork = TRUE)
ver_root <- project_root
analysis_dir <- file.path(ver_root, "data", "processed", "analysis")
network_dir <- file.path(ver_root, "data", "processed", "networks")
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
fit_end <- max(which(panel$date <= as.Date("2024-10-18")))
test_origins <- fit_end:(nrow(X) - 1)
#use the selected equal-neighbour GNAR(2,[1,1]) model
A <- read_matrix_csv(file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv"))
A <- row_normalise(A[site_ids, site_ids, drop = FALSE])
p <- 2
stages <- c(1, 1)
model_order <- "GNAR(2,[1,1])"
#GNAR regression design
design <- build_gnar_design(X[seq_len(fit_end), , drop = FALSE], A, p, stages)
chain_seeds <- c(7101, 7102, 7103, 7104) #run 4 independent MCMC chains
chains <- lapply(
  chain_seeds,
  function(seed) fit_gaussian_chain(
    design$y,
    design$R,
    coefficient_tau2 = 1,
    intercept_tau2 = 100^2,
    n_iter = 5000,
    burn = 1000,
    seed = seed
  )
)
#check MCMC convergence
diagnostics <- chain_diagnostics(chains)
diagnostics$network <- "equal_neighbour"
diagnostics$model_order <- model_order
diagnostics$prior <- "Independent Gaussian: intercept variance 10000; coefficient variance 1"
write.csv(
  diagnostics,
  file.path(table_dir, "bayes_gaussian_gnar_multichain_diagnostics_2024.csv"),
  row.names = FALSE
)
#check stationarity
combined <- combine_chains(chains)
stationarity <- stationarity_summary(combined$theta, A, p, stages, thin = 10)
stationarity$network <- "equal_neighbour"
stationarity$model_order <- model_order
write.csv(
  stationarity,
  file.path(table_dir, "bayes_gaussian_gnar_stationarity_2024.csv"),
  row.names = FALSE
)

set.seed(7199) #use posterior draws to generate test forecasts
forecast_indices <- sample(
  seq_len(nrow(combined$theta)),
  min(4000, nrow(combined$theta))
)
forecast_draws <- list(
  theta = combined$theta[forecast_indices, , drop = FALSE],
  sigma2 = combined$sigma2[forecast_indices]
)
forecast <- posterior_forecast(
  X,
  A,
  forecast_draws,
  p,
  stages,
  test_origins,
  seed = 7200
)
observed <- X[test_origins + 1, , drop = FALSE]
metrics <- forecast_metrics(observed, forecast$mean, forecast$lower, forecast$upper)
#summarise posterior coefficients
coefficient_interval <- apply(combined$theta, 2, quantile, c(0.025, 0.975))
coefficient_table <- data.frame(
  parameter = colnames(combined$theta),
  posterior_mean = colMeans(combined$theta),
  lower_95 = coefficient_interval[1, ],
  upper_95 = coefficient_interval[2, ],
  stringsAsFactors = FALSE
)
write.csv(
  coefficient_table,
  file.path(table_dir, "bayes_gaussian_gnar_multichain_coefficients_2024.csv"),
  row.names = FALSE
)
#summarise & save results
summary_table <- data.frame(
  model = "Bayesian Gaussian GNAR",
  network = "equal_neighbour",
  model_order = model_order,
  prior = "mu~N(0,10000); alpha,beta~N(0,1); sigma2~IG(0.01,0.01)",
  chains = length(chains),
  post_burn_draws = nrow(combined$theta),
  parameters = ncol(design$R),
  rmse = metrics$rmse,
  mae = metrics$mae,
  coverage_95 = metrics$coverage_95,
  interval_width_95 = metrics$interval_width_95,
  max_rhat = max(diagnostics$rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
  max_mcse_mean = max(diagnostics$mcse_mean, na.rm = TRUE),
  stationary_probability = stationarity$stationary_probability,
  stringsAsFactors = FALSE
)
write.csv(
  summary_table,
  file.path(table_dir, "bayes_gaussian_gnar_multichain_2024.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    date = rep(panel$date[test_origins + 1], each = ncol(X)),
    site_id = rep(site_ids, times = length(test_origins)),
    observed = as.vector(t(observed)),
    forecast_mean = as.vector(t(forecast$mean)),
    lower_95 = as.vector(t(forecast$lower)),
    upper_95 = as.vector(t(forecast$upper))
  ),
  file.path(table_dir, "bayes_gaussian_gnar_multichain_forecasts_2024.csv"),
  row.names = FALSE
)
