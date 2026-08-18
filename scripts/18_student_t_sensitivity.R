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

split <- analysis_split_2024(panel$date)
observed_validation <- X[split$validation_origins + 1, , drop = FALSE]
observed_test <- X[split$test_origins + 1, , drop = FALSE]
#load equal-neighbour network.
A <- read_matrix_csv(
  file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv")
)
A <- row_normalise(A[site_ids, site_ids, drop = FALSE])

p <- 2
stages <- c(1, 1)
model_order <- "GNAR(2,[1,1])"
nu_grid <- 3:50 #candidate Student-t degrees of freedom

selection_design <- build_gnar_design(
  X[seq_len(split$selection_train_end), , drop = FALSE],
  A,
  p,
  stages
)
#compare candidate nu values on the validation period
validation_rows <- lapply(seq_along(nu_grid), function(index) {
  nu <- nu_grid[index]
  chain <- fit_student_t_chain(
    selection_design$y,
    selection_design$R,
    nu = nu,
    coefficient_tau2 = 1,
    intercept_tau2 = 100^2,
    n_iter = 3500,
    burn = 900,
    seed = 64000 + nu
  )
  theta_mean <- colMeans(chain$theta)
  forecast <- forecast_point_path(
    X,
    A,
    theta_mean,
    p,
    stages,
    split$validation_origins
  )
  data.frame(
    nu = nu,
    validation_rmspe = rmse(observed_validation - forecast),
    validation_mae = mae(observed_validation - forecast),
    stringsAsFactors = FALSE
  )
})
#rank nu values by validation RMSPE
validation <- do.call(rbind, validation_rows)
validation <- validation[order(validation$validation_rmspe), ]
write.csv(
  validation,
  file.path(table_dir, "bayes_student_t_gnar_nu_validation_3_to_50_2024.csv"),
  row.names = FALSE
)
selected_nu <- validation$nu[1]
final_nu <- 3 #to investigate the effect of a heavy-tailed Student-t likelihood 

final_design <- build_gnar_design(
  X[seq_len(split$fit_end), , drop = FALSE],
  A,
  p,
  stages
)
#run four Student-t Bayesian GNAR chains
chains <- lapply(64101:64104, function(seed) {
	  fit_student_t_chain(
	    final_design$y,
	    final_design$R,
	    nu = final_nu,
    coefficient_tau2 = 1,
    intercept_tau2 = 100^2,
    n_iter = 5000,
    burn = 1000,
    seed = seed
  )
})
#check MCMC convergence
diagnostics <- chain_diagnostics(chains)
diagnostics$network <- "equal_neighbour"
diagnostics$model_order <- model_order
diagnostics$validation_selected_nu <- selected_nu
diagnostics$fitted_nu <- final_nu
diagnostics$likelihood <- "Student-t"
write.csv(
  diagnostics,
  file.path(table_dir, "bayes_student_t_gnar_diagnostics_2024.csv"),
  row.names = FALSE
)
#combine chains and check stationarity
combined <- combine_chains(chains)
stationarity <- stationarity_summary(combined$theta, A, p, stages, thin = 10)
stationarity$network <- "equal_neighbour"
stationarity$model_order <- model_order
stationarity$likelihood <- "Student-t"
stationarity$validation_selected_nu <- selected_nu
stationarity$fitted_nu <- final_nu
write.csv(
  stationarity,
  file.path(table_dir, "bayes_student_t_gnar_stationarity_2024.csv"),
  row.names = FALSE
)
#use posterior draws to generate test forecasts
set.seed(64120)
keep <- sample(seq_len(nrow(combined$theta)), min(4000, nrow(combined$theta)))
forecast <- posterior_forecast_student_t(
  X,
  A,
  list(
	    theta = combined$theta[keep, , drop = FALSE],
	    sigma2 = combined$sigma2[keep],
	    nu = final_nu
  ),
  p,
  stages,
  split$test_origins,
  observed = observed_test,
  seed = 64121
)

metrics <- forecast_metrics(
  observed_test,
  forecast$mean,
  forecast$lower,
  forecast$upper
)
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
  file.path(table_dir, "bayes_student_t_gnar_coefficients_2024.csv"),
  row.names = FALSE
)

coverage <- observed_test >= forecast$lower & observed_test <= forecast$upper
forecast_table <- data.frame(
  date = rep(panel$date[split$test_origins + 1], each = ncol(X)),
  site_id = rep(site_ids, times = length(split$test_origins)),
  observed = as.vector(t(observed_test)),
  forecast_mean = as.vector(t(forecast$mean)),
  lower_95 = as.vector(t(forecast$lower)),
  upper_95 = as.vector(t(forecast$upper)),
  covered_95 = as.vector(t(coverage)),
  crps = as.vector(t(forecast$crps)),
  interval_score_95 = as.vector(t(forecast$interval_score_95)),
  pit = as.vector(t(forecast$pit)),
  stringsAsFactors = FALSE
)
write.csv(
  forecast_table,
  file.path(table_dir, "bayes_student_t_gnar_forecasts_2024.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  model = "Student-t Bayesian GNAR",
  likelihood = "Student-t",
  selected_nu = final_nu,
  validation_selected_nu = selected_nu,
  network = "equal_neighbour",
  model_order = model_order,
  parameters = ncol(final_design$R),
  rmspe = metrics$rmse,
  mae = metrics$mae,
  crps = mean(forecast_table$crps, na.rm = TRUE),
  coverage_95 = metrics$coverage_95,
  interval_width_95 = metrics$interval_width_95,
  max_rhat = max(diagnostics$rhat, na.rm = TRUE),
  min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE),
  stationary_probability = stationarity$stationary_probability,
  stringsAsFactors = FALSE
)
write.csv(
  summary_table,
  file.path(table_dir, "bayes_student_t_gnar_2024.csv"),
  row.names = FALSE
)
