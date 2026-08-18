project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(project_root, "data", "processed", "analysis")
network_dir <- file.path(project_root, "data", "processed", "networks")
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "outputs", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("GNAR", quietly = TRUE)) {
  stop("Package GNAR is required for corbit_plot().")
}
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
if (any(!is.finite(X))) stop("Corbit diagnostics require the complete 2024 panel.")

split <- analysis_split_2024(panel$date)
X_selection <- X[seq_len(split$selection_train_end), , drop = FALSE]
vts <- ts(X_selection, frequency = 1)
colnames(vts) <- site_ids

max_lag <- 12
package_max_stage <- 1
diagnostic_max_stage <- 3

load_plot_network <- function(network_name) {
  if (network_name == "equal_neighbour") {
    A <- read_matrix_csv(
      file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv")
    )
  } else if (network_name == "knn4") {
    A <- read_matrix_csv(
      file.path(network_dir, "aurn_network_knn4_2024_nodes.csv")
    )
    A <- ((A > 0) | (t(A) > 0)) * 1
  } else if (network_name == "threshold150km") {
    A <- read_matrix_csv(
      file.path(network_dir, "aurn_network_threshold150km_2024_nodes.csv")
    )
    A <- (A > 0) * 1
  } else {
    stop("Unknown network: ", network_name)
  }
  A <- A[site_ids, site_ids, drop = FALSE]
  diag(A) <- 0
  A[A < 0 | !is.finite(A)] <- 0
  A
}

#create PNACF corbit plot
make_corbit_plot <- function(A, partial, rectangular_plot, file_stem) {
  gnar_net <- GNAR::matrixtoGNAR(A)
  weight_matrix <- row_normalise(A)
  output_file <- file.path(figure_dir, paste0(file_stem, ".png"))
  png(
    output_file,
    width = ifelse(rectangular_plot == "no", 2400, 2700),
    height = ifelse(rectangular_plot == "no", 2400, 1650),
    res = 300
  )
  on.exit(dev.off(), add = TRUE)
  values <- GNAR::corbit_plot(
    vts = vts,
    net = gnar_net,
    max_lag = max_lag,
    max_stage = package_max_stage,
    weight_matrix = weight_matrix,
    partial = partial,
    rectangular_plot = rectangular_plot,
    viridis_color_option = "magma"
  )
  invisible(values)
}

pnacf_corbit <- make_corbit_plot(
  A = load_plot_network("equal_neighbour"),
  partial = "yes",
  rectangular_plot = "no",
  file_stem = "corbit_pnacf_equal_neighbour_2024_selection"
)

#calculate correlations by lag and network stage.
network_stage_correlation <- function(X_fit, A, lag, stage, partial = FALSE) {
  if (stage == 0) {
    transformed <- X_fit
  } else {
    W_stage <- stage_matrix(A, stage)
    if (all(W_stage == 0)) return(NA_real_)
    transformed <- X_fit %*% t(W_stage)
  }
  target <- as.vector(t(X_fit[(lag + 1):nrow(X_fit), , drop = FALSE]))
  source <- as.vector(t(transformed[seq_len(nrow(X_fit) - lag), , drop = FALSE]))
  keep <- is.finite(target) & is.finite(source)
  if (!partial || lag == 1) {
    return(cor(target[keep], source[keep]))
  }

  controls <- NULL
  for (control_lag in seq_len(lag - 1)) {
    control_values <- as.vector(t(
      X_fit[(lag + 1 - control_lag):(nrow(X_fit) - control_lag), , drop = FALSE]
    ))
    controls <- cbind(controls, control_values)
  }
  keep <- keep & complete.cases(controls)
  if (sum(keep) <= ncol(controls) + 5) return(NA_real_)
  target_residual <- residuals(lm.fit(cbind(1, controls[keep, , drop = FALSE]), target[keep]))
  source_residual <- residuals(lm.fit(cbind(1, controls[keep, , drop = FALSE]), source[keep]))
  cor(target_residual, source_residual)
}

diagnostic_rows <- list()
counter <- 1
for (network_name in c("equal_neighbour", "knn4", "threshold150km")) {
  A <- row_normalise(load_plot_network(network_name))
  for (lag in seq_len(max_lag)) {
    for (stage in 0:diagnostic_max_stage) {
      diagnostic_rows[[counter]] <- data.frame(
        network = network_name,
        lag = lag,
        stage = stage,
        nacf = network_stage_correlation(
          X_selection,
          A,
          lag,
          stage,
          partial = FALSE
        ),
        pnacf = network_stage_correlation(
          X_selection,
          A,
          lag,
          stage,
          partial = TRUE
        )
      )
      counter <- counter + 1
    }
  }
}
diagnostics <- do.call(rbind, diagnostic_rows)
diagnostics$abs_nacf <- abs(diagnostics$nacf)
diagnostics$abs_pnacf <- abs(diagnostics$pnacf)
diagnostics$selection_period_start <- as.character(panel$date[1])
diagnostics$selection_period_end <- as.character(panel$date[split$selection_train_end])
write.csv(
  diagnostics,
  file.path(table_dir, "corbit_network_stage_2024_selection_values.csv"),
  row.names = FALSE
)
write.csv(
  subset(diagnostics, network == "equal_neighbour" & stage == 1),
  file.path(table_dir, "corbit_equal_neighbour_2024_selection_values.csv"),
  row.names = FALSE
)

finite_nacf <- diagnostics[is.finite(diagnostics$abs_nacf), ]
finite_pnacf <- diagnostics[is.finite(diagnostics$abs_pnacf), ]
top_nacf <- finite_nacf[order(-finite_nacf$abs_nacf), ][seq_len(12), ]
top_pnacf <- finite_pnacf[order(-finite_pnacf$abs_pnacf), ][seq_len(12), ]
summary_table <- data.frame(
  diagnostic = c(rep("NACF", nrow(top_nacf)), rep("PNACF", nrow(top_pnacf))),
  rbind(
    top_nacf[, c("network", "lag", "stage", "nacf", "pnacf", "abs_nacf", "abs_pnacf")],
    top_pnacf[, c("network", "lag", "stage", "nacf", "pnacf", "abs_nacf", "abs_pnacf")]
  ),
  row.names = NULL
)
write.csv(
  summary_table,
  file.path(table_dir, "corbit_network_stage_2024_selection_top_values.csv"),
  row.names = FALSE
)
