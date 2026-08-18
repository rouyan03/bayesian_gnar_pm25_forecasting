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

selection_train_end <- max(which(panel$date <= as.Date("2024-08-31")))
validation_end <- max(which(panel$date <= as.Date("2024-10-18")))
test_start <- validation_end + 1
#load candidate networks
load_networks <- function(site_ids) {
  equal <- read_matrix_csv(file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv"))
  equal <- equal[site_ids, site_ids, drop = FALSE]

  distance <- read_matrix_csv(file.path(network_dir, "aurn_station_distance_km_2024_nodes.csv"))
  distance <- distance[site_ids, site_ids, drop = FALSE]
  inverse_raw <- 1 / distance
  diag(inverse_raw) <- 0
  inverse_raw[!is.finite(inverse_raw)] <- 0

  knn <- read_matrix_csv(file.path(network_dir, "aurn_network_knn4_2024_nodes.csv"))
  knn <- knn[site_ids, site_ids, drop = FALSE]
  knn_union <- ((knn > 0) | (t(knn) > 0)) * 1
  diag(knn_union) <- 0

  threshold <- read_matrix_csv(file.path(network_dir, "aurn_network_threshold150km_2024_nodes.csv"))
  threshold <- (threshold[site_ids, site_ids, drop = FALSE] > 0) * 1
  diag(threshold) <- 0

  list(
    equal_neighbour = row_normalise(equal),
    inverse_distance = row_normalise(inverse_raw),
    knn4_undirected_union = row_normalise(knn_union),
    threshold150km = row_normalise(threshold)
  )
}
#candidate GNAR orders considered during validation
orders <- list(
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

fit_ols_gnar <- function(X_fit, A, order) {
  design <- build_gnar_design(X_fit, A, order$p, order$stages)
  if (design$status != "ok") return(list(status = design$status))
  fit <- lm.fit(design$R, design$y)
  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  names(coefficients) <- colnames(design$R)

  residual_vector <- design$y - as.vector(design$R %*% coefficients)
  residual_matrix <- matrix(
    residual_vector,
    nrow = nrow(X_fit) - order$p,
    ncol = ncol(X_fit),
    byrow = TRUE
  )
  covariance <- crossprod(residual_matrix) / nrow(residual_matrix)
  determinant <- determinant(covariance, logarithm = TRUE)
  bic <- if (determinant$sign > 0) {
    as.numeric(determinant$modulus) +
      length(coefficients) * log(nrow(residual_matrix)) / nrow(residual_matrix)
  } else {
    NA_real_
  }

  list(
    status = "ok",
    coefficients = coefficients,
    parameters = length(coefficients),
    bic = bic
  )
}
#produce one-day-ahead forecasts & calculate forecast errors
evaluate_period <- function(X, A, order, coefficients, origins) {
  forecast <- forecast_point_path(X, A, coefficients, order$p, order$stages, origins)
  observed <- X[origins + 1, , drop = FALSE]
  metrics <- forecast_metrics(observed, forecast)
  list(
    forecast = forecast,
    observed = observed,
    rmse = metrics$rmse,
    mae = metrics$mae,
    daily_rmse = apply((observed - forecast)^2, 1, function(x) sqrt(mean(x)))
  )
}

networks <- load_networks(site_ids)
#fit every network order combination on training period
validation_rows <- list()
counter <- 1

for (network_name in names(networks)) {
  A <- networks[[network_name]]
  for (order in orders) {
    fit <- fit_ols_gnar(X[seq_len(selection_train_end), , drop = FALSE], A, order)
    if (fit$status != "ok") {
      validation_rows[[counter]] <- data.frame(
        network = network_name,
        model_order = order$name,
        p = order$p,
        stages = paste(order$stages, collapse = ","),
        parameters = NA_integer_,
        training_bic = NA_real_,
        validation_rmse = NA_real_,
        validation_mae = NA_real_,
        status = fit$status
      )
    } else {
      validation_origins <- selection_train_end:(validation_end - 1)
      evaluated <- evaluate_period(X, A, order, fit$coefficients, validation_origins)
      validation_rows[[counter]] <- data.frame(
        network = network_name,
        model_order = order$name,
        p = order$p,
        stages = paste(order$stages, collapse = ","),
        parameters = fit$parameters,
        training_bic = fit$bic,
        validation_rmse = evaluated$rmse,
        validation_mae = evaluated$mae,
        status = "ok"
      )
    }
    counter <- counter + 1
  }
}
#rank models by validation RMSE
validation_table <- do.call(rbind, validation_rows)
validation_table <- validation_table[
  order(is.na(validation_table$validation_rmse), validation_table$validation_rmse, validation_table$parameters),
]
write.csv(
  validation_table,
  file.path(table_dir, "ols_gnar_selection_validation_2024_corrected.csv"),
  row.names = FALSE
)
#select the best specification and refit it using training + validation data
selected <- validation_table[validation_table$status == "ok", ][1, ]
selected_order <- orders[[match(selected$model_order, vapply(orders, `[[`, "", "name"))]]
selected_A <- networks[[selected$network]]
final_fit <- fit_ols_gnar(X[seq_len(validation_end), , drop = FALSE], selected_A, selected_order)
test_origins <- validation_end:(nrow(X) - 1) #evaluate the selected model on test set
test <- evaluate_period(X, selected_A, selected_order, final_fit$coefficients, test_origins)

test_summary <- data.frame(
  selection_rule = "minimum validation RMSE; parameters break ties",
  network = selected$network,
  model_order = selected$model_order,
  parameters = final_fit$parameters,
  selection_train_start = as.character(panel$date[1]),
  selection_train_end = as.character(panel$date[selection_train_end]),
  validation_start = as.character(panel$date[selection_train_end + 1]),
  validation_end = as.character(panel$date[validation_end]),
  test_start = as.character(panel$date[test_start]),
  test_end = as.character(panel$date[nrow(panel)]),
  validation_rmse = selected$validation_rmse,
  validation_mae = selected$validation_mae,
  test_rmse = test$rmse,
  test_mae = test$mae,
  coefficient_summary = paste(
    names(final_fit$coefficients),
    round(final_fit$coefficients, 5),
    sep = "=",
    collapse = "; "
  ),
  stringsAsFactors = FALSE
)
write.csv(
  test_summary,
  file.path(table_dir, "ols_gnar_selected_test_2024_corrected.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    date = panel$date[test_origins + 1],
    daily_rmse = test$daily_rmse,
    network = selected$network,
    model_order = selected$model_order
  ),
  file.path(table_dir, "ols_gnar_selected_test_daily_rmse_2024_corrected.csv"),
  row.names = FALSE
)
#plot validation RMSE for all candidate models
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  plot_data <- subset(validation_table, status == "ok")
  plot_data$network_label <- c(
    equal_neighbour = "Equal neighbour",
    inverse_distance = "Inverse distance",
    knn4_undirected_union = "4-nearest neighbour",
    threshold150km = "150 km threshold"
  )[plot_data$network]
  plot_data$network_label <- factor(
    plot_data$network_label,
    levels = c(
      "Equal neighbour",
      "Inverse distance",
      "4-nearest neighbour",
      "150 km threshold"
    )
  )
  plot_data$model_order <- factor(
    plot_data$model_order,
    levels = rev(vapply(orders, `[[`, "", "name"))
  )
  selected_tile <- subset(
    plot_data,
    network == selected$network & model_order == selected$model_order
  )
  p <- ggplot(plot_data, aes(x = network_label, y = model_order, fill = validation_rmse)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.3f", validation_rmse)), size = 3) +
    geom_tile(
      data = selected_tile,
      fill = NA,
      colour = "black",
      linewidth = 0.8
    ) +
    scale_fill_gradient2(
      low = "#C9DFEE",
      mid = "#F7F3E8",
      high = "#CC7569",
      midpoint = 3.80
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Validation RMSE"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(
    file.path(figure_dir, "ols_gnar_selection_validation_heatmap_2024_thesis_blue_red.png"),
    p,
    width = 6.2,
    height = 4.4,
    dpi = 300
  )
}
