project_root <- normalizePath(getwd(), mustWork = TRUE)
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "outputs", "figures")

forecast_file <- file.path(
  table_dir,
  "hierarchical_station_intercept_bayes_gnar_forecasts_2024.csv"
)
if (!file.exists(forecast_file)) {
  stop("Missing hierarchical Bayesian GNAR forecasts. Run script 13 first.")
}

forecasts <- read.csv(forecast_file, stringsAsFactors = FALSE)
#group observations into concentration deciles
forecasts$concentration_decile <- cut(
  forecasts$observed,
  breaks = quantile(forecasts$observed, probs = seq(0, 1, 0.1), na.rm = TRUE),
  include.lowest = TRUE,
  labels = 1:10
)
forecasts$concentration_decile <- as.integer(as.character(forecasts$concentration_decile))
#calculate coverage, forecast error and interval width
forecasts$covered_95 <- forecasts$observed >= forecasts$lower_95 &
  forecasts$observed <= forecasts$upper_95
forecasts$forecast_error <- forecasts$forecast_mean - forecasts$observed
forecasts$interval_width <- forecasts$upper_95 - forecasts$lower_95

decile_summary <- aggregate(
  cbind(covered_95, forecast_error, interval_width) ~ concentration_decile,
  data = forecasts,
  FUN = mean
)

write.csv(
  decile_summary,
  file.path(table_dir, "hierarchical_peak_diagnostics_by_decile_2024.csv"),
  row.names = FALSE
)
#plot coverage and mean forecast error by concentration decile
if (requireNamespace("ggplot2", quietly = TRUE)) {
  diagnostic_data <- rbind(
    data.frame(
      concentration_decile = decile_summary$concentration_decile,
      panel = "95% predictive coverage",
      value = decile_summary$covered_95
    ),
    data.frame(
      concentration_decile = decile_summary$concentration_decile,
      panel = "Mean error (forecast - observed)",
      value = decile_summary$forecast_error
    )
  )
  diagnostic_data$panel <- factor(
    diagnostic_data$panel,
    levels = c("95% predictive coverage", "Mean error (forecast - observed)")
  )
  reference_lines <- data.frame(
    panel = factor(
      c("95% predictive coverage", "Mean error (forecast - observed)"),
      levels = levels(diagnostic_data$panel)
    ),
    intercept = c(0.95, 0)
  )

  plot <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(concentration_decile, value)
  ) +
    ggplot2::geom_hline(
      data = reference_lines,
      ggplot2::aes(yintercept = intercept),
      linetype = "dashed",
      colour = "#6B747B",
      linewidth = 0.55
    ) +
    ggplot2::geom_line(colour = "#276F55", linewidth = 1.05) +
    ggplot2::geom_point(colour = "#276F55", size = 2.8) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = 1:10) +
    ggplot2::labs(
      x = "Observed PM2.5 concentration decile",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "#EAF1F4",
        colour = "#C4D2D9"
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = "#1C2C39",
        margin = ggplot2::margin(4, 4, 4, 4)
      ),
      axis.text = ggplot2::element_text(colour = "#333333"),
      panel.spacing = grid::unit(0.55, "lines"),
      plot.margin = ggplot2::margin(5, 8, 5, 5)
    )

  ggplot2::ggsave(
    file.path(figure_dir, "poster_hierarchical_peak_calibration_bias_2024.png"),
    plot,
    width = 8.8,
    height = 6.2,
    dpi = 320
  )
}
