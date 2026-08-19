project_root <- normalizePath(getwd(), mustWork = TRUE)
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "outputs", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
#load hierarchical Bayesian GNAR forecasts
forecast_file <- file.path(
  table_dir,
  "hierarchical_station_intercept_bayes_gnar_forecasts_2024.csv"
)

forecasts <- read.csv(forecast_file, stringsAsFactors = FALSE)
forecasts$date <- as.Date(forecasts$date)

station_rmse <- aggregate(  #calculate RMSE for each station
  (forecast_mean - observed)^2 ~ site_id,
  data = forecasts,
  FUN = function(z) mean(z, na.rm = TRUE)^{1 / 2}
)
names(station_rmse)[2] <- "rmse"
station_rmse <- station_rmse[order(station_rmse$rmse), ]
selected_sites <- unique(station_rmse$site_id[pmax(
  1,
  round(seq(1, nrow(station_rmse), length.out = 6))
)])[seq_len(6)]
plot_data <- subset(forecasts, site_id %in% selected_sites)
#plot posterior means, predictive intervals and observations
interval_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(date, forecast_mean)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lower_95,
      ymax = upper_95,
      fill = "95% predictive interval"
    ),
    alpha = 0.55
  ) +
  ggplot2::geom_line(
    ggplot2::aes(colour = "Posterior mean"),
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = observed, shape = "Observed PM2.5"),
    colour = "#9A3A2F",
    size = 0.9,
    alpha = 0.85
  ) +
  ggplot2::facet_wrap(~site_id, ncol = 2, scales = "free_y") +
  ggplot2::scale_fill_manual(
    values = c("95% predictive interval" = "#C7D9E5"),
    name = NULL
  ) +
  ggplot2::scale_colour_manual(
    values = c("Posterior mean" = "#1D5D7C"),
    name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c("Observed PM2.5" = 16),
    name = NULL
  ) +
  ggplot2::labs(x = NULL, y = "PM2.5") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = grid::unit(0.3, "cm"),
    legend.key.width = grid::unit(1.25, "cm"),
    legend.text = ggplot2::element_text(size = 11.5),
    plot.margin = ggplot2::margin(5, 12, 5, 5)
  )

ggplot2::ggsave(
  file.path(
    figure_dir,
    "poster_hierarchical_bayes_gnar_predictive_intervals_6stations_2024.png"
  ),
  interval_plot,
  width = 8.8,
  height = 8.2,
  dpi = 320
)
