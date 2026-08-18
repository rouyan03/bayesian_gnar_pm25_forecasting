project_root <- normalizePath(getwd(), mustWork = TRUE)
analysis_dir <- file.path(project_root, "data", "processed", "analysis")
table_dir <- file.path(project_root, "outputs", "tables")
figure_dir <- file.path(project_root, "outputs", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

sites <- read.csv(
  file.path(table_dir, "main_2024_site_summary.csv"),
  stringsAsFactors = FALSE
)
panel <- read.csv(
  file.path(analysis_dir, "main_2024_pm25_daily_panel.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
panel$date <- as.Date(panel$date)

site_ids <- setdiff(names(panel), "date")
station_mean <- data.frame(
  site_id = site_ids,
  mean_pm25 = colMeans(panel[, site_ids, drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
)
sites <- merge(sites, station_mean, by = "site_id", all.x = TRUE, sort = FALSE)

uk_outline <- ggplot2::map_data("world", region = c("UK", "Ireland"))

station_map <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = uk_outline,
    ggplot2::aes(long, lat, group = group),
    fill = "#EEF1F2",
    colour = "#AAB4B9",
    linewidth = 0.3
  ) +
  ggplot2::geom_point(
    data = sites,
    ggplot2::aes(longitude, latitude, fill = mean_pm25),
    shape = 21,
    size = 4.1,
    colour = "white",
    stroke = 0.65
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  station_map <- station_map +
    ggrepel::geom_label_repel(
      data = sites,
      ggplot2::aes(longitude, latitude, label = site_id),
      size = 2.35,
      linewidth = 0.12,
      label.padding = grid::unit(0.07, "lines"),
      box.padding = 0.18,
      point.padding = 0.14,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 2026,
      fill = "white",
      alpha = 0.92,
      show.legend = FALSE
    )
} else {
  station_map <- station_map +
    ggplot2::geom_text(
      data = sites,
      ggplot2::aes(longitude, latitude, label = site_id),
      nudge_y = 0.12,
      size = 2.4,
      check_overlap = TRUE
    )
}

station_map <- station_map +
  ggplot2::scale_fill_viridis_c(
    option = "C",
    direction = -1,
    name = expression("Mean PM"[2.5])
  ) +
  ggplot2::coord_quickmap(
    xlim = range(sites$longitude) + c(-0.45, 0.45),
    ylim = range(sites$latitude) + c(-0.25, 0.25),
    expand = FALSE,
    clip = "off"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(3, 6, 3, 3)
  )

ggplot2::ggsave(
  file.path(figure_dir, "poster_station_map_mean_only_2024.png"),
  station_map,
  width = 6.2,
  height = 5.7,
  dpi = 320
)

example_sites <- site_ids[seq_len(min(6, length(site_ids)))]
series_data <- data.frame(
  date = rep(panel$date, times = length(example_sites)),
  site_id = rep(example_sites, each = nrow(panel)),
  pm25 = as.vector(as.matrix(panel[, example_sites, drop = FALSE])),
  stringsAsFactors = FALSE
)
series_plot <- ggplot2::ggplot(series_data, ggplot2::aes(date, pm25)) +
  ggplot2::geom_line(colour = "#1D5D7C", linewidth = 0.28) +
  ggplot2::facet_wrap(~site_id, ncol = 2, scales = "free_y") +
  ggplot2::labs(x = NULL, y = expression(PM[2.5])) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(size = 8)
  )

ggplot2::ggsave(
  file.path(figure_dir, "pm25_daily_2024_example_sites_6.png"),
  series_plot,
  width = 7.0,
  height = 6.0,
  dpi = 300
)
