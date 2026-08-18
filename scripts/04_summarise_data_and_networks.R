project_root <- normalizePath(getwd(), mustWork = TRUE)

processed_dir <- file.path(project_root, "data", "processed", "defra_aurn")
analysis_dir <- file.path(project_root, "data", "processed", "analysis")
network_dir <- file.path(project_root, "data", "processed", "networks")
table_dir <- file.path(project_root, "outputs", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

#check these scripts available
required <- c(
  file.path(processed_dir, "aurn_pm25_panel_summary_2022_2025.csv"),
  file.path(analysis_dir, "main_2024_panel_summary.csv"),
  file.path(analysis_dir, "main_2024_site_completeness.csv"),
  file.path(analysis_dir, "main_2024_train_test_split.csv"),
  file.path(network_dir, "aurn_network_summary_2024_nodes.csv")
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required files. Run scripts 01, 02, and 03 first:\n", paste(missing, collapse = "\n"))
}
#load data, split and network summaries
panel_22_25 <- read.csv(file.path(processed_dir, "aurn_pm25_panel_summary_2022_2025.csv"))
main_summary <- read.csv(file.path(analysis_dir, "main_2024_panel_summary.csv"), stringsAsFactors = FALSE)
site_completeness <- read.csv(file.path(analysis_dir, "main_2024_site_completeness.csv"), stringsAsFactors = FALSE)
split <- read.csv(file.path(analysis_dir, "main_2024_train_test_split.csv"), stringsAsFactors = FALSE)
network_summary <- read.csv(file.path(network_dir, "aurn_network_summary_2024_nodes.csv"), stringsAsFactors = FALSE)

dataset_table <- data.frame(
  dataset = c(
    "main_2024",
    paste0("expanded_", panel_22_25$year)
  ),
  role = c(
    "main model comparison",
    rep("robustness/external validation input", nrow(panel_22_25))
  ),
  nodes = c(
    as.numeric(main_summary$value[main_summary$quantity == "nodes"]),
    panel_22_25$nodes
  ),
  daily_time_points = c(
    as.numeric(main_summary$value[main_summary$quantity == "daily_time_points"]),
    panel_22_25$daily_rows
  ),
  daily_cells = c(
    as.numeric(main_summary$value[main_summary$quantity == "total_daily_cells"]),
    panel_22_25$daily_cells
  ),
  daily_missing = c(
    as.numeric(main_summary$value[main_summary$quantity == "missing_daily_cells"]),
    panel_22_25$daily_missing
  ),
  daily_completeness = c(
    as.numeric(main_summary$value[main_summary$quantity == "overall_daily_completeness"]),
    1 - panel_22_25$daily_missing / panel_22_25$daily_cells
  )
)

site_table <- site_completeness[order(site_completeness$site_id), c(
  "site_id",
  "site_name",
  "environment_type",
  "latitude",
  "longitude",
  "daily_observed",
  "daily_missing",
  "daily_completeness"
)]

network_table <- network_summary[order(network_summary$network), ]

write.csv(dataset_table, file.path(table_dir, "dataset_summary_2022_2025_and_main_2024.csv"), row.names = FALSE)
write.csv(split, file.path(table_dir, "main_2024_train_test_split.csv"), row.names = FALSE)
write.csv(site_table, file.path(table_dir, "main_2024_site_summary.csv"), row.names = FALSE)
write.csv(network_table, file.path(table_dir, "network_summary_2024_nodes.csv"), row.names = FALSE)
