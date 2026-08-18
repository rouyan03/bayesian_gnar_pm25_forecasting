project_root <- normalizePath(getwd(), mustWork = TRUE)

processed_dir <- file.path(project_root, "data", "processed", "defra_aurn")
analysis_dir <- file.path(project_root, "data", "processed", "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

#load 2024 daily panel and selected 30 stations
panel_file <- file.path(processed_dir, "aurn_pm25_daily_panel_2024.csv")
sites_file <- file.path(processed_dir, "aurn_pm25_selected_sites_fixed_2024_nodes.csv")

if (!file.exists(panel_file)) {
  stop("Missing 2024 daily panel. Run final_code/scripts/01_extract_aurn_data.R first.")
}
if (!file.exists(sites_file)) {
  stop("Missing selected-site file. Run final_code/scripts/01_extract_aurn_data.R first.")
}

panel <- read.csv(panel_file, check.names = FALSE, stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
sites <- read.csv(sites_file, stringsAsFactors = FALSE)
site_ids <- sites$site_id

#check all selected stations are present
missing_sites <- setdiff(site_ids, names(panel))
if (length(missing_sites) > 0) {
  stop("These sites are missing: ", paste(missing_sites, collapse = ", "))
}

panel <- panel[, c("date", site_ids), drop = FALSE]
panel <- panel[order(panel$date), , drop = FALSE]

X <- as.matrix(panel[, site_ids, drop = FALSE])
storage.mode(X) <- "numeric"

#define the training, validation and test dataset
split <- data.frame(
  split = c("training", "validation", "test"),
  start_date = c("2024-01-01", "2024-09-01", "2024-10-19"),
  end_date = c("2024-08-31", "2024-10-18", "2024-12-31"),
  rows = c(
    sum(panel$date >= as.Date("2024-01-01") & panel$date <= as.Date("2024-08-31")),
    sum(panel$date >= as.Date("2024-09-01") & panel$date <= as.Date("2024-10-18")),
    sum(panel$date >= as.Date("2024-10-19") & panel$date <= as.Date("2024-12-31"))
  )
)

site_missing <- data.frame(
  site_id = site_ids,
  daily_observed = colSums(!is.na(X)),
  daily_missing = colSums(is.na(X)),
  daily_total = nrow(X),
  daily_completeness = colSums(!is.na(X)) / nrow(X),
  stringsAsFactors = FALSE
)
site_missing <- merge(site_missing, sites, by = "site_id", all.x = TRUE)
#summarise 2024 panel
panel_summary <- data.frame(
  quantity = c(
    "year",
    "nodes",
    "daily_time_points",
    "training_days",
    "validation_days",
    "test_days",
    "first_date",
    "last_date",
    "total_daily_cells",
    "missing_daily_cells",
    "overall_daily_completeness",
    "mean_pm25",
    "sd_pm25",
    "min_pm25",
    "max_pm25"
  ),
  value = c(
    2024,
    length(site_ids),
    nrow(panel),
    split$rows[split$split == "training"],
    split$rows[split$split == "validation"],
    split$rows[split$split == "test"],
    as.character(panel$date[1]),
    as.character(panel$date[nrow(panel)]),
    length(X),
    sum(is.na(X)),
    mean(!is.na(X)),
    mean(X, na.rm = TRUE),
    sd(as.vector(X), na.rm = TRUE),
    min(X, na.rm = TRUE),
    max(X, na.rm = TRUE)
  )
)
#save the analysis data & summary files
write.csv(panel, file.path(analysis_dir, "main_2024_pm25_daily_panel.csv"), row.names = FALSE)
write.csv(split, file.path(analysis_dir, "main_2024_train_test_split.csv"), row.names = FALSE)
write.csv(site_missing, file.path(analysis_dir, "main_2024_site_completeness.csv"), row.names = FALSE)
write.csv(panel_summary, file.path(analysis_dir, "main_2024_panel_summary.csv"), row.names = FALSE)
