project_root <- normalizePath(getwd(), mustWork = TRUE)

years <- 2022:2025
selection_year <- 2024
max_nodes <- 30
min_hourly_completeness <- 0.90
pm25_column_pattern <- "PM<sub>2.5</sub> particulate matter \\(Hourly measured\\)"

#create folders for AURN files
raw_dir <- file.path(project_root, "data_raw", "defra_aurn", "raw_csv")
metadata_dir <- file.path(project_root, "data_raw", "defra_aurn", "metadata")
processed_dir <- file.path(project_root, "data", "processed", "defra_aurn")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

network_url <- "https://uk-air.defra.gov.uk/networks/network-info?view=aurn"
network_html <- file.path(metadata_dir, "aurn_network_info.html")

download_if_missing <- function(url, dest) {
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    ok <- tryCatch({
      download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) warning("Failed to download: ", url)
  }
  file.exists(dest) && file.info(dest)$size > 0
}

#to find where the actual AURN data table begins
find_data_header <- function(csv_file) {
  first_lines <- readLines(csv_file, n = 20, warn = FALSE)
  header_line <- grep("^Date,time,", first_lines)[1]
  if (is.na(header_line)) stop("Could not find Date,time header in ", basename(csv_file))
  header_line - 1
}

#extract site ids and names 
extract_site_options <- function(html_file) {
  html <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  matches <- gregexpr("<option value=\"[^\"]+\">[^<]+</option>", html, perl = TRUE)
  options <- regmatches(html, matches)[[1]]
  data.frame(
    site_id = sub(".*value=\"([^\"]+)\".*", "\\1", options),
    site_name = sub(".*\">([^<]+)</option>.*", "\\1", options),
    stringsAsFactors = FALSE
  )
}

#extract coordinates & environment type for each site
parse_site_metadata <- function(site_id) {
  dest <- file.path(metadata_dir, paste0(site_id, "_site_info.html"))

  url <- paste0("https://uk-air.defra.gov.uk/networks/aurn-site-info?site_id=", site_id)
  if (!download_if_missing(url, dest)) {
    return(data.frame(site_id = site_id, latitude = NA_real_, longitude = NA_real_, environment_type = NA_character_))
  }

  text <- paste(readLines(dest, warn = FALSE), collapse = " ")
  latlon <- regexpr("Latitude/Longitude:</strong>\\s*[-0-9.]+,\\s*[-0-9.]+", text, perl = TRUE)
  if (latlon[1] > 0) {
    value <- regmatches(text, latlon)
    nums <- regmatches(value, gregexpr("[-0-9.]+", value, perl = TRUE))[[1]]
    latitude <- as.numeric(nums[1])
    longitude <- as.numeric(nums[2])
  } else {
    latitude <- NA_real_
    longitude <- NA_real_
  }

  env_match <- regexpr("Environment Type:</strong>\\s*<a[^>]*>[^<]+", text, perl = TRUE)
  environment_type <- if (env_match[1] > 0) trimws(sub(".*>", "", regmatches(text, env_match))) else NA_character_

  data.frame(
    site_id = site_id,
    latitude = latitude,
    longitude = longitude,
    environment_type = environment_type,
    stringsAsFactors = FALSE
  )
}

#download hourly pm2.5 observations for one site and year
read_site_pollutant <- function(site_id, site_name, year) {
  csv_file <- file.path(raw_dir, paste0(site_id, "_", year, ".csv"))

  url <- paste0("https://uk-air.defra.gov.uk/datastore/data_files/site_data/", site_id, "_", year, ".csv?v=1")
  if (!download_if_missing(url, csv_file)) return(NULL)

  parse_error <- NULL
  dat <- tryCatch(
    read.csv(
      csv_file,
      skip = find_data_header(csv_file),
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    ),
    error = function(e) {
      parse_error <<- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(parse_error)) {
    warning("Could not parse ", basename(csv_file), ": ", parse_error)
  }
  if (is.null(dat) || nrow(dat) == 0 || !"Date" %in% names(dat) || !"time" %in% names(dat)) return(NULL)

  dat <- dat[trimws(dat$Date) != "" & !is.na(dat$Date), , drop = FALSE]
  pollutant_col <- grep(pm25_column_pattern, names(dat), value = TRUE)
  if (length(pollutant_col) == 0) return(NULL)

  timestamp <- as.POSIXct(paste(dat$Date, dat$time), format = "%d-%m-%Y %H:%M", tz = "GMT")
  value <- suppressWarnings(as.numeric(dat[[pollutant_col[1]]]))
  status_col <- which(names(dat) == pollutant_col[1]) + 1
  status <- if (status_col <= ncol(dat)) dat[[status_col]] else NA_character_

  out <- data.frame(
    timestamp = timestamp,
    site_id = site_id,
    site_name = site_name,
    pm25 = value,
    status = status,
    stringsAsFactors = FALSE
  )
  out[!is.na(out$timestamp), , drop = FALSE]
}

make_wide_panels <- function(pm25_long, selected_ids, year) {
  hourly_wide <- reshape(
    pm25_long[, c("timestamp", "site_id", "pm25")],
    idvar = "timestamp",
    timevar = "site_id",
    direction = "wide"
  )
  names(hourly_wide) <- sub("^pm25\\.", "", names(hourly_wide))
  hourly_wide <- hourly_wide[order(hourly_wide$timestamp), , drop = FALSE]
  for (site_id in selected_ids) {
    if (!site_id %in% names(hourly_wide)) hourly_wide[[site_id]] <- NA_real_
  }
  hourly_wide <- hourly_wide[, c("timestamp", selected_ids), drop = FALSE]
  #aggregate hourly observations to daily mean PM2.5
  daily <- pm25_long[pm25_long$site_id %in% selected_ids, , drop = FALSE]
  daily$date <- as.Date(daily$timestamp - 1, tz = "GMT")
  daily_mean <- aggregate(pm25 ~ date + site_id, data = daily, FUN = function(x) mean(x, na.rm = TRUE))
  daily_wide <- reshape(daily_mean, idvar = "date", timevar = "site_id", direction = "wide")
  names(daily_wide) <- sub("^pm25\\.", "", names(daily_wide))
  daily_wide <- daily_wide[order(daily_wide$date), , drop = FALSE]
  for (site_id in selected_ids) {
    if (!site_id %in% names(daily_wide)) daily_wide[[site_id]] <- NA_real_
  }
  daily_wide <- daily_wide[, c("date", selected_ids), drop = FALSE]

  write.csv(hourly_wide, file.path(processed_dir, paste0("aurn_pm25_hourly_panel_", year, ".csv")), row.names = FALSE)
  write.csv(daily_wide, file.path(processed_dir, paste0("aurn_pm25_daily_panel_", year, ".csv")), row.names = FALSE)
  #completeness summary
  data.frame(
    year = year,
    hourly_rows = nrow(hourly_wide),
    daily_rows = nrow(daily_wide),
    nodes = length(selected_ids),
    hourly_cells = nrow(hourly_wide) * length(selected_ids),
    hourly_missing = sum(is.na(hourly_wide[, selected_ids, drop = FALSE])),
    daily_cells = nrow(daily_wide) * length(selected_ids),
    daily_missing = sum(is.na(daily_wide[, selected_ids, drop = FALSE]))
  )
}

if (!download_if_missing(network_url, network_html)) stop("Could not obtain AURN network site list.")

sites <- extract_site_options(network_html)
write.csv(sites, file.path(metadata_dir, "aurn_current_sites.csv"), row.names = FALSE)

metadata <- do.call(rbind, lapply(sites$site_id, parse_site_metadata))
metadata <- merge(sites, metadata, by = "site_id", all.x = TRUE)
write.csv(metadata, file.path(metadata_dir, "aurn_current_site_metadata.csv"), row.names = FALSE)
#download hourly PM2.5 data for every available AURN site and year
yearly_long <- list()
for (year in years) {
  pm25_list <- vector("list", nrow(sites))
  for (i in seq_len(nrow(sites))) {
    pm25_list[[i]] <- read_site_pollutant(sites$site_id[i], sites$site_name[i], year)
  }
  pm25_long <- do.call(rbind, pm25_list[!vapply(pm25_list, is.null, logical(1))])
  if (is.null(pm25_long) || nrow(pm25_long) == 0) stop("No PM2.5 data extracted for ", year)
  write.csv(pm25_long, file.path(raw_dir, paste0("aurn_pm25_hourly_long_", year, ".csv")), row.names = FALSE)
  yearly_long[[as.character(year)]] <- pm25_long
}
#use 2024 completeness to choose one fixed set of monitoring stations
selection_long <- yearly_long[[as.character(selection_year)]]
selection_wide <- reshape(
  selection_long[, c("timestamp", "site_id", "pm25")],
  idvar = "timestamp",
  timevar = "site_id",
  direction = "wide"
)
names(selection_wide) <- sub("^pm25\\.", "", names(selection_wide))
selection_wide <- selection_wide[order(selection_wide$timestamp), , drop = FALSE]
#calculate hourly completeness for each site in 2024
site_columns <- setdiff(names(selection_wide), "timestamp")
completeness <- data.frame(
  site_id = site_columns,
  observed_hours = colSums(!is.na(selection_wide[site_columns])),
  total_hours = nrow(selection_wide),
  completeness = colSums(!is.na(selection_wide[site_columns])) / nrow(selection_wide),
  stringsAsFactors = FALSE
)
completeness <- merge(completeness, metadata, by = "site_id", all.x = TRUE)
completeness <- completeness[order(-completeness$completeness), ]
write.csv(completeness, file.path(metadata_dir, paste0("aurn_pm25_completeness_", selection_year, ".csv")), row.names = FALSE)
#keep 30 sites with at least 90% completeness 
selected <- completeness[
  completeness$completeness >= min_hourly_completeness &
    !is.na(completeness$latitude) &
    !is.na(completeness$longitude),
]
selected <- selected[seq_len(min(max_nodes, nrow(selected))), ]
selected_ids <- selected$site_id
write.csv(selected, file.path(processed_dir, paste0("aurn_pm25_selected_sites_", selection_year, ".csv")), row.names = FALSE)
write.csv(selected, file.path(processed_dir, "aurn_pm25_selected_sites_fixed_2024_nodes.csv"), row.names = FALSE)

#build the same 30station panel for each year
panel_summaries <- do.call(rbind, lapply(years, function(year) {
  make_wide_panels(yearly_long[[as.character(year)]], selected_ids, year)
}))
write.csv(panel_summaries, file.path(processed_dir, "aurn_pm25_panel_summary_2022_2025.csv"), row.names = FALSE)
