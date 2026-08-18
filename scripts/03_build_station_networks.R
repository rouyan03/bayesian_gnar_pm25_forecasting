project_root <- normalizePath(getwd(), mustWork = TRUE)

processed_dir <- file.path(project_root, "data", "processed", "defra_aurn")
network_dir <- file.path(project_root, "data", "processed", "networks")
dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)

sites_file <- file.path(processed_dir, "aurn_pm25_selected_sites_fixed_2024_nodes.csv")
if (!file.exists(sites_file)) {
  stop("Missing selected-site file. Run final_code/scripts/01_extract_aurn_data.R first.")
}

sites <- read.csv(sites_file, stringsAsFactors = FALSE)
site_ids <- sites$site_id

#great-circle distance between the stations
haversine_km <- function(lat1, lon1, lat2, lon2) {
  radius <- 6371 #earth's radius
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * radius * atan2(sqrt(a), sqrt(1 - a))
}

row_normalise <- function(A) {
  row_sums <- rowSums(A)
  W <- A
  W[] <- 0
  keep <- row_sums > 0
  W[keep, ] <- A[keep, , drop = FALSE] / row_sums[keep]
  W
}

#summary in chapter 2
network_stats <- function(A, name) {
  nonzero <- A > 0
  diag(nonzero) <- FALSE
  degrees <- rowSums(nonzero)
  data.frame(
    network = name,
    nodes = nrow(A),
    directed_nonzero_edges = sum(nonzero),
    density_directed_no_self = sum(nonzero) / (nrow(A) * (nrow(A) - 1)),
    min_out_degree = min(degrees),
    mean_out_degree = mean(degrees),
    max_out_degree = max(degrees),
    row_normalised = all(abs(rowSums(A)[rowSums(A) > 0] - 1) < 1e-8),
    stringsAsFactors = FALSE
  )
}

n <- nrow(sites)
distance <- matrix(0, n, n, dimnames = list(site_ids, site_ids))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    distance[i, j] <- haversine_km(
      sites$latitude[i], sites$longitude[i],
      sites$latitude[j], sites$longitude[j]
    )
  }
}
#inv distance
inverse_distance <- 1 / distance
diag(inverse_distance) <- 0
inverse_distance[!is.finite(inverse_distance)] <- 0
inverse_distance <- row_normalise(inverse_distance)

#knn4
k <- 4
knn4 <- matrix(0, n, n, dimnames = list(site_ids, site_ids))
for (i in seq_len(n)) {
  nearest <- order(distance[i, ])[2:(k + 1)]
  knn4[i, nearest] <- 1
}
knn4 <- row_normalise(knn4)

#150km
threshold_km <- 150
threshold150km <- ifelse(distance > 0 & distance <= threshold_km, 1, 0)
threshold150km <- row_normalise(threshold150km)
dimnames(threshold150km) <- list(site_ids, site_ids)

#equal neighbour
equal_neighbour <- matrix(1, n, n, dimnames = list(site_ids, site_ids))
diag(equal_neighbour) <- 0
equal_neighbour <- row_normalise(equal_neighbour)

write.csv(distance, file.path(network_dir, "aurn_station_distance_km_2024_nodes.csv"))
write.csv(inverse_distance, file.path(network_dir, "aurn_network_inverse_distance_2024_nodes.csv"))
write.csv(knn4, file.path(network_dir, "aurn_network_knn4_2024_nodes.csv"))
write.csv(threshold150km, file.path(network_dir, "aurn_network_threshold150km_2024_nodes.csv"))
write.csv(equal_neighbour, file.path(network_dir, "aurn_network_equal_neighbour_2024_nodes.csv"))

network_summary <- do.call(rbind, list(
  network_stats(inverse_distance, "inverse_distance"),
  network_stats(knn4, "knn4"),
  network_stats(threshold150km, "threshold150km"),
  network_stats(equal_neighbour, "equal_neighbour")
))
write.csv(network_summary, file.path(network_dir, "aurn_network_summary_2024_nodes.csv"), row.names = FALSE)
