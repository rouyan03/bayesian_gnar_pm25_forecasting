# Data notes

## Main 2024 panel

`data/processed/analysis/main_2024_pm25_daily_panel.csv` contains the complete
daily PM2.5 panel used for the main 2024 experiment. It has 30 AURN stations
and 366 days.

## 2022-2025 panels

`data/processed/defra_aurn/aurn_pm25_daily_panel_2022.csv` to
`aurn_pm25_daily_panel_2025.csv` contain daily PM2.5 panels used for the
external-year validation. The same 30-station cohort is used.

## Networks

The station network matrices are stored in `data/processed/networks`.

- `aurn_network_equal_neighbour_2024_nodes.csv`
- `aurn_network_inverse_distance_2024_nodes.csv`
- `aurn_network_knn4_2024_nodes.csv`
- `aurn_network_threshold150km_2024_nodes.csv`

Rows are station origins and columns are station neighbours. Diagonal entries
are zero. Non-empty rows are row-normalised.

