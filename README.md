# Bayesian GNAR PM2.5 Forecasting

This repository contains the R code for my thesis analysis on
Bayesian GNAR forecasting of daily PM2.5 concentrations from UK Automatic Urban
and Rural Network (AURN) monitoring stations.

The code covers:
- construction of the 2024 daily PM2.5 panel and station networks;
- GNAR network and order selection;
- naive, AR(1), unrestricted VAR and ridge VAR benchmarks;
- Bayesian Gaussian GNAR and hierarchical Bayesian GNAR;
- synthetic-data validation of the Bayesian GNAR sampler;
- missing-data robustness;
- calibration and upper-tail diagnostics;
- two modifications to improve underestimation of high concentration PM2.5 (Student-t and peak-weighted GNAR);
- external-year validation using 2022-2024 for training and 2025 for testing.

## Folder structure

```text
final_code/
├── R/                  Shared modelling and evaluation functions
├── scripts/            Numbered analysis scripts
├── data/processed/     Processed PM2.5 panels and network matrices
├── outputs/tables/     Tables used in the thesis
├── outputs/figures/    Figures used in the thesis
└── docs/               Result mapping and data notes
```

## Running The Workflow

From the repository root, run:

```bash
Rscript run_all.R
```

The processed AURN daily panels are included under `data/processed`. If these
processed files are missing, `run_all.R` automatically runs
`scripts/01_extract_aurn_data.R` before the remaining analysis scripts.
