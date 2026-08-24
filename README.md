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


