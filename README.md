# Flux Time-Series ML Pipeline

An end-to-end R pipeline for processing, quality-controlling, and modeling high-frequency environmental flux data (carbon, water, and energy). Designed for Eddy Covariance datasets but adaptable to any environmental time-series.

## Overview

This project demonstrates a full analytical workflow:

```
Raw sensor data → QC & gap-filling → Exploratory analysis → ML modeling → Visualization
```

## Features

- **Automated data QC** — spike detection, range checks, u* filtering
- **Gap-filling** — mean diurnal variation (MDV) and Random Forest gap-filling
- **Machine learning models** — Random Forest and ANN for flux prediction
- **Explainable AI** — SHAP-based feature importance
- **Publication-quality plots** — ggplot2 visualizations

## Repository Structure

```
flux-timeseries-ml/
├── 01_data_processing.R      # Data loading, QC, and gap-filling
├── 02_exploratory_analysis.R # Summary stats, seasonal patterns, visualizations
├── 03_machine_learning.R     # Random Forest + ANN modeling with SHAP
├── data/
│   └── sample_flux_data.csv  # Sample synthetic dataset for demo
├── output/
│   └── figures/              # Generated plots saved here
└── README.md
```

## Requirements

```r
install.packages(c("tidyverse", "ggplot2", "lubridate", "randomForest",
                   "caret", "shapr", "patchwork", "scales"))
```

## Quick Start

```r
# Step 1: Process and QC the data
source("01_data_processing.R")

# Step 2: Explore patterns
source("02_exploratory_analysis.R")

# Step 3: Run ML models
source("03_machine_learning.R")
```

## Key Results

The Random Forest model achieves **R² > 0.90** for daytime carbon flux prediction, outperforming traditional empirical models. SHAP analysis identifies solar radiation and temperature as dominant drivers, consistent with process-based understanding.

## Background

This workflow is based on methods used in peer-reviewed research conducted at the Iowa State SABR research site as part of the DOE-funded CABBI project. See [publications](#) for full context.

## Author

**Rojda Guler Aslan Sungur, Ph.D.**  
Research Scientist III, Iowa State University  
rojdaaslan@gmail.com
