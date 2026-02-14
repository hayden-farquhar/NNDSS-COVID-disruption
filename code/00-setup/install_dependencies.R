#!/usr/bin/env Rscript
# ============================================================================
# Install all R package dependencies for NNDSS COVID disruption analysis
# Run this script once before executing the analysis pipeline.
# Tested with R 4.4.1
# ============================================================================

packages <- c(
  # Tidyverse ecosystem (data manipulation, plotting, I/O)
  "tidyverse",     # dplyr, ggplot2, tidyr, readr, stringr, forcats, lubridate, purrr
  "readxl",        # Excel imports
  "scales",        # Axis formatting
  "patchwork",     # Multi-panel figure composition

  # Web scraping / acquisition
  "httr",          # HTTP requests
  "rvest",         # HTML parsing
  "jsonlite",      # JSON parsing
  "cli",           # Console formatting

  # Time series
  "forecast",      # ARIMA modelling
  "tseries",       # Stationarity tests

  # Statistical modelling
  "nlme",          # GLS with AR(1) correlation

  "lme4",          # Mixed-effects models
  "lmtest",        # Coefficient tests (coeftest)
  "sandwich",      # Heteroscedasticity / autocorrelation robust SEs
  "car",           # VIF diagnostics
  "broom",         # Model tidying
  "broom.mixed",   # Tidying for lme4 models
  "mgcv",          # Generalised additive models

  # Clustering
  "cluster",       # k-medoids (PAM)
  "factoextra",    # Silhouette plots

  # Survival analysis
  "survival",      # Cox PH, Kaplan-Meier
  "survminer",     # KM plot formatting

  # Bootstrap
  "boot"           # Bootstrapped confidence intervals
)

cat("Checking and installing", length(packages), "packages...\n\n")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(pkg, "OK\n")
  }
}

cat("\nAll dependencies installed.\n")
cat("R version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
