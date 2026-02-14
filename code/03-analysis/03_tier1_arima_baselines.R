#!/usr/bin/env Rscript
# =============================================================================
# 03_tier1_arima_baselines.R
#
# Tier 1 baselines: Seasonal ARIMA on monthly data for 3 diseases with
# sub-annual resolution (influenza, meningococcal, salmonella).
# Pneumococcal excluded — record-level data is annual only.
#
# Inputs:
#   data/processed/tier1_monthly_national.csv
#
# Outputs:
#   data/processed/tier1_arima_forecasts.csv       — monthly expected + O/E
#   data/processed/tier1_arima_diagnostics.csv      — model summary per disease
#   data/processed/tier1_arima_annual_summary.csv   — annualised for cross-val
#   output/figures/arima_diagnostics_{disease}.png   — 3 diagnostic plots
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forecast)
  library(tseries)
  library(ggplot2)
})

outdir  <- "data/processed"
figdir  <- "output/figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Tier 1: ARIMA Baselines ===\n\n")

# -----------------------------------------------------------------------------
# Load monthly national data
# -----------------------------------------------------------------------------
monthly <- read_csv(file.path(outdir, "tier1_monthly_national.csv"),
                    show_col_types = FALSE)

diseases <- unique(monthly$disease)
cat("Diseases available:", paste(diseases, collapse = "; "), "\n\n")

# Baseline period: Jan 2015 – Feb 2020 (62 months)
baseline_start_year  <- 2015
baseline_start_month <- 1
baseline_end_year    <- 2020
baseline_end_month   <- 2

# Forecast period: Mar 2020 – Dec 2024
forecast_start_year  <- 2020
forecast_start_month <- 3

# =============================================================================
# Fit ARIMA per disease
# =============================================================================
all_forecasts   <- list()
all_diagnostics <- list()

for (dis in diseases) {
  cat("--- Fitting:", dis, "---\n")

  dis_data <- monthly %>%
    filter(disease == dis) %>%
    arrange(year, month)

  # Extract baseline period
  bl <- dis_data %>%
    filter((year > baseline_start_year |
              (year == baseline_start_year & month >= baseline_start_month)) &
           (year < baseline_end_year |
              (year == baseline_end_year & month <= baseline_end_month)))

  cat("  Baseline months:", nrow(bl), "\n")

  if (nrow(bl) < 24) {
    cat("  SKIP: insufficient baseline data\n\n")
    next
  }

  # Create time series object
  y <- ts(bl$count, start = c(baseline_start_year, baseline_start_month),
          frequency = 12)

  # Fit seasonal ARIMA with Box-Cox transformation for variance stabilisation
  fit <- auto.arima(y, seasonal = TRUE, stepwise = FALSE,
                    approximation = FALSE, lambda = "auto",
                    trace = FALSE)

  cat("  Model:", capture.output(print(fit))[2], "\n")
  cat("  Lambda:", round(fit$lambda, 3), "\n")

  # --- Diagnostics ---
  resid <- residuals(fit)

  # Ljung-Box test (no residual autocorrelation)
  lb_test <- Box.test(resid, lag = min(24, length(resid) - 1), type = "Ljung-Box")

  # Shapiro-Wilk normality test on residuals
  sw_test <- if (length(resid) >= 3 & length(resid) <= 5000) {
    shapiro.test(resid)
  } else {
    list(statistic = NA, p.value = NA)
  }

  # ADF stationarity test on residuals
  adf_test <- tryCatch(
    adf.test(resid),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  diagnostics <- tibble(
    disease           = dis,
    arima_order       = paste(fit$arma[1], fit$arma[6], fit$arma[2], sep = ","),
    seasonal_order    = paste(fit$arma[3], fit$arma[7], fit$arma[4], sep = ","),
    lambda            = fit$lambda,
    aic               = fit$aic,
    bic               = fit$bic,
    ljung_box_stat    = lb_test$statistic,
    ljung_box_p       = lb_test$p.value,
    ljung_box_pass    = lb_test$p.value > 0.05,
    shapiro_wilk_stat = sw_test$statistic,
    shapiro_wilk_p    = sw_test$p.value,
    adf_stat          = adf_test$statistic,
    adf_p             = adf_test$p.value,
    baseline_months   = nrow(bl)
  )

  all_diagnostics[[dis]] <- diagnostics

  cat("  Ljung-Box p =", round(lb_test$p.value, 4),
      if (lb_test$p.value > 0.05) "(PASS)" else "(FAIL)", "\n")

  # --- Forecast ---
  # Calculate months to forecast: Mar 2020 to end of available data
  last_year  <- max(dis_data$year)
  last_month <- max(dis_data$month[dis_data$year == last_year])

  # Total months from Mar 2020 to end of data
  h <- (last_year - forecast_start_year) * 12 +
    (last_month - forecast_start_month) + 1
  cat("  Forecasting", h, "months (Mar 2020 to",
      paste0(last_year, "/", last_month), ")\n")

  fc <- forecast(fit, h = h, level = c(80, 95))

  # Build forecast dataframe
  fc_months <- seq(
    from = as.Date(paste(forecast_start_year, forecast_start_month, 1, sep = "-")),
    by = "month",
    length.out = h
  )

  fc_df <- tibble(
    disease          = dis,
    year             = as.integer(format(fc_months, "%Y")),
    month            = as.integer(format(fc_months, "%m")),
    expected         = as.numeric(fc$mean),
    expected_lo80    = as.numeric(fc$lower[, 1]),
    expected_hi80    = as.numeric(fc$upper[, 1]),
    expected_lo95    = as.numeric(fc$lower[, 2]),
    expected_hi95    = as.numeric(fc$upper[, 2])
  )

  # Floor at 0 (counts can't be negative)
  fc_df <- fc_df %>%
    mutate(across(starts_with("expected"), ~ pmax(0, .)))

  # Join observed data
  fc_df <- fc_df %>%
    left_join(
      dis_data %>% select(disease, year, month, count),
      by = c("disease", "year", "month")
    ) %>%
    rename(observed = count) %>%
    mutate(
      oe_ratio    = observed / expected,
      pct_change  = (observed - expected) / expected * 100,
      below_pi95  = observed < expected_lo95,
      above_pi95  = observed > expected_hi95,
      period = case_when(
        year < 2022 ~ "covid_acute",
        year == 2022 ~ "transition",
        TRUE ~ "post_covid"
      )
    )

  # Also include the baseline period for completeness
  bl_df <- tibble(
    disease     = dis,
    year        = bl$year,
    month       = bl$month,
    expected    = as.numeric(fitted(fit)),
    observed    = bl$count,
    period      = "baseline"
  ) %>%
    mutate(
      oe_ratio   = observed / expected,
      pct_change = (observed - expected) / expected * 100
    )

  all_forecasts[[dis]] <- bind_rows(bl_df, fc_df)

  # --- Diagnostic plots ---
  dis_clean <- gsub("[^a-zA-Z0-9]", "_", tolower(dis))

  png(file.path(figdir, paste0("arima_diagnostics_", dis_clean, ".png")),
      width = 12, height = 10, units = "in", res = 300)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # 1. Observed vs Fitted + Forecast
  plot(fc, main = paste(dis, "- ARIMA Forecast"),
       xlab = "Year", ylab = "Monthly count")

  # Overlay observed post-baseline data
  obs_post <- dis_data %>%
    filter(year >= forecast_start_year |
             (year == forecast_start_year & month >= forecast_start_month))
  if (nrow(obs_post) > 0) {
    obs_ts <- ts(obs_post$count,
                 start = c(forecast_start_year, forecast_start_month),
                 frequency = 12)
    lines(obs_ts, col = "red", lwd = 1.5)
    legend("topleft", legend = c("Forecast", "Observed"),
           col = c("blue", "red"), lwd = 1.5, cex = 0.8)
  }

  # 2. Residual ACF
  Acf(resid, main = "ACF of Residuals")

  # 3. Residual histogram
  hist(resid, breaks = 30, main = "Residual Distribution",
       xlab = "Residual", col = "lightblue", border = "white")

  # 4. QQ plot
  qqnorm(resid, main = "Normal Q-Q Plot")
  qqline(resid, col = "red")

  dev.off()
  cat("  Saved: arima_diagnostics_", dis_clean, ".png\n\n")
}

# =============================================================================
# Combine and save
# =============================================================================
forecasts <- bind_rows(all_forecasts) %>%
  mutate(across(c(oe_ratio, pct_change), ~ if_else(is.infinite(.) | is.nan(.), NA_real_, .)))

diagnostics <- bind_rows(all_diagnostics)

# Annualised summary for Tier 2 cross-validation
annual_summary <- forecasts %>%
  filter(!is.na(observed), period != "baseline") %>%
  group_by(disease, year) %>%
  summarise(
    observed_annual    = sum(observed),
    expected_annual    = sum(expected),
    expected_lo95_annual = sum(expected_lo95, na.rm = TRUE),
    expected_hi95_annual = sum(expected_hi95, na.rm = TRUE),
    oe_ratio_annual    = sum(observed) / sum(expected),
    months_available   = n(),
    .groups = "drop"
  ) %>%
  filter(months_available == 12)  # Only complete years

# =============================================================================
# Save outputs
# =============================================================================
write_csv(forecasts, file.path(outdir, "tier1_arima_forecasts.csv"))
cat("Saved: tier1_arima_forecasts.csv (", nrow(forecasts), "rows )\n")

write_csv(diagnostics, file.path(outdir, "tier1_arima_diagnostics.csv"))
cat("Saved: tier1_arima_diagnostics.csv (", nrow(diagnostics), "rows )\n")

write_csv(annual_summary, file.path(outdir, "tier1_arima_annual_summary.csv"))
cat("Saved: tier1_arima_annual_summary.csv (", nrow(annual_summary), "rows )\n")

# Quick validation
cat("\n=== Quick Validation ===\n")

diagnostics %>%
  select(disease, arima_order, seasonal_order, lambda,
         ljung_box_p, ljung_box_pass) %>%
  print()

cat("\nAnnual O/E ratios (ARIMA-based):\n")
annual_summary %>%
  select(disease, year, observed_annual, expected_annual, oe_ratio_annual) %>%
  filter(year %in% c(2020, 2021, 2022, 2023, 2024)) %>%
  print(n = 30)

cat("\n=== Tier 1 ARIMA baselines complete ===\n")
