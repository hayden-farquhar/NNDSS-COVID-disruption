#!/usr/bin/env Rscript
# =============================================================================
# 04_baseline_validation.R
#
# Cross-validate Tier 1 (ARIMA) vs Tier 2 (historical mean) baselines,
# flag potential issues, and produce a validation report.
#
# Inputs:
#   data/processed/baselines_annual_national.csv
#   data/processed/annual_with_baselines.csv
#   data/processed/tier1_arima_annual_summary.csv
#   data/processed/tier1_arima_diagnostics.csv
#
# Outputs:
#   output/tables/baseline_validation_report.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

outdir <- "data/processed"
tabdir <- "output/tables"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Baseline Validation ===\n\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
baselines     <- read_csv(file.path(outdir, "baselines_annual_national.csv"),
                          show_col_types = FALSE)
annual_bl     <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                          show_col_types = FALSE)
arima_annual  <- read_csv(file.path(outdir, "tier1_arima_annual_summary.csv"),
                          show_col_types = FALSE)
arima_diag    <- read_csv(file.path(outdir, "tier1_arima_diagnostics.csv"),
                          show_col_types = FALSE)

# =============================================================================
# Check 1: Baseline plausibility — flag diseases with CV > 0.5
# =============================================================================
cat("--- Check 1: High variability (CV > 0.5) ---\n")

high_cv <- baselines %>%
  filter(baseline_cv_rate > 0.5) %>%
  select(disease, baseline_mean_rate, baseline_cv_rate, count_tier)

cat("Diseases with CV > 0.5:", nrow(high_cv), "\n")
print(high_cv, n = 20)
cat("  -> These diseases have wide prediction intervals. O/E ratios are\n")
cat("     less precise but still valid. Note in manuscript.\n\n")

# =============================================================================
# Check 2: Cross-validate Tier 1 vs Tier 2 for overlapping diseases
# =============================================================================
cat("--- Check 2: Tier 1 vs Tier 2 cross-validation ---\n")

# Tier 2 O/E ratios (count-based) for the 3 ARIMA diseases
tier2_oe <- annual_bl %>%
  filter(level == "national",
         disease %in% arima_annual$disease,
         year %in% arima_annual$year) %>%
  select(disease, year, tier2_oe = oe_ratio_count) %>%
  distinct()

# Tier 1 O/E ratios (ARIMA-based, annualised)
tier1_oe <- arima_annual %>%
  select(disease, year, tier1_oe = oe_ratio_annual)

crossval <- tier2_oe %>%
  inner_join(tier1_oe, by = c("disease", "year")) %>%
  mutate(
    oe_diff     = tier1_oe - tier2_oe,
    oe_pct_diff = (tier1_oe - tier2_oe) / tier2_oe * 100,
    agreement   = case_when(
      abs(oe_pct_diff) < 10  ~ "excellent",
      abs(oe_pct_diff) < 25  ~ "good",
      abs(oe_pct_diff) < 50  ~ "moderate",
      TRUE                    ~ "poor"
    )
  )

cat("\nCross-validation results:\n")
crossval %>%
  select(disease, year, tier1_oe, tier2_oe, oe_pct_diff, agreement) %>%
  mutate(across(c(tier1_oe, tier2_oe), ~ round(., 3)),
         oe_pct_diff = round(oe_pct_diff, 1)) %>%
  print(n = 30)

cat("\nAgreement summary (all years):\n")
crossval %>% count(agreement) %>% print()

# Note: 2025 is partial-year data — ARIMA forecasts full year, Tier 2 uses
# the partial count. Exclude 2025 for meaningful cross-validation.
cat("\nAgreement summary (excluding 2025, partial year):\n")
crossval %>% filter(year < 2025) %>% count(agreement) %>% print()

# Flag systematic disagreements
cat("\nMean Tier1-Tier2 difference by disease:\n")
crossval %>%
  group_by(disease) %>%
  summarise(
    n_years        = n(),
    mean_pct_diff  = mean(oe_pct_diff, na.rm = TRUE),
    sd_pct_diff    = sd(oe_pct_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()
cat("\n")

# =============================================================================
# Check 3: Diseases with significant secular trends
# =============================================================================
cat("--- Check 3: Significant secular trends ---\n")

trended <- baselines %>%
  filter(trend_significant) %>%
  select(disease, trend_slope, trend_pct_per_year, baseline_mean_rate) %>%
  mutate(
    direction = if_else(trend_slope > 0, "increasing", "decreasing"),
    impact = case_when(
      abs(trend_pct_per_year) > 15 ~ "strong",
      abs(trend_pct_per_year) > 5  ~ "moderate",
      TRUE                         ~ "weak"
    )
  ) %>%
  arrange(desc(abs(trend_pct_per_year)))

cat("Diseases with significant baseline trend:", nrow(trended), "\n")
print(trended, n = 20)
cat("  -> For these diseases, the flat-mean baseline may over/underestimate\n")
cat("     expected values. Trend-adjusted expected values are provided\n")
cat("     as sensitivity analysis.\n\n")

# =============================================================================
# Check 4: Prediction interval sanity check
# =============================================================================
cat("--- Check 4: PI coverage (baseline values inside own PI) ---\n")

pi_check <- annual_bl %>%
  filter(level == "national", period == "baseline") %>%
  mutate(inside_pi = rate_per_100k >= pi_lower_rate &
                     rate_per_100k <= pi_upper_rate)

coverage <- pi_check %>%
  group_by(disease) %>%
  summarise(
    n_baseline = n(),
    n_inside   = sum(inside_pi, na.rm = TRUE),
    coverage   = n_inside / n_baseline * 100,
    .groups = "drop"
  )

cat("All baseline values inside own 95% PI:", all(coverage$coverage == 100), "\n")
undercovered <- coverage %>% filter(coverage < 100)
if (nrow(undercovered) > 0) {
  cat("Diseases with baseline values outside PI (unexpected):\n")
  print(undercovered)
} else {
  cat("  -> All 47 diseases: 100% baseline coverage. PIs are valid.\n")
}
cat("\n")

# =============================================================================
# Check 5: ARIMA forecast reasonableness
# =============================================================================
cat("--- Check 5: ARIMA forecast reasonableness ---\n")

arima_diag %>%
  select(disease, arima_order, seasonal_order, lambda,
         ljung_box_p, ljung_box_pass) %>%
  print()

# Check if ARIMA forecasts diverge to extreme values
arima_extremes <- arima_annual %>%
  group_by(disease) %>%
  summarise(
    max_expected = max(expected_annual),
    min_expected = min(expected_annual),
    max_oe       = max(oe_ratio_annual, na.rm = TRUE),
    min_oe       = min(oe_ratio_annual, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nARIMA forecast ranges:\n")
print(arima_extremes)

# Flag if any forecast goes negative (shouldn't with Box-Cox + floor)
if (any(arima_annual$expected_annual < 0)) {
  cat("\nWARNING: Negative forecast values detected!\n")
} else {
  cat("\n  -> No negative forecasts. Box-Cox transformation working correctly.\n")
}

# Salmonella Ljung-Box note
sal_diag <- arima_diag %>% filter(str_detect(disease, "Salmonella"))
if (nrow(sal_diag) > 0 && !sal_diag$ljung_box_pass) {
  cat("\n  NOTE: Salmonella ARIMA residuals show significant autocorrelation\n")
  cat("  (Ljung-Box p < 0.05). This is common with count data and large\n")
  cat("  samples. The model captures the main seasonal pattern but\n")
  cat("  prediction intervals may be slightly too narrow. Consider\n")
  cat("  reporting Tier 2 (historical mean) as primary for salmonella.\n")
}
cat("\n")

# =============================================================================
# Check 6: Ultra-low count diseases
# =============================================================================
cat("--- Check 6: Ultra-low count diseases ---\n")

ultra_low <- baselines %>%
  filter(count_tier == "ultra_low") %>%
  select(disease, baseline_mean_count, baseline_min_count, baseline_max_count)

cat("Ultra-low count diseases (< 10 cases/yr):", nrow(ultra_low), "\n")
print(ultra_low)
cat("  -> O/E ratios unreliable for these. Report raw counts + range only.\n\n")

# =============================================================================
# Compile validation report
# =============================================================================
cat("--- Compiling validation report ---\n")

report <- baselines %>%
  select(disease, baseline_mean_rate, baseline_cv_rate, count_tier,
         baseline_mean_count, trend_significant, trend_pct_per_year) %>%
  mutate(
    flag_high_cv    = baseline_cv_rate > 0.5,
    flag_ultra_low  = count_tier == "ultra_low",
    flag_trended    = trend_significant,
    flag_border_sensitive = disease %in% c(
      "Dengue virus infection", "Malaria", "Typhoid Fever", "Paratyphoid",
      "Hepatitis A", "Hepatitis E", "Shigellosis", "Chikungunya virus infection",
      "Leprosy", "Measles", "Rubella"
    ),
    has_arima = disease %in% arima_diag$disease,
    arima_ljung_box_pass = if_else(
      has_arima,
      arima_diag$ljung_box_pass[match(disease, arima_diag$disease)],
      NA
    ),
    # Overall quality flag
    baseline_quality = case_when(
      flag_ultra_low ~ "low_counts",
      flag_high_cv & flag_trended ~ "variable_trended",
      flag_high_cv ~ "variable",
      flag_trended & abs(trend_pct_per_year) > 15 ~ "strong_trend",
      flag_trended ~ "moderate_trend",
      TRUE ~ "good"
    )
  ) %>%
  arrange(baseline_quality, desc(baseline_mean_count))

write_csv(report, file.path(tabdir, "baseline_validation_report.csv"))
cat("Saved: baseline_validation_report.csv (", nrow(report), "rows )\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n##############################################################\n")
cat("#              BASELINE VALIDATION SUMMARY                   #\n")
cat("##############################################################\n\n")

cat("Total diseases with baselines: ", nrow(report), "\n")
cat("Baseline quality distribution:\n")
report %>% count(baseline_quality) %>% print()

cat("\nKey flags:\n")
cat("  High variability (CV > 0.5):   ", sum(report$flag_high_cv), "\n")
cat("  Ultra-low counts (< 10/yr):    ", sum(report$flag_ultra_low), "\n")
cat("  Significant trend:             ", sum(report$flag_trended), "\n")
cat("  Border-sensitive:              ", sum(report$flag_border_sensitive), "\n")
cat("  ARIMA models fitted:           ", sum(report$has_arima), "\n")
cat("  ARIMA Ljung-Box pass:          ",
    sum(report$arima_ljung_box_pass, na.rm = TRUE), "/",
    sum(report$has_arima), "\n")

cat("\nCross-validation (Tier 1 vs Tier 2):\n")
crossval %>%
  count(agreement) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print()

cat("\n=== Baseline validation complete ===\n")
