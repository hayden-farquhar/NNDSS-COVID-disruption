#!/usr/bin/env Rscript
# =============================================================================
# 05_disruption_metrics.R
#
# Compute aggregate disruption metrics for all 47 analysis-suitable diseases:
#   - Period-level mean O/E ratios and percent change
#   - Wilcoxon signed-rank tests for significance
#   - Cumulative deficit (missing cases) and excess
#   - Recovery trajectory classification
#   - Year-by-year disruption summary
#
# Inputs:
#   data/processed/annual_with_baselines.csv
#   data/processed/baselines_annual_national.csv
#   data/processed/tier1_arima_forecasts.csv
#
# Outputs:
#   data/processed/disruption_metrics_national.csv  — 47 diseases, period metrics
#   data/processed/disruption_year_by_year.csv      — 47 diseases × years
#   data/processed/disruption_monthly_tier1.csv     — 3 diseases, monthly detail
#   output/tables/disruption_summary_table.csv      — Table 2 for manuscript
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

cat("=== Phase 3: Disruption Metrics ===\n\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
abl       <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                      show_col_types = FALSE)
baselines <- read_csv(file.path(outdir, "baselines_annual_national.csv"),
                      show_col_types = FALSE)
arima     <- read_csv(file.path(outdir, "tier1_arima_forecasts.csv"),
                      show_col_types = FALSE)

# National level only for aggregate metrics
national <- abl %>% filter(level == "national")

# =============================================================================
# 1. Year-by-year disruption (already in annual_with_baselines, but add extras)
# =============================================================================
cat("--- Year-by-year disruption ---\n")

year_by_year <- national %>%
  select(disease, year, period, rate_per_100k, count,
         transmission_mode, border_sensitive, vaccine_preventable,
         count_tier, baseline_mean_rate, baseline_mean_count,
         oe_ratio, pct_change, below_pi, above_pi,
         oe_ratio_count, pct_change_count) %>%
  filter(period != "baseline") %>%
  arrange(disease, year)

cat("Year-by-year rows:", nrow(year_by_year), "\n")

# =============================================================================
# 2. Period-level aggregate metrics
# =============================================================================
cat("\n--- Period-level aggregate metrics ---\n")

# Helper: safe Wilcoxon test (handles ties and small samples)
safe_wilcox <- function(x, y) {
  if (length(x) < 2 | length(y) < 2) return(NA_real_)
  tryCatch(
    wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

period_metrics <- national %>%
  group_by(disease, transmission_mode, border_sensitive,
           vaccine_preventable, count_tier) %>%
  summarise(
    # Baseline stats (from the data directly)
    baseline_mean_rate  = mean(rate_per_100k[period == "baseline"], na.rm = TRUE),
    baseline_mean_count = mean(count[period == "baseline"], na.rm = TRUE),
    baseline_n = sum(period == "baseline"),

    # COVID acute (2020-2021)
    covid_mean_rate     = mean(rate_per_100k[period == "covid_acute"], na.rm = TRUE),
    covid_mean_count    = mean(count[period == "covid_acute"], na.rm = TRUE),
    covid_n = sum(period == "covid_acute"),
    covid_oe_rate       = covid_mean_rate / baseline_mean_rate,
    covid_pct_change    = (covid_mean_rate - baseline_mean_rate) /
                            baseline_mean_rate * 100,

    # Transition (2022)
    transition_rate     = mean(rate_per_100k[period == "transition"], na.rm = TRUE),
    transition_count    = mean(count[period == "transition"], na.rm = TRUE),
    transition_oe_rate  = transition_rate / baseline_mean_rate,
    transition_pct      = (transition_rate - baseline_mean_rate) /
                            baseline_mean_rate * 100,

    # Post-COVID (2023+)
    post_mean_rate      = mean(rate_per_100k[period == "post_covid"], na.rm = TRUE),
    post_mean_count     = mean(count[period == "post_covid"], na.rm = TRUE),
    post_n = sum(period == "post_covid"),
    post_oe_rate        = post_mean_rate / baseline_mean_rate,
    post_pct_change     = (post_mean_rate - baseline_mean_rate) /
                            baseline_mean_rate * 100,

    # Wilcoxon tests (baseline vs COVID, baseline vs post-COVID)
    wilcox_p_covid      = safe_wilcox(
      rate_per_100k[period == "covid_acute"],
      rate_per_100k[period == "baseline"]
    ),
    wilcox_p_post       = safe_wilcox(
      rate_per_100k[period == "post_covid"],
      rate_per_100k[period == "baseline"]
    ),

    # Maximum single-year reduction
    max_reduction_year  = year[which.min(oe_ratio[period != "baseline"])],
    max_reduction_oe    = min(oe_ratio[period != "baseline"], na.rm = TRUE),

    # Maximum single-year excess
    max_excess_year     = year[which.max(oe_ratio[period != "baseline"])],
    max_excess_oe       = max(oe_ratio[period != "baseline"], na.rm = TRUE),

    .groups = "drop"
  )

# --- Cumulative deficit and excess (count-based, national) ---
cumulative <- national %>%
  filter(period != "baseline", !is.na(count), !is.na(baseline_mean_count)) %>%
  mutate(
    deficit = pmax(0, baseline_mean_count - count),  # Missing cases
    excess  = pmax(0, count - baseline_mean_count)   # Extra cases
  ) %>%
  group_by(disease) %>%
  summarise(
    cumulative_deficit_total  = sum(deficit),
    cumulative_excess_total   = sum(excess),
    cumulative_deficit_covid  = sum(deficit[period == "covid_acute"]),
    cumulative_excess_covid   = sum(excess[period == "covid_acute"]),
    cumulative_deficit_post   = sum(deficit[period %in% c("transition", "post_covid")]),
    cumulative_excess_post    = sum(excess[period %in% c("transition", "post_covid")]),
    net_cumulative            = cumulative_excess_total - cumulative_deficit_total,
    .groups = "drop"
  )

period_metrics <- period_metrics %>%
  left_join(cumulative, by = "disease")

# --- Recovery trajectory classification ---
period_metrics <- period_metrics %>%
  mutate(
    # Based on most recent post-COVID O/E ratio
    trajectory = case_when(
      count_tier == "ultra_low" ~ "insufficient_data",
      post_oe_rate > 1.2  ~ "overshoot",
      post_oe_rate >= 0.9 ~ "returned",
      post_oe_rate >= 0.5 ~ "partial_recovery",
      TRUE                ~ "sustained_suppression"
    ),
    # Significance flags
    covid_significant    = !is.na(wilcox_p_covid) & wilcox_p_covid < 0.05,
    post_significant     = !is.na(wilcox_p_post) & wilcox_p_post < 0.05,
    # Direction of disruption
    covid_direction = case_when(
      covid_oe_rate < 0.9  ~ "decreased",
      covid_oe_rate > 1.1  ~ "increased",
      TRUE                  ~ "unchanged"
    )
  )

# Handle Inf/NaN
period_metrics <- period_metrics %>%
  mutate(across(where(is.numeric),
                ~ if_else(is.infinite(.) | is.nan(.), NA_real_, .)))

cat("Period metrics computed for", nrow(period_metrics), "diseases\n")

# =============================================================================
# 3. Monthly disruption detail for Tier 1 diseases
# =============================================================================
cat("\n--- Monthly disruption for Tier 1 diseases ---\n")

monthly_disruption <- arima %>%
  filter(!is.na(observed), period != "baseline") %>%
  mutate(
    deficit = pmax(0, expected - observed),
    excess  = pmax(0, observed - expected)
  ) %>%
  select(disease, year, month, period, observed, expected,
         expected_lo95, expected_hi95, oe_ratio, pct_change,
         below_pi95, above_pi95, deficit, excess)

cat("Monthly disruption rows:", nrow(monthly_disruption), "\n")

# Monthly cumulative deficit
monthly_cumulative <- monthly_disruption %>%
  arrange(disease, year, month) %>%
  group_by(disease) %>%
  mutate(
    cum_deficit = cumsum(deficit),
    cum_excess  = cumsum(excess),
    cum_net     = cum_excess - cum_deficit
  ) %>%
  ungroup()

# =============================================================================
# 4. Summary table for manuscript (Table 2)
# =============================================================================
cat("\n--- Building manuscript Table 2 ---\n")

table2 <- period_metrics %>%
  select(
    disease, transmission_mode, count_tier,
    baseline_mean_rate,
    covid_oe_rate, covid_pct_change, wilcox_p_covid, covid_significant,
    transition_oe_rate,
    post_oe_rate, post_pct_change, wilcox_p_post, post_significant,
    max_reduction_year, max_reduction_oe,
    trajectory, border_sensitive, vaccine_preventable,
    cumulative_deficit_total, cumulative_excess_total, net_cumulative
  ) %>%
  arrange(transmission_mode, covid_oe_rate)

# =============================================================================
# 5. Save outputs
# =============================================================================
write_csv(period_metrics, file.path(outdir, "disruption_metrics_national.csv"))
cat("\nSaved: disruption_metrics_national.csv (", nrow(period_metrics), "rows )\n")

write_csv(year_by_year, file.path(outdir, "disruption_year_by_year.csv"))
cat("Saved: disruption_year_by_year.csv (", nrow(year_by_year), "rows )\n")

write_csv(monthly_cumulative, file.path(outdir, "disruption_monthly_tier1.csv"))
cat("Saved: disruption_monthly_tier1.csv (", nrow(monthly_cumulative), "rows )\n")

write_csv(table2, file.path(tabdir, "disruption_summary_table.csv"))
cat("Saved: disruption_summary_table.csv (", nrow(table2), "rows )\n")

# =============================================================================
# 6. Summary output
# =============================================================================
cat("\n##############################################################\n")
cat("#              DISRUPTION METRICS SUMMARY                    #\n")
cat("##############################################################\n\n")

cat("COVID-acute disruption direction:\n")
period_metrics %>% count(covid_direction) %>% print()

cat("\nSignificant COVID disruption (Wilcoxon p < 0.05):",
    sum(period_metrics$covid_significant, na.rm = TRUE), "of",
    nrow(period_metrics), "diseases\n")

cat("\nRecovery trajectory classification:\n")
period_metrics %>% count(trajectory) %>% arrange(desc(n)) %>% print()

cat("\nTop 10 most disrupted (lowest COVID O/E):\n")
period_metrics %>%
  filter(count_tier != "ultra_low") %>%
  select(disease, transmission_mode, covid_oe_rate, covid_pct_change,
         trajectory) %>%
  arrange(covid_oe_rate) %>%
  head(10) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print(width = 120)

cat("\nTop 10 least disrupted / increased:\n")
period_metrics %>%
  filter(count_tier != "ultra_low") %>%
  select(disease, transmission_mode, covid_oe_rate, covid_pct_change,
         trajectory) %>%
  arrange(desc(covid_oe_rate)) %>%
  head(10) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print(width = 120)

cat("\nMean disruption by transmission mode:\n")
period_metrics %>%
  filter(count_tier != "ultra_low") %>%
  group_by(transmission_mode) %>%
  summarise(
    n = n(),
    mean_covid_oe = mean(covid_oe_rate, na.rm = TRUE),
    sd_covid_oe   = sd(covid_oe_rate, na.rm = TRUE),
    mean_post_oe  = mean(post_oe_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_covid_oe) %>%
  mutate(across(where(is.numeric) & !c(n), ~ round(., 3))) %>%
  print()

cat("\nCumulative deficit — top 10 by absolute missing cases:\n")
period_metrics %>%
  select(disease, cumulative_deficit_total, cumulative_excess_total,
         net_cumulative) %>%
  arrange(desc(cumulative_deficit_total)) %>%
  head(10) %>%
  mutate(across(where(is.numeric), ~ round(.))) %>%
  print()

cat("\n=== Disruption metrics complete ===\n")
