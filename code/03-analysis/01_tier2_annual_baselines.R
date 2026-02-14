#!/usr/bin/env Rscript
# =============================================================================
# 01_tier2_annual_baselines.R
#
# Tier 2 baselines: historical range from 2015-2019 annual data for all 47
# analysis-suitable diseases. Computes mean, SD, prediction intervals,
# baseline trends, O/E ratios, and count tiers.
#
# Inputs:
#   data/processed/nndss_annual_analysis.csv
#   data/processed/disease_metadata.csv
#
# Outputs:
#   data/processed/baselines_annual_national.csv     — 47 rows, one per disease
#   data/processed/baselines_annual_by_state.csv     — 47 × 8 state-level
#   data/processed/annual_with_baselines.csv         — full time series + O/E
#   data/processed/trend_adjusted_expected.csv       — trend-adjusted subset
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(broom)
})

outdir <- "data/processed"

cat("=== Tier 2: Annual Baselines ===\n\n")

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
annual <- read_csv(file.path(outdir, "nndss_annual_analysis.csv"),
                   show_col_types = FALSE)
meta   <- read_csv(file.path(outdir, "disease_metadata.csv"),
                   show_col_types = FALSE)

# Filter to analysis-suitable diseases
suitable_diseases <- meta %>%
  filter(analysis_suitable == TRUE) %>%
  pull(disease)

cat("Analysis-suitable diseases:", length(suitable_diseases), "\n")

# =============================================================================
# 2. National baselines (AUS level, 2015-2019)
# =============================================================================
cat("\n--- Computing national baselines ---\n")

baseline_national <- annual %>%
  filter(disease %in% suitable_diseases,
         state == "AUS",
         period == "baseline") %>%
  group_by(disease) %>%
  summarise(
    baseline_n_years     = n(),
    baseline_mean_rate   = mean(rate_per_100k, na.rm = TRUE),
    baseline_sd_rate     = sd(rate_per_100k, na.rm = TRUE),
    baseline_cv_rate     = baseline_sd_rate / baseline_mean_rate,
    baseline_median_rate = median(rate_per_100k, na.rm = TRUE),
    baseline_min_rate    = min(rate_per_100k, na.rm = TRUE),
    baseline_max_rate    = max(rate_per_100k, na.rm = TRUE),
    baseline_mean_count  = mean(count, na.rm = TRUE),
    baseline_sd_count    = sd(count, na.rm = TRUE),
    baseline_median_count = median(count, na.rm = TRUE),
    baseline_min_count   = min(count, na.rm = TRUE),
    baseline_max_count   = max(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # 95% prediction interval: mean ± t(0.975, n-1) * SD * sqrt(1 + 1/n)
    t_crit = qt(0.975, baseline_n_years - 1),
    pi_lower_rate = pmax(0, baseline_mean_rate -
                           t_crit * baseline_sd_rate *
                           sqrt(1 + 1 / baseline_n_years)),
    pi_upper_rate = baseline_mean_rate +
                      t_crit * baseline_sd_rate *
                      sqrt(1 + 1 / baseline_n_years),
    pi_lower_count = pmax(0, baseline_mean_count -
                            t_crit * baseline_sd_count *
                            sqrt(1 + 1 / baseline_n_years)),
    pi_upper_count = baseline_mean_count +
                       t_crit * baseline_sd_count *
                       sqrt(1 + 1 / baseline_n_years),
    # Count tiers
    count_tier = case_when(
      baseline_mean_count < 10  ~ "ultra_low",
      baseline_mean_count < 50  ~ "low",
      baseline_mean_count < 200 ~ "moderate",
      TRUE                      ~ "high"
    )
  ) %>%
  select(-t_crit)

# --- Baseline trends (linear rate ~ year) ---
cat("Computing baseline trends...\n")

trends <- annual %>%
  filter(disease %in% suitable_diseases,
         state == "AUS",
         period == "baseline") %>%
  group_by(disease) %>%
  do({
    mod <- lm(rate_per_100k ~ year, data = .)
    tidy_mod <- tidy(mod) %>% filter(term == "year")
    tibble(
      trend_slope       = tidy_mod$estimate,
      trend_se          = tidy_mod$std.error,
      trend_p_value     = tidy_mod$p.value,
      trend_significant = tidy_mod$p.value < 0.05,
      # Percent change per year relative to mean
      trend_pct_per_year = tidy_mod$estimate / mean(.$rate_per_100k) * 100
    )
  }) %>%
  ungroup()

baseline_national <- baseline_national %>%
  left_join(trends, by = "disease")

cat("Diseases with significant trend:", sum(baseline_national$trend_significant), "\n")
baseline_national %>%
  filter(trend_significant) %>%
  select(disease, baseline_mean_rate, trend_slope, trend_pct_per_year, trend_p_value) %>%
  arrange(trend_p_value) %>%
  print(n = 20)

# =============================================================================
# 3. Trend-adjusted expected values (sensitivity analysis)
# =============================================================================
cat("\n--- Computing trend-adjusted expected values ---\n")

trended_diseases <- baseline_national %>%
  filter(trend_significant) %>%
  pull(disease)

if (length(trended_diseases) > 0) {
  trend_adjusted <- annual %>%
    filter(disease %in% trended_diseases,
           state == "AUS",
           period == "baseline") %>%
    group_by(disease) %>%
    do({
      mod <- lm(rate_per_100k ~ year, data = .)
      new_data <- tibble(year = 2020:2025)
      preds <- predict(mod, newdata = new_data, interval = "prediction")
      new_data %>%
        mutate(
          expected_rate_trend  = pmax(0, preds[, "fit"]),
          expected_lower_trend = pmax(0, preds[, "lwr"]),
          expected_upper_trend = pmax(0, preds[, "upr"])
        )
    }) %>%
    ungroup()

  # Join observed data
  trend_adjusted <- trend_adjusted %>%
    left_join(
      annual %>%
        filter(state == "AUS") %>%
        select(disease, year, rate_per_100k, count, period),
      by = c("disease", "year")
    ) %>%
    mutate(
      oe_ratio_trend = rate_per_100k / expected_rate_trend,
      pct_change_trend = (rate_per_100k - expected_rate_trend) /
                           expected_rate_trend * 100
    )

  write_csv(trend_adjusted, file.path(outdir, "trend_adjusted_expected.csv"))
  cat("Saved: trend_adjusted_expected.csv (", nrow(trend_adjusted), "rows )\n")
} else {
  cat("No diseases with significant trends — skipping trend-adjusted file\n")
}

# =============================================================================
# 4. State-level baselines (8 states, exclude NSW for varicella)
# =============================================================================
cat("\n--- Computing state-level baselines ---\n")

varicella_diseases <- c("Varicella zoster (chickenpox)",
                        "Varicella zoster (shingles)",
                        "Varicella zoster (unspecified)")

baseline_state <- annual %>%
  filter(disease %in% suitable_diseases,
         state != "AUS",
         period == "baseline") %>%
  # Exclude NSW for varicella (not notifiable)
  filter(!(disease %in% varicella_diseases & state == "NSW")) %>%
  group_by(disease, state) %>%
  summarise(
    baseline_n_years     = n(),
    baseline_mean_rate   = mean(rate_per_100k, na.rm = TRUE),
    baseline_sd_rate     = sd(rate_per_100k, na.rm = TRUE),
    baseline_cv_rate     = if_else(baseline_mean_rate > 0,
                                   baseline_sd_rate / baseline_mean_rate,
                                   NA_real_),
    baseline_median_rate = median(rate_per_100k, na.rm = TRUE),
    baseline_min_rate    = min(rate_per_100k, na.rm = TRUE),
    baseline_max_rate    = max(rate_per_100k, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    t_crit = qt(0.975, baseline_n_years - 1),
    pi_lower_rate = pmax(0, baseline_mean_rate -
                           t_crit * baseline_sd_rate *
                           sqrt(1 + 1 / baseline_n_years)),
    pi_upper_rate = baseline_mean_rate +
                      t_crit * baseline_sd_rate *
                      sqrt(1 + 1 / baseline_n_years)
  ) %>%
  select(-t_crit)

cat("State-level baselines:", nrow(baseline_state), "rows\n")
cat("States:", paste(sort(unique(baseline_state$state)), collapse = ", "), "\n")

# =============================================================================
# 5. Join baselines to full time series → O/E ratios
# =============================================================================
cat("\n--- Building annual_with_baselines ---\n")

# National level
annual_national <- annual %>%
  filter(disease %in% suitable_diseases, state == "AUS") %>%
  left_join(
    baseline_national %>%
      select(disease, baseline_mean_rate, baseline_sd_rate, baseline_cv_rate,
             pi_lower_rate, pi_upper_rate,
             baseline_mean_count, pi_lower_count, pi_upper_count,
             count_tier, trend_significant, trend_slope, trend_pct_per_year),
    by = "disease"
  ) %>%
  mutate(
    oe_ratio = rate_per_100k / baseline_mean_rate,
    pct_change = (rate_per_100k - baseline_mean_rate) / baseline_mean_rate * 100,
    below_pi = rate_per_100k < pi_lower_rate,
    above_pi = rate_per_100k > pi_upper_rate,
    oe_ratio_count = count / baseline_mean_count,
    pct_change_count = (count - baseline_mean_count) / baseline_mean_count * 100,
    below_pi_count = count < pi_lower_count,
    above_pi_count = count > pi_upper_count,
    level = "national"
  )

# State level
annual_state <- annual %>%
  filter(disease %in% suitable_diseases, state != "AUS") %>%
  # Exclude NSW for varicella
  filter(!(disease %in% varicella_diseases & state == "NSW")) %>%
  left_join(
    baseline_state %>%
      select(disease, state, baseline_mean_rate, baseline_sd_rate,
             baseline_cv_rate, pi_lower_rate, pi_upper_rate),
    by = c("disease", "state")
  ) %>%
  mutate(
    oe_ratio = rate_per_100k / baseline_mean_rate,
    pct_change = (rate_per_100k - baseline_mean_rate) / baseline_mean_rate * 100,
    below_pi = rate_per_100k < pi_lower_rate,
    above_pi = rate_per_100k > pi_upper_rate,
    # No state-level counts available — leave as NA
    oe_ratio_count = NA_real_,
    pct_change_count = NA_real_,
    below_pi_count = NA,
    above_pi_count = NA,
    level = "state",
    # Carry national count_tier and trend info through for reference
    count_tier = NA_character_,
    trend_significant = NA,
    trend_slope = NA_real_,
    trend_pct_per_year = NA_real_,
    baseline_mean_count = NA_real_,
    pi_lower_count = NA_real_,
    pi_upper_count = NA_real_
  )

annual_with_baselines <- bind_rows(annual_national, annual_state)

# Handle Inf/NaN O/E ratios (division by zero for ultra-low diseases)
annual_with_baselines <- annual_with_baselines %>%
  mutate(
    across(c(oe_ratio, pct_change, oe_ratio_count, pct_change_count),
           ~ if_else(is.infinite(.) | is.nan(.), NA_real_, .))
  )

# =============================================================================
# 6. Save outputs
# =============================================================================
write_csv(baseline_national, file.path(outdir, "baselines_annual_national.csv"))
cat("Saved: baselines_annual_national.csv (", nrow(baseline_national), "rows )\n")

write_csv(baseline_state, file.path(outdir, "baselines_annual_by_state.csv"))
cat("Saved: baselines_annual_by_state.csv (", nrow(baseline_state), "rows )\n")

write_csv(annual_with_baselines, file.path(outdir, "annual_with_baselines.csv"))
cat("Saved: annual_with_baselines.csv (", nrow(annual_with_baselines), "rows )\n")

# =============================================================================
# 7. Quick validation
# =============================================================================
cat("\n=== Quick Validation ===\n")

# Influenza 2020 O/E should be ~0.07 (93% reduction)
flu_2020 <- annual_with_baselines %>%
  filter(disease == "Influenza (laboratory confirmed)",
         year == 2020, level == "national")
cat("Influenza 2020 O/E ratio:", round(flu_2020$oe_ratio, 3),
    "(expected ~0.07-0.13)\n")
cat("Influenza 2020 pct_change:", round(flu_2020$pct_change, 1), "%\n")

# Ross River 2020 should be ~2.0 (increase)
rr_2020 <- annual_with_baselines %>%
  filter(disease == "Ross River virus infection",
         year == 2020, level == "national")
cat("Ross River 2020 O/E ratio:", round(rr_2020$oe_ratio, 3),
    "(expected ~1-3, epidemic year)\n")

# Count tier distribution
cat("\nCount tier distribution:\n")
baseline_national %>% count(count_tier) %>% print()

# High CV diseases
cat("\nHigh variability diseases (CV > 0.5):\n")
baseline_national %>%
  filter(baseline_cv_rate > 0.5) %>%
  select(disease, baseline_mean_rate, baseline_cv_rate, count_tier) %>%
  arrange(desc(baseline_cv_rate)) %>%
  print(n = 20)

cat("\n=== Tier 2 baselines complete ===\n")
