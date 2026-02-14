#!/usr/bin/env Rscript
# =============================================================================
# 07_recovery_trajectories.R
#
# Detailed post-pandemic recovery trajectory analysis:
#   - Year-by-year recovery tracking with trajectory refinement
#   - Time-to-recovery estimates (year first exceeded 90% of baseline)
#   - Trajectory patterns by transmission mode, VPD, border sensitivity
#   - Annual trajectory evolution (how classification changes over time)
#   - Monthly recovery dynamics for Tier 1 diseases
#
# Inputs:
#   data/processed/annual_with_baselines.csv
#   data/processed/disruption_metrics_national.csv
#   data/processed/tier1_arima_forecasts.csv
#
# Outputs:
#   data/processed/recovery_trajectories.csv       — 47 diseases, detailed
#   data/processed/recovery_year_by_year.csv       — annual trajectory tracking
#   data/processed/recovery_monthly_tier1.csv      — monthly for 3 diseases
#   output/tables/recovery_summary_table.csv       — manuscript table
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

cat("=== Phase 4: Recovery Trajectory Analysis ===\n\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
abl     <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                    show_col_types = FALSE)
metrics <- read_csv(file.path(outdir, "disruption_metrics_national.csv"),
                    show_col_types = FALSE)
arima   <- read_csv(file.path(outdir, "tier1_arima_forecasts.csv"),
                    show_col_types = FALSE)

national <- abl %>% filter(level == "national")

# =============================================================================
# 1. Year-by-year recovery tracking
# =============================================================================
cat("--- Year-by-year recovery tracking ---\n")

recovery_yby <- national %>%
  filter(period %in% c("covid_acute", "transition", "post_covid")) %>%
  select(disease, year, period, rate_per_100k, count,
         transmission_mode, border_sensitive, vaccine_preventable,
         count_tier, baseline_mean_rate, oe_ratio, pct_change) %>%
  group_by(disease) %>%
  mutate(
    # Year-specific trajectory classification
    year_trajectory = case_when(
      count_tier == "ultra_low" ~ "insufficient_data",
      oe_ratio > 1.2  ~ "overshoot",
      oe_ratio >= 0.9 ~ "returned",
      oe_ratio >= 0.5 ~ "partial_recovery",
      TRUE            ~ "sustained_suppression"
    ),
    # Was the disease suppressed during COVID acute?
    was_suppressed = any(oe_ratio[period == "covid_acute"] < 0.9, na.rm = TRUE),
    # Track recovery from nadir
    nadir_oe = min(oe_ratio[period == "covid_acute"], na.rm = TRUE),
    nadir_year = year[period == "covid_acute"][which.min(oe_ratio[period == "covid_acute"])],
    # Percent recovery from nadir (0% = still at nadir, 100% = back to baseline)
    recovery_from_nadir = if_else(
      nadir_oe < 0.9,
      pmin(100, (oe_ratio - nadir_oe) / (1 - nadir_oe) * 100),
      NA_real_
    )
  ) %>%
  ungroup()

cat("Recovery year-by-year rows:", nrow(recovery_yby), "\n")

# =============================================================================
# 2. Time-to-recovery estimates
# =============================================================================
cat("\n--- Time-to-recovery estimates ---\n")

# For diseases that were suppressed, find first year O/E >= 0.9
# Use a simple split-apply approach to avoid summarise column issues
suppressed_diseases <- recovery_yby %>%
  filter(was_suppressed, count_tier != "ultra_low") %>%
  distinct(disease) %>%
  pull(disease)

ttr_list <- lapply(suppressed_diseases, function(dis) {
  d <- recovery_yby %>% filter(disease == dis) %>% arrange(year)
  covid_rows <- d %>% filter(period == "covid_acute")
  nadir_idx  <- which.min(covid_rows$oe_ratio)
  n_year     <- covid_rows$year[nadir_idx]
  n_oe       <- covid_rows$oe_ratio[nadir_idx]

  post_nadir <- d %>% filter(year > n_year)
  rec_rows   <- post_nadir %>% filter(oe_ratio >= 0.9)
  full_rows  <- post_nadir %>% filter(oe_ratio >= 1.0)

  tibble(
    disease            = dis,
    nadir_year         = n_year,
    nadir_oe           = n_oe,
    recovery_year      = if (nrow(rec_rows) > 0) min(rec_rows$year) else NA_integer_,
    full_recovery_year = if (nrow(full_rows) > 0) min(full_rows$year) else NA_integer_,
    latest_year        = max(d$year),
    latest_oe          = d$oe_ratio[d$year == max(d$year)]
  )
})

time_to_recovery <- bind_rows(ttr_list) %>%
  mutate(
    years_to_recovery = recovery_year - nadir_year,
    years_to_full     = full_recovery_year - nadir_year,
    recovered_by_2024 = !is.na(recovery_year) & recovery_year <= 2024,
    still_suppressed  = is.na(recovery_year)
  )

cat("Diseases that were suppressed during COVID:", nrow(time_to_recovery), "\n")
cat("Recovered (O/E >= 0.9) by 2024:",
    sum(time_to_recovery$recovered_by_2024, na.rm = TRUE), "\n")
cat("Still suppressed (O/E < 0.9):",
    sum(time_to_recovery$still_suppressed), "\n")

cat("\nTime-to-recovery summary:\n")
time_to_recovery %>%
  filter(!is.na(years_to_recovery)) %>%
  summarise(
    n = n(),
    mean_years = mean(years_to_recovery, na.rm = TRUE),
    median_years = median(years_to_recovery, na.rm = TRUE),
    min_years = min(years_to_recovery, na.rm = TRUE),
    max_years = max(years_to_recovery, na.rm = TRUE)
  ) %>%
  print()

cat("\nStill suppressed diseases:\n")
time_to_recovery %>%
  filter(still_suppressed) %>%
  select(disease, nadir_year, nadir_oe, latest_oe) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print(n = 20)

# =============================================================================
# 3. Detailed trajectory classification with refinements
# =============================================================================
cat("\n--- Refined trajectory classification ---\n")

trajectories <- metrics %>%
  select(disease, transmission_mode, border_sensitive, vaccine_preventable,
         count_tier, covid_oe_rate, covid_pct_change,
         post_oe_rate, post_pct_change, trajectory,
         cumulative_deficit_total, cumulative_excess_total, net_cumulative) %>%
  left_join(
    time_to_recovery %>%
      select(disease, nadir_year, nadir_oe, recovery_year,
             years_to_recovery, recovered_by_2024, still_suppressed),
    by = "disease"
  ) %>%
  mutate(
    # Refined trajectory with more nuance
    trajectory_refined = case_when(
      count_tier == "ultra_low" ~ "insufficient_data",
      # Overshoot subcategories
      post_oe_rate > 2.0 ~ "major_overshoot",
      post_oe_rate > 1.2 ~ "moderate_overshoot",
      # Returned
      post_oe_rate >= 0.9 & post_oe_rate <= 1.2 ~ "returned",
      # Partial recovery subcategories
      post_oe_rate >= 0.7 & covid_oe_rate < 0.9 ~ "partial_nearing",
      post_oe_rate >= 0.5 & covid_oe_rate < 0.9 ~ "partial_slow",
      post_oe_rate >= 0.5 ~ "partial_recovery",
      # Sustained suppression
      TRUE ~ "sustained_suppression"
    ),
    # Was there a post-COVID peak (any single year > 1.5)?
    had_rebound_peak = disease %in% (
      recovery_yby %>%
        filter(oe_ratio > 1.5, period %in% c("transition", "post_covid")) %>%
        pull(disease) %>% unique()
    ),
    # Immunity debt signal: suppressed during COVID + overshoot post-COVID
    immunity_debt_signal = covid_oe_rate < 0.5 & post_oe_rate > 1.2
  )

cat("Refined trajectory distribution:\n")
trajectories %>%
  count(trajectory_refined) %>%
  arrange(desc(n)) %>%
  print()

cat("\nImmunity debt signal diseases:\n")
trajectories %>%
  filter(immunity_debt_signal) %>%
  select(disease, transmission_mode, covid_oe_rate, post_oe_rate,
         cumulative_deficit_total, cumulative_excess_total) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  print(n = 20)

# =============================================================================
# 4. Trajectory patterns by disease characteristics
# =============================================================================
cat("\n--- Trajectory patterns by group ---\n")

# By transmission mode
cat("By transmission mode:\n")
trajectories %>%
  filter(count_tier != "ultra_low") %>%
  count(transmission_mode, trajectory) %>%
  pivot_wider(names_from = trajectory, values_from = n, values_fill = 0) %>%
  print()

# By VPD status
cat("\nBy vaccine-preventable status:\n")
trajectories %>%
  filter(count_tier != "ultra_low") %>%
  count(vaccine_preventable, trajectory) %>%
  pivot_wider(names_from = trajectory, values_from = n, values_fill = 0) %>%
  print()

# By border sensitivity
cat("\nBy border sensitivity:\n")
trajectories %>%
  filter(count_tier != "ultra_low") %>%
  count(border_sensitive, trajectory) %>%
  pivot_wider(names_from = trajectory, values_from = n, values_fill = 0) %>%
  print()

# =============================================================================
# 5. Monthly recovery dynamics for Tier 1 diseases
# =============================================================================
cat("\n--- Monthly recovery dynamics (Tier 1) ---\n")

monthly_recovery <- arima %>%
  filter(!is.na(observed), period != "baseline") %>%
  mutate(
    date = as.Date(paste(year, month, 1, sep = "-")),
    month_trajectory = case_when(
      oe_ratio > 1.2 ~ "overshoot",
      oe_ratio >= 0.9 ~ "returned",
      oe_ratio >= 0.5 ~ "partial",
      TRUE ~ "suppressed"
    )
  ) %>%
  group_by(disease) %>%
  mutate(
    # Cumulative deficit and excess from COVID onset
    cum_deficit = cumsum(pmax(0, expected - observed)),
    cum_excess  = cumsum(pmax(0, observed - expected)),
    cum_net     = cum_excess - cum_deficit,
    # Rolling 3-month average O/E
    oe_roll3 = zoo::rollmean(oe_ratio, k = 3, fill = NA, align = "right")
  ) %>%
  ungroup()

# First month of sustained recovery (3-month rolling O/E >= 0.9)
monthly_first_recovery <- monthly_recovery %>%
  filter(!is.na(oe_roll3), oe_roll3 >= 0.9) %>%
  group_by(disease) %>%
  slice_min(date, n = 1) %>%
  ungroup() %>%
  select(disease, first_recovery_date = date, first_recovery_oe = oe_roll3)

cat("First month of sustained recovery (rolling 3m O/E >= 0.9):\n")
monthly_first_recovery %>% print()

# =============================================================================
# 6. Save outputs
# =============================================================================
write_csv(trajectories, file.path(outdir, "recovery_trajectories.csv"))
cat("\nSaved: recovery_trajectories.csv (", nrow(trajectories), "rows )\n")

write_csv(recovery_yby, file.path(outdir, "recovery_year_by_year.csv"))
cat("Saved: recovery_year_by_year.csv (", nrow(recovery_yby), "rows )\n")

write_csv(monthly_recovery, file.path(outdir, "recovery_monthly_tier1.csv"))
cat("Saved: recovery_monthly_tier1.csv (", nrow(monthly_recovery), "rows )\n")

# Manuscript table
recovery_table <- trajectories %>%
  select(disease, transmission_mode, count_tier,
         covid_oe_rate, post_oe_rate,
         trajectory, trajectory_refined,
         nadir_year, nadir_oe, recovery_year, years_to_recovery,
         cumulative_deficit_total, cumulative_excess_total, net_cumulative,
         immunity_debt_signal) %>%
  arrange(transmission_mode, trajectory_refined)

write_csv(recovery_table, file.path(tabdir, "recovery_summary_table.csv"))
cat("Saved: recovery_summary_table.csv\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n##############################################################\n")
cat("#          RECOVERY TRAJECTORY SUMMARY                       #\n")
cat("##############################################################\n\n")

cat("Trajectory distribution (refined):\n")
trajectories %>%
  filter(count_tier != "ultra_low") %>%
  count(trajectory_refined) %>%
  arrange(desc(n)) %>%
  print()

cat("\nDiseases showing immunity debt signal:", sum(trajectories$immunity_debt_signal), "\n")

cat("\nMost complete recoveries (border-sensitive):\n")
trajectories %>%
  filter(border_sensitive, !is.na(recovery_year)) %>%
  select(disease, covid_oe_rate, recovery_year, post_oe_rate) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print()

cat("\n=== Recovery trajectory analysis complete ===\n")
