#!/usr/bin/env Rscript
# =============================================================================
# 08_immunity_debt.R
#
# Immunity debt analysis — quantify cumulative "missing" cases during COVID
# and compare to post-COVID excess. Monthly resolution for Tier 1, annual
# for all 47 diseases.
#
# Candidate diseases for immunity debt: influenza, pertussis, pneumococcal,
# meningococcal, rotavirus, cryptosporidiosis (+ any showing overshoot)
#
# Inputs:
#   data/processed/annual_with_baselines.csv
#   data/processed/disruption_metrics_national.csv
#   data/processed/tier1_arima_forecasts.csv
#   data/processed/recovery_trajectories.csv
#
# Outputs:
#   data/processed/immunity_debt_annual.csv         — all diseases, annual
#   data/processed/immunity_debt_monthly.csv        — Tier 1, monthly detail
#   output/tables/immunity_debt_summary.csv         — manuscript table
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

outdir <- "data/processed"
tabdir <- "output/tables"

cat("=== Immunity Debt Analysis ===\n\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
abl         <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                        show_col_types = FALSE)
metrics     <- read_csv(file.path(outdir, "disruption_metrics_national.csv"),
                        show_col_types = FALSE)
arima       <- read_csv(file.path(outdir, "tier1_arima_forecasts.csv"),
                        show_col_types = FALSE)
trajectories <- read_csv(file.path(outdir, "recovery_trajectories.csv"),
                         show_col_types = FALSE)

national <- abl %>% filter(level == "national")

# =============================================================================
# 1. Annual immunity debt for all diseases (count-based)
# =============================================================================
cat("--- Annual immunity debt (all diseases) ---\n")

annual_debt <- national %>%
  filter(period != "baseline", !is.na(count), !is.na(baseline_mean_count)) %>%
  mutate(
    expected_count = baseline_mean_count,
    deficit = pmax(0, expected_count - count),
    excess  = pmax(0, count - expected_count)
  ) %>%
  select(disease, year, period, transmission_mode,
         count_tier, border_sensitive, vaccine_preventable,
         count, expected_count, deficit, excess) %>%
  arrange(disease, year) %>%
  group_by(disease) %>%
  mutate(
    cum_deficit = cumsum(deficit),
    cum_excess  = cumsum(excess),
    cum_net     = cum_excess - cum_deficit,
    # Debt ratio: cumulative excess / cumulative deficit
    # >1 means more excess than deficit (debt "overpaid")
    # <1 means deficit not yet recovered
    debt_ratio = cum_excess / pmax(1, cum_deficit)
  ) %>%
  ungroup()

# By-disease summary
annual_debt_summary <- annual_debt %>%
  group_by(disease, transmission_mode, count_tier,
           border_sensitive, vaccine_preventable) %>%
  summarise(
    # COVID-acute period deficit
    deficit_covid_acute  = sum(deficit[period == "covid_acute"]),
    excess_covid_acute   = sum(excess[period == "covid_acute"]),
    # Transition + post-COVID
    deficit_post         = sum(deficit[period %in% c("transition", "post_covid")]),
    excess_post          = sum(excess[period %in% c("transition", "post_covid")]),
    # Totals
    total_deficit        = sum(deficit),
    total_excess         = sum(excess),
    net_cumulative       = sum(excess) - sum(deficit),
    # Debt ratio (total)
    debt_ratio           = sum(excess) / pmax(1, sum(deficit)),
    # Peak single-year deficit (most missing cases in one year)
    peak_deficit_year    = year[which.max(deficit)],
    peak_deficit_count   = max(deficit),
    # Peak single-year excess
    peak_excess_year     = if (any(excess > 0)) year[which.max(excess)] else NA_integer_,
    peak_excess_count    = max(excess),
    .groups = "drop"
  ) %>%
  mutate(
    # Debt classification
    debt_status = case_when(
      total_deficit < 10    ~ "minimal_deficit",
      debt_ratio > 1.5      ~ "overpaid",       # Excess far exceeds deficit
      debt_ratio > 1.0      ~ "repaid",          # Roughly balanced
      debt_ratio > 0.5      ~ "partial_repayment",
      debt_ratio > 0.0      ~ "minimal_repayment",
      TRUE                  ~ "no_debt"          # Never had excess
    )
  ) %>%
  arrange(desc(total_deficit))

cat("Annual debt summary:\n")
annual_debt_summary %>%
  filter(total_deficit > 100) %>%
  select(disease, total_deficit, total_excess, net_cumulative,
         debt_ratio, debt_status) %>%
  mutate(across(c(total_deficit, total_excess, net_cumulative), ~ round(.)),
         debt_ratio = round(debt_ratio, 2)) %>%
  print(n = 20)

cat("\nDebt status distribution:\n")
annual_debt_summary %>%
  filter(count_tier != "ultra_low") %>%
  count(debt_status) %>%
  print()

# =============================================================================
# 2. Monthly immunity debt for Tier 1 diseases
# =============================================================================
cat("\n--- Monthly immunity debt (Tier 1) ---\n")

monthly_debt <- arima %>%
  filter(period != "baseline", !is.na(observed)) %>%
  mutate(
    date = as.Date(paste(year, month, 1, sep = "-")),
    deficit = pmax(0, expected - observed),
    excess  = pmax(0, observed - expected)
  ) %>%
  arrange(disease, year, month) %>%
  group_by(disease) %>%
  mutate(
    cum_deficit = cumsum(deficit),
    cum_excess  = cumsum(excess),
    cum_net     = cum_excess - cum_deficit,
    debt_ratio  = cum_excess / pmax(1, cum_deficit)
  ) %>%
  ungroup()

# Monthly summary per disease
cat("Monthly debt by disease:\n")
monthly_debt %>%
  group_by(disease) %>%
  summarise(
    total_deficit   = sum(deficit),
    total_excess    = sum(excess),
    net             = sum(excess) - sum(deficit),
    debt_ratio      = sum(excess) / pmax(1, sum(deficit)),
    # When did debt start being repaid? (first month of excess)
    first_excess_month = {
      exc <- which(excess > 0)
      if (length(exc) > 0) as.character(date[exc[1]]) else NA_character_
    },
    # When did cumulative net turn positive?
    net_positive_month = {
      pos <- which(cum_net > 0)
      if (length(pos) > 0) as.character(date[pos[1]]) else NA_character_
    },
    .groups = "drop"
  ) %>%
  mutate(across(c(total_deficit, total_excess, net), ~ round(.)),
         debt_ratio = round(debt_ratio, 2)) %>%
  print()

# =============================================================================
# 3. Key immunity debt candidate deep dives
# =============================================================================
cat("\n--- Immunity debt candidate diseases ---\n")

# Diseases with clear immunity debt signal
debt_candidates <- trajectories %>%
  filter(immunity_debt_signal | disease %in% c(
    "Influenza (laboratory confirmed)",
    "Pertussis",
    "Pneumococcal disease (invasive)",
    "Meningococcal disease (invasive)",
    "Rotavirus",
    "Cryptosporidiosis"
  )) %>%
  pull(disease) %>%
  unique()

cat("Debt candidate diseases:", length(debt_candidates), "\n")

# Detailed annual profile for these diseases
candidate_profiles <- national %>%
  filter(disease %in% debt_candidates) %>%
  select(disease, year, period, count, baseline_mean_count, oe_ratio) %>%
  mutate(
    expected = baseline_mean_count,
    deficit  = pmax(0, expected - count),
    excess   = pmax(0, count - expected)
  ) %>%
  group_by(disease) %>%
  mutate(
    cum_deficit = cumsum(deficit),
    cum_excess  = cumsum(excess),
    cum_net     = cum_excess - cum_deficit
  ) %>%
  ungroup()

cat("\nCandidate annual profiles:\n")
for (dis in debt_candidates) {
  cat("\n  ", dis, ":\n")
  d <- candidate_profiles %>% filter(disease == dis)
  d %>%
    select(year, period, count, expected, deficit, excess, cum_net) %>%
    mutate(across(where(is.numeric), ~ round(.))) %>%
    print(n = 12)
}

# =============================================================================
# 4. Save outputs
# =============================================================================
write_csv(annual_debt, file.path(outdir, "immunity_debt_annual.csv"))
cat("\nSaved: immunity_debt_annual.csv (", nrow(annual_debt), "rows )\n")

write_csv(monthly_debt, file.path(outdir, "immunity_debt_monthly.csv"))
cat("Saved: immunity_debt_monthly.csv (", nrow(monthly_debt), "rows )\n")

# Manuscript summary table
debt_table <- annual_debt_summary %>%
  filter(disease %in% debt_candidates | total_deficit > 1000) %>%
  select(disease, transmission_mode,
         deficit_covid_acute, excess_post,
         total_deficit, total_excess, net_cumulative,
         debt_ratio, debt_status,
         peak_deficit_year, peak_excess_year) %>%
  arrange(desc(total_deficit))

write_csv(debt_table, file.path(tabdir, "immunity_debt_summary.csv"))
cat("Saved: immunity_debt_summary.csv\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n##############################################################\n")
cat("#             IMMUNITY DEBT SUMMARY                          #\n")
cat("##############################################################\n\n")

cat("Key findings:\n")
flu <- annual_debt_summary %>%
  filter(disease == "Influenza (laboratory confirmed)")
cat("  Influenza: deficit", round(flu$total_deficit), "cases,",
    "excess", round(flu$total_excess), "cases,",
    "ratio", round(flu$debt_ratio, 2), "\n")

per <- annual_debt_summary %>%
  filter(disease == "Pertussis")
cat("  Pertussis: deficit", round(per$total_deficit), "cases,",
    "excess", round(per$total_excess), "cases,",
    "ratio", round(per$debt_ratio, 2), "\n")

cat("\nDebt repayment status (all diseases with deficit > 100):\n")
annual_debt_summary %>%
  filter(total_deficit > 100) %>%
  count(debt_status) %>%
  print()

cat("\n=== Immunity debt analysis complete ===\n")
