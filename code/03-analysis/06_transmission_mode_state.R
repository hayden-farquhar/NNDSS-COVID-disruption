#!/usr/bin/env Rscript
# =============================================================================
# 06_transmission_mode_state.R
#
# Transmission mode comparison, state-level heterogeneity analysis,
# and border-sensitive disease analysis.
#
# Inputs:
#   data/processed/annual_with_baselines.csv
#   data/processed/disruption_metrics_national.csv
#   data/processed/baselines_annual_by_state.csv
#
# Outputs:
#   data/processed/disruption_by_transmission_mode.csv
#   data/processed/disruption_by_state.csv
#   data/processed/border_sensitive_analysis.csv
#   output/tables/transmission_mode_summary.csv
#   output/tables/state_heterogeneity_table.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(broom)
})

outdir <- "data/processed"
tabdir <- "output/tables"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Transmission Mode & State Analysis ===\n\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
abl       <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                      show_col_types = FALSE)
metrics   <- read_csv(file.path(outdir, "disruption_metrics_national.csv"),
                      show_col_types = FALSE)

# =============================================================================
# 1. Transmission mode comparison
# =============================================================================
cat("--- Transmission mode comparison ---\n")

# Exclude ultra-low count diseases (unreliable O/E)
metrics_reliable <- metrics %>%
  filter(count_tier != "ultra_low")

tm_summary <- metrics_reliable %>%
  group_by(transmission_mode) %>%
  summarise(
    n_diseases          = n(),
    # COVID-acute disruption
    mean_covid_oe       = mean(covid_oe_rate, na.rm = TRUE),
    median_covid_oe     = median(covid_oe_rate, na.rm = TRUE),
    sd_covid_oe         = sd(covid_oe_rate, na.rm = TRUE),
    min_covid_oe        = min(covid_oe_rate, na.rm = TRUE),
    max_covid_oe        = max(covid_oe_rate, na.rm = TRUE),
    mean_covid_pct      = mean(covid_pct_change, na.rm = TRUE),
    n_decreased         = sum(covid_direction == "decreased", na.rm = TRUE),
    n_increased         = sum(covid_direction == "increased", na.rm = TRUE),
    n_unchanged         = sum(covid_direction == "unchanged", na.rm = TRUE),
    # Post-COVID recovery
    mean_post_oe        = mean(post_oe_rate, na.rm = TRUE),
    median_post_oe      = median(post_oe_rate, na.rm = TRUE),
    sd_post_oe          = sd(post_oe_rate, na.rm = TRUE),
    # Trajectory counts
    n_overshoot         = sum(trajectory == "overshoot", na.rm = TRUE),
    n_returned          = sum(trajectory == "returned", na.rm = TRUE),
    n_partial           = sum(trajectory == "partial_recovery", na.rm = TRUE),
    n_suppressed        = sum(trajectory == "sustained_suppression", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_covid_oe)

cat("\nTransmission mode disruption summary:\n")
tm_summary %>%
  select(transmission_mode, n_diseases, mean_covid_oe, median_covid_oe,
         sd_covid_oe, mean_post_oe) %>%
  mutate(across(where(is.numeric) & !c(n_diseases), ~ round(., 3))) %>%
  print()

# --- Kruskal-Wallis test: does transmission mode predict COVID disruption? ---
kw_test <- kruskal.test(covid_oe_rate ~ transmission_mode,
                        data = metrics_reliable)
cat("\nKruskal-Wallis test (COVID O/E ~ transmission mode):\n")
cat("  chi-squared =", round(kw_test$statistic, 2),
    ", df =", kw_test$parameter,
    ", p =", format.pval(kw_test$p.value, digits = 3), "\n")

# --- Pairwise Wilcoxon between transmission modes ---
if (kw_test$p.value < 0.05) {
  pw <- pairwise.wilcox.test(metrics_reliable$covid_oe_rate,
                             metrics_reliable$transmission_mode,
                             p.adjust.method = "BH",
                             exact = FALSE)
  cat("\nPairwise Wilcoxon (BH-adjusted):\n")
  print(round(pw$p.value, 3))
}

# --- Linear model: disruption ~ transmission mode + covariates ---
cat("\n--- Transmission mode regression ---\n")

lm_fit <- lm(covid_oe_rate ~ transmission_mode + log1p(baseline_mean_rate) +
               border_sensitive + vaccine_preventable,
             data = metrics_reliable)

cat("\nRegression: COVID O/E ~ transmission_mode + log(baseline_rate) + border + vaccine\n")
tidy(lm_fit) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  print(n = 15)
cat("  R-squared:", round(summary(lm_fit)$r.squared, 3),
    "  Adj R-squared:", round(summary(lm_fit)$adj.r.squared, 3), "\n")

# Also test post-COVID recovery
kw_post <- kruskal.test(post_oe_rate ~ transmission_mode,
                        data = metrics_reliable)
cat("\nKruskal-Wallis test (post-COVID O/E ~ transmission mode):\n")
cat("  chi-squared =", round(kw_post$statistic, 2),
    ", df =", kw_post$parameter,
    ", p =", format.pval(kw_post$p.value, digits = 3), "\n")

# =============================================================================
# 2. State-level heterogeneity
# =============================================================================
cat("\n--- State-level heterogeneity ---\n")

state_data <- abl %>%
  filter(level == "state", !is.na(baseline_mean_rate), baseline_mean_rate > 0)

# Compute state-period disruption metrics
state_disruption <- state_data %>%
  group_by(disease, state, transmission_mode, border_sensitive) %>%
  summarise(
    baseline_mean_rate = mean(rate_per_100k[period == "baseline"], na.rm = TRUE),
    covid_mean_rate    = mean(rate_per_100k[period == "covid_acute"], na.rm = TRUE),
    post_mean_rate     = mean(rate_per_100k[period == "post_covid"], na.rm = TRUE),
    covid_oe = covid_mean_rate / baseline_mean_rate,
    post_oe  = post_mean_rate / baseline_mean_rate,
    covid_pct = (covid_mean_rate - baseline_mean_rate) / baseline_mean_rate * 100,
    .groups = "drop"
  ) %>%
  mutate(across(c(covid_oe, post_oe, covid_pct),
                ~ if_else(is.infinite(.) | is.nan(.), NA_real_, .)))

cat("State disruption rows:", nrow(state_disruption), "\n")

# State-level summary across all diseases
state_summary <- state_disruption %>%
  group_by(state) %>%
  summarise(
    n_diseases      = n(),
    mean_covid_oe   = mean(covid_oe, na.rm = TRUE),
    median_covid_oe = median(covid_oe, na.rm = TRUE),
    mean_post_oe    = mean(post_oe, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_covid_oe)

cat("\nState-level mean disruption (all diseases):\n")
state_summary %>%
  mutate(across(where(is.numeric) & !c(n_diseases), ~ round(., 3))) %>%
  print()

# State disruption by transmission mode
state_tm <- state_disruption %>%
  group_by(state, transmission_mode) %>%
  summarise(
    n = n(),
    mean_covid_oe   = mean(covid_oe, na.rm = TRUE),
    median_covid_oe = median(covid_oe, na.rm = TRUE),
    mean_post_oe    = mean(post_oe, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nState × transmission mode disruption (respiratory only):\n")
state_tm %>%
  filter(transmission_mode == "respiratory") %>%
  mutate(across(where(is.numeric) & !c(n), ~ round(., 3))) %>%
  arrange(mean_covid_oe) %>%
  print()

# --- Test: VIC vs WA respiratory disruption ---
cat("\n--- VIC vs WA respiratory disease comparison ---\n")
vic_resp <- state_disruption %>%
  filter(state == "VIC", transmission_mode == "respiratory") %>%
  pull(covid_oe)
wa_resp <- state_disruption %>%
  filter(state == "WA", transmission_mode == "respiratory") %>%
  pull(covid_oe)

if (length(vic_resp) >= 3 & length(wa_resp) >= 3) {
  vw_test <- wilcox.test(vic_resp, wa_resp, exact = FALSE)
  cat("VIC median COVID O/E:", round(median(vic_resp, na.rm = TRUE), 3), "\n")
  cat("WA  median COVID O/E:", round(median(wa_resp, na.rm = TRUE), 3), "\n")
  cat("Wilcoxon p =", format.pval(vw_test$p.value, digits = 3), "\n")
} else {
  cat("Insufficient data for VIC vs WA comparison\n")
}

# State-level COVID O/E for key diseases (table for manuscript)
key_diseases <- c("Influenza (laboratory confirmed)", "Pertussis",
                  "Salmonellosis", "Campylobacteriosis",
                  "Gonococcal infection", "Ross River virus infection",
                  "Dengue virus infection", "Tuberculosis")

state_table <- state_disruption %>%
  filter(disease %in% key_diseases) %>%
  select(disease, state, covid_oe, post_oe) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# =============================================================================
# 3. Border-sensitive disease analysis
# =============================================================================
cat("\n--- Border-sensitive disease analysis ---\n")

border_diseases <- metrics %>%
  filter(border_sensitive == TRUE) %>%
  pull(disease)

cat("Border-sensitive diseases:", length(border_diseases), "\n")
cat(paste(border_diseases, collapse = "\n"), "\n\n")

border_analysis <- abl %>%
  filter(level == "national", disease %in% border_diseases) %>%
  select(disease, year, period, rate_per_100k, count,
         baseline_mean_rate, oe_ratio, pct_change,
         transmission_mode) %>%
  arrange(disease, year)

# Summarise by period
border_summary <- border_analysis %>%
  group_by(disease, transmission_mode) %>%
  summarise(
    baseline_mean = mean(rate_per_100k[period == "baseline"], na.rm = TRUE),
    # 2020: borders largely closed from March
    rate_2020     = rate_per_100k[year == 2020],
    oe_2020       = oe_ratio[year == 2020],
    # 2021: borders still closed
    rate_2021     = rate_per_100k[year == 2021],
    oe_2021       = oe_ratio[year == 2021],
    # 2022: borders reopened Feb-March
    rate_2022     = rate_per_100k[year == 2022],
    oe_2022       = oe_ratio[year == 2022],
    # Most recent
    rate_2024     = rate_per_100k[year == 2024],
    oe_2024       = oe_ratio[year == 2024],
    # Mean COVID O/E
    covid_oe      = mean(oe_ratio[period == "covid_acute"], na.rm = TRUE),
    # Recovery status
    post_oe       = mean(oe_ratio[period == "post_covid"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    collapse_pct_2020 = (1 - oe_2020) * 100,
    recovery_status = case_when(
      post_oe > 0.9 ~ "recovered",
      post_oe > 0.5 ~ "partial",
      TRUE           ~ "suppressed"
    )
  ) %>%
  arrange(covid_oe)

cat("Border-sensitive disease disruption:\n")
border_summary %>%
  select(disease, baseline_mean, oe_2020, oe_2021, oe_2022, oe_2024,
         covid_oe, recovery_status) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print(width = 120)

# Compare border-sensitive vs non-border diseases
cat("\nBorder-sensitive vs non-border disruption:\n")
metrics_reliable %>%
  group_by(border_sensitive) %>%
  summarise(
    n = n(),
    mean_covid_oe = mean(covid_oe_rate, na.rm = TRUE),
    sd_covid_oe   = sd(covid_oe_rate, na.rm = TRUE),
    mean_post_oe  = mean(post_oe_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric) & !c(n), ~ round(., 3))) %>%
  print()

border_test <- wilcox.test(
  metrics_reliable$covid_oe_rate[metrics_reliable$border_sensitive == TRUE],
  metrics_reliable$covid_oe_rate[metrics_reliable$border_sensitive == FALSE],
  exact = FALSE
)
cat("Wilcoxon border vs non-border p =",
    format.pval(border_test$p.value, digits = 3), "\n")

# =============================================================================
# 4. Vaccine-preventable disease analysis
# =============================================================================
cat("\n--- Vaccine-preventable disease analysis ---\n")

metrics_reliable %>%
  group_by(vaccine_preventable) %>%
  summarise(
    n = n(),
    mean_covid_oe = mean(covid_oe_rate, na.rm = TRUE),
    sd_covid_oe   = sd(covid_oe_rate, na.rm = TRUE),
    mean_post_oe  = mean(post_oe_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric) & !c(n), ~ round(., 3))) %>%
  print()

vax_test <- wilcox.test(
  metrics_reliable$covid_oe_rate[metrics_reliable$vaccine_preventable == TRUE],
  metrics_reliable$covid_oe_rate[metrics_reliable$vaccine_preventable == FALSE],
  exact = FALSE
)
cat("Wilcoxon VPD vs non-VPD p =",
    format.pval(vax_test$p.value, digits = 3), "\n")

# =============================================================================
# 5. Save outputs
# =============================================================================
write_csv(tm_summary, file.path(tabdir, "transmission_mode_summary.csv"))
cat("\nSaved: transmission_mode_summary.csv\n")

write_csv(state_disruption, file.path(outdir, "disruption_by_state.csv"))
cat("Saved: disruption_by_state.csv (", nrow(state_disruption), "rows )\n")

write_csv(state_table, file.path(tabdir, "state_heterogeneity_table.csv"))
cat("Saved: state_heterogeneity_table.csv\n")

write_csv(border_summary, file.path(outdir, "border_sensitive_analysis.csv"))
cat("Saved: border_sensitive_analysis.csv\n")

cat("\n=== Transmission mode & state analysis complete ===\n")
