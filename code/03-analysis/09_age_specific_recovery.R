#!/usr/bin/env Rscript
# =============================================================================
# 09_age_specific_recovery.R
#
# Age-stratified disruption and recovery for 4 diseases with record-level data:
#   - Influenza, Meningococcal, Pneumococcal, Salmonella
#
# Tests hypothesis: children show greater immunity debt rebound than adults.
#
# Inputs:
#   data/processed/tier1_annual_by_age.csv
#
# Outputs:
#   data/processed/age_specific_disruption.csv   — age × disease × year
#   output/tables/age_disruption_table.csv       — summary for manuscript
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

outdir <- "data/processed"
tabdir <- "output/tables"

cat("=== Age-Specific Recovery Analysis ===\n\n")

# -----------------------------------------------------------------------------
# Load and harmonise age data
# -----------------------------------------------------------------------------
age_data <- read_csv(file.path(outdir, "tier1_annual_by_age.csv"),
                     show_col_types = FALSE)

cat("Raw data:\n")
cat("  Rows:", nrow(age_data), "\n")
cat("  Diseases:", paste(unique(age_data$disease), collapse = "; "), "\n")
cat("  Age groups:", paste(sort(unique(age_data$age_group)), collapse = ", "), "\n\n")

# Harmonise age groups into broader categories
# Pneumococcal uses "< 2 yo" and "02-04"; others use "00-04"
age_data <- age_data %>%
  mutate(
    age_broad = case_when(
      age_group %in% c("00-04", "< 2 yo", "02-04") ~ "0-4",
      age_group %in% c("05-09", "10-14")           ~ "5-14",
      age_group %in% c("15-19", "20-24")           ~ "15-24",
      age_group %in% c("25-29", "30-34", "35-39",
                        "40-44", "45-49")           ~ "25-49",
      age_group %in% c("50-54", "55-59", "60-64")  ~ "50-64",
      age_group %in% c("65-69", "70-74", "75-79",
                        "80-84", "85+")             ~ "65+",
      age_group == "Unknown"                        ~ "Unknown",
      TRUE                                          ~ "Other"
    ),
    period = case_when(
      year >= 2015 & year <= 2019 ~ "baseline",
      year >= 2020 & year <= 2021 ~ "covid_acute",
      year == 2022                ~ "transition",
      year >= 2023                ~ "post_covid",
      TRUE                        ~ "pre_baseline"
    )
  )

# Aggregate to broad age groups
age_broad <- age_data %>%
  filter(age_broad != "Unknown", year >= 2015) %>%
  group_by(disease, year, period, age_broad) %>%
  summarise(count = sum(count), .groups = "drop")

cat("Broad age groups:", paste(sort(unique(age_broad$age_broad)), collapse = ", "), "\n")

# =============================================================================
# 1. Baseline by age group (2015-2019)
# =============================================================================
cat("\n--- Computing age-specific baselines ---\n")

age_baselines <- age_broad %>%
  filter(period == "baseline") %>%
  group_by(disease, age_broad) %>%
  summarise(
    baseline_mean = mean(count, na.rm = TRUE),
    baseline_sd   = sd(count, na.rm = TRUE),
    .groups = "drop"
  )

# Join baselines and compute O/E
age_disruption <- age_broad %>%
  left_join(age_baselines, by = c("disease", "age_broad")) %>%
  mutate(
    oe_ratio   = count / baseline_mean,
    pct_change = (count - baseline_mean) / baseline_mean * 100
  ) %>%
  mutate(across(c(oe_ratio, pct_change),
                ~ if_else(is.infinite(.) | is.nan(.), NA_real_, .)))

# =============================================================================
# 2. Age-specific period summaries
# =============================================================================
cat("--- Age-specific period summaries ---\n\n")

age_period <- age_disruption %>%
  filter(period != "baseline") %>%
  group_by(disease, age_broad, period) %>%
  summarise(
    mean_count   = mean(count, na.rm = TRUE),
    mean_oe      = mean(oe_ratio, na.rm = TRUE),
    mean_pct     = mean(pct_change, na.rm = TRUE),
    .groups = "drop"
  )

# Children (0-4) vs adults (25-49) for each disease
cat("Children (0-4) vs Adults (25-49) COVID-acute disruption:\n")
age_period %>%
  filter(period == "covid_acute",
         age_broad %in% c("0-4", "25-49")) %>%
  select(disease, age_broad, mean_oe) %>%
  pivot_wider(names_from = age_broad, values_from = mean_oe,
              names_prefix = "oe_") %>%
  mutate(across(starts_with("oe_"), ~ round(., 3))) %>%
  print()

cat("\nChildren (0-4) vs Adults (25-49) post-COVID recovery:\n")
age_period %>%
  filter(period == "post_covid",
         age_broad %in% c("0-4", "25-49")) %>%
  select(disease, age_broad, mean_oe) %>%
  pivot_wider(names_from = age_broad, values_from = mean_oe,
              names_prefix = "oe_") %>%
  mutate(across(starts_with("oe_"), ~ round(., 3))) %>%
  print()

# =============================================================================
# 3. Full age × period profile per disease
# =============================================================================
cat("\n--- Full age profiles ---\n")

for (dis in unique(age_disruption$disease)) {
  cat("\n  ", dis, ":\n")
  age_disruption %>%
    filter(disease == dis, period %in% c("covid_acute", "post_covid")) %>%
    group_by(age_broad, period) %>%
    summarise(mean_oe = mean(oe_ratio, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = mean_oe) %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    arrange(age_broad) %>%
    print()
}

# =============================================================================
# 4. Statistical test: do children show greater rebound?
# =============================================================================
cat("\n--- Testing age-differential immunity debt ---\n")

# For each disease, compare post-COVID O/E between children (0-4) and adults
age_tests <- age_disruption %>%
  filter(period == "post_covid",
         age_broad %in% c("0-4", "25-49", "65+")) %>%
  select(disease, year, age_broad, oe_ratio) %>%
  pivot_wider(names_from = age_broad, values_from = oe_ratio)

for (dis in unique(age_tests$disease)) {
  d <- age_tests %>% filter(disease == dis)
  cat("\n  ", dis, ":\n")
  cat("    0-4 mean O/E:", round(mean(d$`0-4`, na.rm = TRUE), 3), "\n")
  cat("    25-49 mean O/E:", round(mean(d$`25-49`, na.rm = TRUE), 3), "\n")
  cat("    65+ mean O/E:", round(mean(d$`65+`, na.rm = TRUE), 3), "\n")
  # Can't do formal test with n=2-3 years, just report descriptively
}

# =============================================================================
# 5. Save outputs
# =============================================================================
write_csv(age_disruption, file.path(outdir, "age_specific_disruption.csv"))
cat("\nSaved: age_specific_disruption.csv (", nrow(age_disruption), "rows )\n")

# Manuscript table: age × period O/E for each disease
age_table <- age_period %>%
  select(disease, age_broad, period, mean_oe) %>%
  mutate(mean_oe = round(mean_oe, 3)) %>%
  pivot_wider(names_from = period, values_from = mean_oe)

write_csv(age_table, file.path(tabdir, "age_disruption_table.csv"))
cat("Saved: age_disruption_table.csv\n")

cat("\n=== Age-specific recovery analysis complete ===\n")
