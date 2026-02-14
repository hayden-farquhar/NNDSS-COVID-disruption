#!/usr/bin/env Rscript
# =============================================================================
# 01_data_quality_checks.R
#
# Comprehensive data quality checks on parsed NNDSS datasets.
# Generates a quality report and flags issues for manual review.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

cat("Loading datasets...\n")
rates     <- read_csv("data/processed/nndss_rates_by_state_year.csv", show_col_types = FALSE)
counts    <- read_csv("data/processed/nndss_counts_national_annual.csv", show_col_types = FALSE)
caveats   <- read_csv("data/processed/nndss_data_caveats.csv", show_col_types = FALSE)

cat("\n")
cat("##############################################################\n")
cat("#               NNDSS DATA QUALITY REPORT                    #\n")
cat("##############################################################\n\n")

# =============================================================================
# 1. Disease name reconciliation: rates vs counts
# =============================================================================
cat("=== 1. DISEASE NAME RECONCILIATION ===\n\n")

rate_diseases  <- sort(unique(rates$disease))
count_diseases <- sort(unique(counts$disease))

in_rates_only  <- setdiff(rate_diseases, count_diseases)
in_counts_only <- setdiff(count_diseases, rate_diseases)
in_both        <- intersect(rate_diseases, count_diseases)

cat("Diseases in rate data:", length(rate_diseases), "\n")
cat("Diseases in count data:", length(count_diseases), "\n")
cat("Matched:", length(in_both), "\n")

if (length(in_rates_only) > 0) {
  cat("\nIn RATES only (not in counts):\n")
  for (d in in_rates_only) cat("  -", d, "\n")
}
if (length(in_counts_only) > 0) {
  cat("\nIn COUNTS only (not in rates):\n")
  for (d in in_counts_only) cat("  -", d, "\n")
}

# =============================================================================
# 2. Year coverage by disease (rates data)
# =============================================================================
cat("\n=== 2. YEAR COVERAGE BY DISEASE (rates, AUS only) ===\n\n")

coverage <- rates %>%
  filter(state == "AUS", !is.na(rate_per_100k)) %>%
  group_by(disease) %>%
  summarise(
    min_year = min(year),
    max_year = max(year),
    n_years = n(),
    .groups = "drop"
  ) %>%
  arrange(min_year, disease)

# Diseases with data starting after 2000
late_start <- coverage %>% filter(min_year > 2000)
cat("Diseases with data starting AFTER 2000 (may not have full baseline):\n")
if (nrow(late_start) > 0) {
  print(late_start, n = 30)
} else {
  cat("  None\n")
}

# Diseases with data ending before 2024
early_end <- coverage %>% filter(max_year < 2024)
cat("\nDiseases with data ending BEFORE 2024:\n")
if (nrow(early_end) > 0) {
  print(early_end, n = 30)
} else {
  cat("  None\n")
}

# Diseases with insufficient baseline (fewer than 3 years in 2015-2019)
baseline_check <- rates %>%
  filter(state == "AUS", year >= 2015, year <= 2019, !is.na(rate_per_100k)) %>%
  group_by(disease) %>%
  summarise(baseline_years = n(), .groups = "drop") %>%
  filter(baseline_years < 5)

cat("\nDiseases with INCOMPLETE baseline (< 5 years in 2015-2019, AUS):\n")
if (nrow(baseline_check) > 0) {
  print(baseline_check, n = 30)
} else {
  cat("  All diseases have full 5-year baselines\n")
}

# =============================================================================
# 3. Missing data patterns (rates, by state and year)
# =============================================================================
cat("\n=== 3. MISSING DATA PATTERNS ===\n\n")

# Focus on analysis period (2015-2025)
missing_analysis <- rates %>%
  filter(year >= 2015, year <= 2025) %>%
  group_by(disease, state) %>%
  summarise(
    n_total = n(),
    n_missing = sum(is.na(rate_per_100k)),
    pct_missing = round(n_missing / n_total * 100, 1),
    .groups = "drop"
  ) %>%
  filter(n_missing > 0)

cat("Disease-state combos with NAs in 2015-2025:\n")
if (nrow(missing_analysis) > 0) {
  # Summarise by disease
  missing_summary <- missing_analysis %>%
    group_by(disease) %>%
    summarise(
      states_with_missing = paste(state, collapse = ", "),
      total_missing_cells = sum(n_missing),
      .groups = "drop"
    ) %>%
    arrange(desc(total_missing_cells))
  print(missing_summary, n = 40)
} else {
  cat("  No missing values in analysis period\n")
}

# =============================================================================
# 4. Zero and very low rates (potential reporting gaps vs genuine zeros)
# =============================================================================
cat("\n=== 4. ZERO / VERY LOW RATES IN ANALYSIS PERIOD ===\n\n")

zero_rates <- rates %>%
  filter(year >= 2015, year <= 2025, state == "AUS", rate_per_100k == 0) %>%
  select(disease, year, rate_per_100k) %>%
  arrange(disease, year)

cat("AUS national rate = 0 (2015-2025):\n")
if (nrow(zero_rates) > 0) {
  # Summarise
  zero_summary <- zero_rates %>%
    group_by(disease) %>%
    summarise(
      zero_years = paste(year, collapse = ", "),
      n_zero = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(n_zero))
  print(zero_summary, n = 30)
} else {
  cat("  None\n")
}

# =============================================================================
# 5. Outlier detection (year-on-year changes > 300%)
# =============================================================================
cat("\n=== 5. LARGE YEAR-ON-YEAR CHANGES (>300%, AUS, rate > 1) ===\n\n")

yoy <- rates %>%
  filter(state == "AUS", !is.na(rate_per_100k)) %>%
  arrange(disease, year) %>%
  group_by(disease) %>%
  mutate(
    prev_rate = lag(rate_per_100k),
    pct_change = round((rate_per_100k - prev_rate) / prev_rate * 100, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(pct_change), abs(pct_change) > 300,
         rate_per_100k > 1 | prev_rate > 1) %>%
  select(disease, year, prev_rate, rate_per_100k, pct_change) %>%
  arrange(desc(abs(pct_change)))

cat("Extreme year-on-year changes:\n")
if (nrow(yoy) > 0) {
  print(yoy, n = 30)
} else {
  cat("  None found\n")
}

# =============================================================================
# 6. Cross-check: rate data vs count data for 2019
# =============================================================================
cat("\n=== 6. CROSS-CHECK: Rate vs Count data (2019, AUS) ===\n\n")

rate_2019 <- rates %>%
  filter(year == 2019, state == "AUS", !is.na(rate_per_100k)) %>%
  select(disease, rate_per_100k)

count_2019 <- counts %>%
  filter(year == 2019, !is.na(count)) %>%
  select(disease, count)

# Use approximate ABS population for 2019 (~25.4M)
cross <- inner_join(rate_2019, count_2019, by = "disease") %>%
  mutate(
    implied_rate = round(count / 25400000 * 100000, 1),
    rate_diff = round(rate_per_100k - implied_rate, 1),
    pct_diff = round(rate_diff / implied_rate * 100, 1)
  ) %>%
  filter(abs(pct_diff) > 10) %>%
  arrange(desc(abs(pct_diff)))

cat("Diseases where rate and implied rate differ by >10% (2019):\n")
if (nrow(cross) > 0) {
  print(cross %>% select(disease, rate_per_100k, count, implied_rate, pct_diff), n = 20)
} else {
  cat("  All consistent within 10%\n")
}

# =============================================================================
# 7. Caveats affecting analysis period
# =============================================================================
cat("\n=== 7. DATA CAVEATS AFFECTING 2015-2025 ===\n\n")

# Parse year ranges from caveats
relevant_caveats <- caveats %>%
  mutate(
    # Extract all 4-digit years mentioned
    years_mentioned = str_extract_all(notes, "[0-9]{4}")
  ) %>%
  rowwise() %>%
  mutate(
    max_year_mentioned = if (length(years_mentioned) > 0) max(as.integer(years_mentioned)) else NA_integer_
  ) %>%
  ungroup() %>%
  # Keep caveats that mention years >= 2010 or are general ("nationally notifiable")
  filter(is.na(max_year_mentioned) | max_year_mentioned >= 2010 |
         str_detect(notes, "(?i)nationally notifiable"))

cat("Caveats potentially relevant to analysis period:\n")
if (nrow(relevant_caveats) > 0) {
  for (i in seq_len(nrow(relevant_caveats))) {
    cat(sprintf("  [%s | %s] %s\n",
                relevant_caveats$disease[i],
                relevant_caveats$jurisdiction[i],
                relevant_caveats$notes[i]))
  }
} else {
  cat("  None found\n")
}

# =============================================================================
# 8. Summary: diseases suitable for analysis
# =============================================================================
cat("\n=== 8. ANALYSIS SUITABILITY SUMMARY ===\n\n")

suitability <- rates %>%
  filter(state == "AUS") %>%
  group_by(disease) %>%
  summarise(
    min_year = min(year[!is.na(rate_per_100k)], na.rm = TRUE),
    max_year = max(year[!is.na(rate_per_100k)], na.rm = TRUE),
    baseline_years = sum(year >= 2015 & year <= 2019 & !is.na(rate_per_100k)),
    covid_years = sum(year >= 2020 & year <= 2021 & !is.na(rate_per_100k)),
    post_years = sum(year >= 2022 & year <= 2025 & !is.na(rate_per_100k)),
    rate_2019 = rate_per_100k[year == 2019],
    rate_2020 = rate_per_100k[year == 2020],
    .groups = "drop"
  ) %>%
  mutate(
    has_full_baseline = baseline_years == 5,
    has_covid_data = covid_years == 2,
    has_post_data = post_years >= 2,
    suitable = has_full_baseline & has_covid_data & has_post_data
  )

cat("Diseases with FULL baseline + COVID + post-COVID data:", sum(suitability$suitable), "\n")
cat("Diseases missing some periods:", sum(!suitability$suitable), "\n\n")

unsuitable <- suitability %>%
  filter(!suitable) %>%
  select(disease, min_year, baseline_years, covid_years, post_years)

if (nrow(unsuitable) > 0) {
  cat("Diseases with INCOMPLETE period coverage:\n")
  print(unsuitable, n = 30)
}

cat("\n--- Complete list of analysis-suitable diseases ---\n")
suitable_diseases <- suitability %>%
  filter(suitable) %>%
  mutate(
    pct_change_2020 = round((rate_2020 - rate_2019) / rate_2019 * 100, 1)
  ) %>%
  arrange(pct_change_2020) %>%
  select(disease, rate_2019, rate_2020, pct_change_2020)

print(suitable_diseases, n = 70)

cat("\n##############################################################\n")
cat("#               END OF QUALITY REPORT                        #\n")
cat("##############################################################\n")
