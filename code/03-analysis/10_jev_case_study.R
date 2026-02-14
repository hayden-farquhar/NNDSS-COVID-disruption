#!/usr/bin/env Rscript
# =============================================================================
# 10_jev_case_study.R
#
# Japanese Encephalitis Virus (JEV) case study:
#   - Document the 2022 emergence in mainland Australia
#   - Historical context (sporadic before 2022, first mainland outbreak)
#   - Geographic/state distribution
#   - Context: emergence during pandemic disruption period
#
# JEV was excluded from the 47 analysis-suitable diseases (insufficient
# baseline), but is epidemiologically significant as a novel pathogen
# emergence during the pandemic period.
#
# Inputs:
#   data/processed/nndss_annual_all_years.csv      — full historical data
#   data/processed/nndss_annual_analysis.csv       — 2015-2025 state-level
#
# Outputs:
#   data/processed/jev_case_study.csv              — JEV annual profile
#   output/tables/jev_summary.csv                  — manuscript insert
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

outdir <- "data/processed"
tabdir <- "output/tables"

cat("=== JEV Case Study ===\n\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
all_years <- read_csv(file.path(outdir, "nndss_annual_all_years.csv"),
                      show_col_types = FALSE)
analysis  <- read_csv(file.path(outdir, "nndss_annual_analysis.csv"),
                      show_col_types = FALSE)

# =============================================================================
# 1. Historical JEV profile (national)
# =============================================================================
cat("--- JEV historical profile ---\n")

jev_national <- all_years %>%
  filter(str_detect(disease, "(?i)japanese"),
         state == "AUS") %>%
  select(disease, year, rate_per_100k, count, period) %>%
  arrange(year)

cat("JEV notifications (Australia), all available years:\n")
jev_national %>%
  filter(count > 0 | year >= 2015) %>%
  mutate(rate_per_100k = round(rate_per_100k, 3)) %>%
  print(n = 40)

cat("\nTotal JEV cases:\n")
jev_national %>%
  summarise(
    total_all_years = sum(count, na.rm = TRUE),
    pre_2022        = sum(count[year < 2022], na.rm = TRUE),
    year_2022       = sum(count[year == 2022], na.rm = TRUE),
    post_2022       = sum(count[year > 2022], na.rm = TRUE),
    years_with_cases = sum(count > 0, na.rm = TRUE)
  ) %>%
  print()

# =============================================================================
# 2. State-level distribution (2015-2025)
# =============================================================================
cat("\n--- JEV by state (2015-2025) ---\n")

jev_state <- analysis %>%
  filter(str_detect(disease, "(?i)japanese"),
         state != "AUS",
         year >= 2015) %>%
  select(disease, year, state, rate_per_100k) %>%
  pivot_wider(names_from = state, values_from = rate_per_100k,
              values_fill = 0)

cat("JEV state rates (per 100,000):\n")
jev_state %>%
  mutate(across(where(is.numeric) & !c(year), ~ round(., 3))) %>%
  print(n = 15)

# Which states reported JEV in 2022?
jev_2022_states <- analysis %>%
  filter(str_detect(disease, "(?i)japanese"),
         year == 2022,
         state != "AUS",
         rate_per_100k > 0) %>%
  select(state, rate_per_100k) %>%
  arrange(desc(rate_per_100k))

cat("\nStates with JEV cases in 2022:\n")
jev_2022_states %>%
  mutate(rate_per_100k = round(rate_per_100k, 3)) %>%
  print()

# =============================================================================
# 3. Context: other emerging/re-emerging vector-borne diseases
# =============================================================================
cat("\n--- Context: other vector-borne diseases during pandemic ---\n")

# Compare JEV emergence to other arboviruses
arbovirus_comparison <- all_years %>%
  filter(state == "AUS",
         transmission_mode == "vector_borne" |
           str_detect(disease, "(?i)japanese"),
         year >= 2015, year <= 2025) %>%
  select(disease, year, rate_per_100k, count) %>%
  group_by(disease) %>%
  mutate(
    baseline_mean = mean(count[year >= 2015 & year <= 2019], na.rm = TRUE)
  ) %>%
  ungroup()

cat("Vector-borne disease comparison (2022 counts):\n")
arbovirus_comparison %>%
  filter(year == 2022) %>%
  select(disease, count, baseline_mean) %>%
  mutate(
    oe_ratio = round(count / pmax(1, baseline_mean), 3),
    baseline_mean = round(baseline_mean, 1)
  ) %>%
  arrange(desc(count)) %>%
  print(n = 15)

# =============================================================================
# 4. Narrative summary
# =============================================================================
cat("\n--- JEV Case Study Narrative ---\n\n")

total_2022 <- jev_national %>% filter(year == 2022) %>% pull(count)
pre_total  <- jev_national %>% filter(year < 2022) %>%
  summarise(sum(count, na.rm = TRUE)) %>% pull()
n_states_2022 <- nrow(jev_2022_states)

cat("Japanese Encephalitis Virus (JEV) Case Study\n\n")
cat("Prior to 2022, JEV in Australia was rare and confined to the Torres Strait\n")
cat("Islands and Far North Queensland. A total of", pre_total,
    "cases were notified\n")
cat("in the entire pre-2022 period.\n\n")
cat("In 2022, during the pandemic transition period, Australia experienced its\n")
cat("first mainland JEV outbreak with", total_2022, "notifications across",
    n_states_2022, "states/territories.\n")
cat("This represented an unprecedented geographic expansion of JEV into temperate\n")
cat("southeastern Australia, associated with flooding events and porcine hosts.\n\n")
cat("The emergence coincided with:\n")
cat("  - Easing of COVID-19 restrictions and border reopenings\n")
cat("  - La Nina-associated flooding in eastern Australia\n")
cat("  - Possible reduced vector surveillance during pandemic disruption\n\n")
cat("Post-2022, JEV notifications returned to near-zero (2023: 0 cases,\n")
cat("2024:", jev_national %>% filter(year == 2024) %>% pull(count),
    "case(s)),\n")
cat("suggesting the 2022 outbreak may have been a transient incursion.\n")
cat("However, the 2025 count of",
    jev_national %>% filter(year == 2025) %>% pull(count),
    "cases suggests ongoing low-level activity.\n")

# =============================================================================
# 5. Save outputs
# =============================================================================
write_csv(jev_national, file.path(outdir, "jev_case_study.csv"))
cat("\nSaved: jev_case_study.csv\n")

jev_summary <- tibble(
  metric = c("Total cases (all years)", "Cases pre-2022",
             "Cases in 2022", "States affected 2022",
             "Cases 2023", "Cases 2024", "Cases 2025",
             "Historical annual mean (2015-2019)",
             "2022 O/E ratio"),
  value = c(
    sum(jev_national$count, na.rm = TRUE),
    pre_total,
    total_2022,
    n_states_2022,
    jev_national %>% filter(year == 2023) %>% pull(count),
    jev_national %>% filter(year == 2024) %>% pull(count),
    jev_national %>% filter(year == 2025) %>% pull(count),
    round(mean(jev_national$count[jev_national$year >= 2015 &
                                    jev_national$year <= 2019],
               na.rm = TRUE), 1),
    round(total_2022 / pmax(1, mean(jev_national$count[
      jev_national$year >= 2015 & jev_national$year <= 2019],
      na.rm = TRUE)), 1)
  )
)

write_csv(jev_summary, file.path(tabdir, "jev_summary.csv"))
cat("Saved: jev_summary.csv\n")

cat("\n=== JEV case study complete ===\n")
