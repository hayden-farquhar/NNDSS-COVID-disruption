#!/usr/bin/env Rscript
# =============================================================================
# 02_tier1_parse_record_level.R
#
# Parse 4 record-level XLSX files (influenza, meningococcal, pneumococcal,
# salmonella) into monthly national/state counts for ARIMA modelling.
#
# Inputs:
#   data/raw/record_level/*.xlsx (4 files)
#
# Outputs:
#   data/processed/tier1_monthly_national.csv  — disease × year × month
#   data/processed/tier1_monthly_by_state.csv  — disease × year × month × state
#   data/processed/tier1_annual_by_age.csv     — disease × year × age_group
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(lubridate)
})

rawdir <- "data/raw/record_level"
outdir <- "data/processed"

# Standardised state names for output
state_map <- c(
  "ACT" = "ACT", "NSW" = "NSW", "NT" = "NT",
  "Qld" = "QLD", "QLD" = "QLD",
  "SA" = "SA", "Tas" = "TAS", "TAS" = "TAS",
  "Vic" = "VIC", "VIC" = "VIC", "WA" = "WA",
  "Vic/Tas" = "VIC_TAS"
)

# =============================================================================
# 1. Influenza (4 sheets, weekly, Excel serial dates)
# =============================================================================
cat("=== Parsing Influenza ===\n")

flu_file <- file.path(rawdir, "nndss-public-dataset-influenza-laboratory-confirmed.xlsx")
flu_sheets <- excel_sheets(flu_file)
cat("Sheets:", paste(flu_sheets, collapse = ", "), "\n")

flu_all <- list()
for (sheet in flu_sheets) {
  cat("  Reading sheet:", sheet, "...")
  df <- read_excel(flu_file, sheet = sheet, skip = 1,
                   col_names = FALSE, col_types = "text",
                   .name_repair = "minimal")

  # Assign column names from what we know about the structure
  names(df) <- c("week_ending", "state", "age_group", "sex",
                 "indigenous_status", "type_subtype")[1:ncol(df)]

  # Remove embedded header rows (where state == "State")
  df <- df %>% filter(state != "State", !is.na(state))

  # Convert Excel serial dates to proper dates
  df <- df %>%
    mutate(
      week_serial = suppressWarnings(as.numeric(week_ending)),
      date = as.Date(week_serial, origin = "1899-12-30"),
      year = year(date),
      month = month(date)
    ) %>%
    filter(!is.na(date))

  cat(" rows:", nrow(df), " date range:", as.character(min(df$date)),
      "to", as.character(max(df$date)), "\n")
  flu_all[[sheet]] <- df
}

flu <- bind_rows(flu_all) %>%
  mutate(
    disease = "Influenza (laboratory confirmed)",
    state_std = recode(state, !!!state_map)
  )
cat("Total influenza records:", nrow(flu), "\n\n")

# =============================================================================
# 2. Meningococcal (1 sheet, monthly, Year + Month text)
# =============================================================================
cat("=== Parsing Meningococcal ===\n")

men_file <- file.path(rawdir, "nndss-public-dataset-meningococcal-disease-invasive.xlsx")
men_sheets <- excel_sheets(men_file)
cat("Sheet:", men_sheets[1], "\n")

men <- read_excel(men_file, sheet = 1, skip = 1,
                  col_names = FALSE, col_types = "text",
                  .name_repair = "minimal")

# Assign known column names
men_cols <- c("year", "month_name", "state", "serotype", "age_group", "sex",
              "indigenous_status", "vaccination_status", "total_vaccinations",
              paste0("vax_", 1:15))
names(men) <- men_cols[1:ncol(men)]

# Remove embedded header row
men <- men %>% filter(year != "Year", !is.na(year))

# Parse month name to number
month_lookup <- setNames(1:12, month.name)

men <- men %>%
  mutate(
    disease = "Meningococcal disease (invasive)",
    year = as.integer(year),
    month = month_lookup[month_name],
    state_std = recode(state, !!!state_map)
  ) %>%
  filter(!is.na(month), !is.na(year))

cat("Meningococcal records:", nrow(men),
    " years:", min(men$year), "-", max(men$year), "\n\n")

# =============================================================================
# 3. Pneumococcal (1 sheet, annual only — no month column)
# =============================================================================
cat("=== Parsing Pneumococcal ===\n")

pneu_file <- file.path(rawdir, "nndss-public-dataset-pneumococcal-disease-invasive.xlsx")
pneu_sheets <- excel_sheets(pneu_file)
cat("Sheet:", pneu_sheets[1], "\n")

pneu <- read_excel(pneu_file, sheet = 1, skip = 1,
                   col_names = FALSE, col_types = "text",
                   .name_repair = "minimal")

# Column structure: Year, Serotype, State, Age group, Sex, Indigenous status,
# Clinical category, Other clinical category, Vaccination history,
# Total vaccinations, + vaccine brand columns
pneu_cols <- c("year", "serotype", "state", "age_group", "sex",
               "indigenous_status", "clinical_category",
               "other_clinical_category", "vaccination_history",
               "total_vaccinations", paste0("vax_", 1:5))
names(pneu) <- pneu_cols[1:ncol(pneu)]

pneu <- pneu %>%
  filter(year != "Year", !is.na(year)) %>%
  mutate(
    disease = "Pneumococcal disease (invasive)",
    year = as.integer(year),
    state_std = recode(state, !!!state_map)
  )

cat("Pneumococcal records:", nrow(pneu),
    " years:", min(pneu$year), "-", max(pneu$year), "\n")
cat("  NOTE: No month column — annual only. Cannot build monthly time series.\n")
cat("  Will aggregate to annual counts for age-stratified analysis.\n\n")

# =============================================================================
# 4. Salmonella (1 sheet, weekly, DD/MM/YYYY string)
# =============================================================================
cat("=== Parsing Salmonella ===\n")

sal_file <- file.path(rawdir, "nndss-public-dataset-salmonella.xlsx")
sal_sheets <- excel_sheets(sal_file)
cat("Sheet:", sal_sheets[1], "\n")

# Salmonella has proper headers (no title row to skip)
sal <- read_excel(sal_file, sheet = 1, col_types = "text")

names(sal) <- c("week_ending", "state", "age_group", "sex", "serovar")[1:ncol(sal)]

sal <- sal %>%
  mutate(
    disease = "Salmonellosis",
    # Date is DD/MM/YYYY character string
    date = as.Date(week_ending, format = "%d/%m/%Y"),
    year = year(date),
    month = month(date),
    state_std = recode(state, !!!state_map)
  ) %>%
  filter(!is.na(date))

cat("Salmonella records:", nrow(sal),
    " years:", min(sal$year), "-", max(sal$year), "\n\n")

# =============================================================================
# 5. Aggregate to monthly national counts (3 diseases with temporal data)
# =============================================================================
cat("=== Aggregating to monthly counts ===\n")

# Influenza: weekly → monthly
flu_monthly <- flu %>%
  group_by(disease, year, month) %>%
  summarise(count = n(), .groups = "drop")

# Meningococcal: already has year + month
men_monthly <- men %>%
  group_by(disease, year, month) %>%
  summarise(count = n(), .groups = "drop")

# Salmonella: weekly → monthly
sal_monthly <- sal %>%
  group_by(disease, year, month) %>%
  summarise(count = n(), .groups = "drop")

# Combine (3 diseases — pneumococcal excluded, no monthly data)
monthly_national <- bind_rows(flu_monthly, men_monthly, sal_monthly) %>%
  arrange(disease, year, month)

# Fill in missing months with zero
complete_grid <- monthly_national %>%
  group_by(disease) %>%
  reframe(
    year  = rep(min(year):max(year), each = 12),
    month = rep(1:12, times = max(year) - min(year) + 1)
  )

monthly_national <- complete_grid %>%
  left_join(monthly_national, by = c("disease", "year", "month")) %>%
  mutate(count = replace_na(count, 0))

cat("Monthly national rows:", nrow(monthly_national), "\n")
cat("Diseases:", paste(unique(monthly_national$disease), collapse = "; "), "\n")

# Quick sanity: influenza 2019 annual total should be ~313k
flu_2019 <- monthly_national %>%
  filter(disease == "Influenza (laboratory confirmed)", year == 2019) %>%
  summarise(total = sum(count))
cat("Influenza 2019 total from record-level:", flu_2019$total,
    "(Power BI reported: 313463)\n")

# =============================================================================
# 6. Aggregate to monthly by state
# =============================================================================
cat("\n--- Monthly by state ---\n")

flu_monthly_state <- flu %>%
  group_by(disease, year, month, state = state_std) %>%
  summarise(count = n(), .groups = "drop")

men_monthly_state <- men %>%
  group_by(disease, year, month, state = state_std) %>%
  summarise(count = n(), .groups = "drop")

sal_monthly_state <- sal %>%
  group_by(disease, year, month, state = state_std) %>%
  summarise(count = n(), .groups = "drop")

monthly_by_state <- bind_rows(flu_monthly_state, men_monthly_state,
                              sal_monthly_state) %>%
  arrange(disease, state, year, month)

cat("Monthly by state rows:", nrow(monthly_by_state), "\n")
cat("States:", paste(sort(unique(monthly_by_state$state)), collapse = ", "), "\n")

# =============================================================================
# 7. Annual counts by age group (all 4 diseases, for Phase 4)
# =============================================================================
cat("\n--- Annual by age group ---\n")

flu_age <- flu %>%
  group_by(disease, year, age_group) %>%
  summarise(count = n(), .groups = "drop")

men_age <- men %>%
  group_by(disease, year, age_group) %>%
  summarise(count = n(), .groups = "drop")

pneu_age <- pneu %>%
  group_by(disease, year, age_group) %>%
  summarise(count = n(), .groups = "drop")

sal_age <- sal %>%
  group_by(disease, year, age_group) %>%
  summarise(count = n(), .groups = "drop")

annual_by_age <- bind_rows(flu_age, men_age, pneu_age, sal_age) %>%
  arrange(disease, year, age_group)

cat("Annual by age rows:", nrow(annual_by_age), "\n")
cat("Diseases:", paste(unique(annual_by_age$disease), collapse = "; "), "\n")

# =============================================================================
# 8. Save outputs
# =============================================================================
write_csv(monthly_national, file.path(outdir, "tier1_monthly_national.csv"))
cat("\nSaved: tier1_monthly_national.csv (", nrow(monthly_national), "rows )\n")

write_csv(monthly_by_state, file.path(outdir, "tier1_monthly_by_state.csv"))
cat("Saved: tier1_monthly_by_state.csv (", nrow(monthly_by_state), "rows )\n")

write_csv(annual_by_age, file.path(outdir, "tier1_annual_by_age.csv"))
cat("Saved: tier1_annual_by_age.csv (", nrow(annual_by_age), "rows )\n")

# Summary table
cat("\n=== Record-level data summary ===\n")
tibble(
  disease = c("Influenza", "Meningococcal", "Pneumococcal", "Salmonella"),
  total_records = c(nrow(flu), nrow(men), nrow(pneu), nrow(sal)),
  temporal_resolution = c("weekly", "monthly", "annual only", "weekly"),
  monthly_ts_available = c(TRUE, TRUE, FALSE, TRUE),
  year_range = c(
    paste(range(flu$year), collapse = "-"),
    paste(range(men$year), collapse = "-"),
    paste(range(pneu$year), collapse = "-"),
    paste(range(sal$year), collapse = "-")
  )
) %>% print()

cat("\n=== Tier 1 parsing complete ===\n")
