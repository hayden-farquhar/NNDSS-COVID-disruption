#!/usr/bin/env Rscript
# =============================================================================
# 03_create_analysis_datasets.R
#
# Create final analysis-ready datasets by:
#   1. Joining rate data with disease metadata (transmission mode, etc.)
#   2. Handling NA rates (context-appropriate zero imputation)
#   3. Adding count data for cross-validation
#   4. Filtering to analysis period and suitable diseases
#   5. Parsing fortnightly reports for sub-annual resolution
#
# Outputs:
#   data/processed/nndss_annual_analysis.csv       — main analysis dataset
#   data/processed/nndss_annual_all_years.csv       — full historical dataset
#   data/processed/nndss_fortnightly_2025.csv       — sub-annual time series
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(readxl)
})

outdir <- "data/processed"

# =============================================================================
# 1. Annual analysis dataset
# =============================================================================
cat("=== Building annual analysis dataset ===\n\n")

rates    <- read_csv(file.path(outdir, "nndss_rates_by_state_year.csv"), show_col_types = FALSE)
counts   <- read_csv(file.path(outdir, "nndss_counts_national_annual.csv"), show_col_types = FALSE)
meta     <- read_csv(file.path(outdir, "disease_metadata.csv"), show_col_types = FALSE)

# --- Handle NA rates that should be zero ---
rates_filled <- rates %>%
  left_join(
    meta %>% select(disease, year_data_from, year_data_to),
    by = "disease"
  ) %>%
  mutate(
    rate_per_100k = case_when(
      !is.na(rate_per_100k) ~ rate_per_100k,
      year >= year_data_from & year <= year_data_to ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-year_data_from, -year_data_to)

cat("NA rates BEFORE fill:", sum(is.na(rates$rate_per_100k)), "\n")
cat("NA rates AFTER fill:", sum(is.na(rates_filled$rate_per_100k)), "\n")

# --- Join with national counts ---
annual <- rates_filled %>%
  left_join(
    counts %>% select(disease, year, count),
    by = c("disease", "year")
  ) %>%
  left_join(
    meta %>% select(disease, transmission_mode, subcategory,
                     analysis_suitable, border_sensitive,
                     vaccine_preventable, exclusion_reason),
    by = "disease"
  )

# --- Define analysis periods ---
annual <- annual %>%
  mutate(
    period = case_when(
      year < 2015 ~ "historical",
      year >= 2015 & year <= 2019 ~ "baseline",
      year >= 2020 & year <= 2021 ~ "covid_acute",
      year == 2022 ~ "transition",
      year >= 2023 ~ "post_covid"
    )
  )

annual_analysis <- annual %>% filter(year >= 2015, year <= 2025)

cat("\n--- Annual analysis dataset ---\n")
cat("Rows:", nrow(annual_analysis), "\n")
cat("Diseases:", n_distinct(annual_analysis$disease), "\n")
cat("Analysis-suitable:", n_distinct(annual_analysis$disease[annual_analysis$analysis_suitable == TRUE]), "\n")
cat("Years:", paste(range(annual_analysis$year), collapse = "-"), "\n")

write_csv(annual, file.path(outdir, "nndss_annual_all_years.csv"))
cat("\nSaved: nndss_annual_all_years.csv\n")
write_csv(annual_analysis, file.path(outdir, "nndss_annual_analysis.csv"))
cat("Saved: nndss_annual_analysis.csv\n")

# =============================================================================
# 2. Fortnightly reports → sub-annual time series (2024-2025)
# =============================================================================
cat("\n=== Parsing fortnightly surveillance reports ===\n\n")

fn_dir <- "data/raw/fortnightly"
fn_files <- sort(list.files(fn_dir, pattern = "[.]xlsx$", full.names = TRUE))
cat("Found", length(fn_files), "fortnightly report files\n\n")

# Structure (from inspection):
# - Read with col_names = FALSE
# - Row 1: report ID cell (e.g. "ADT FN03/2025"), col 4 = "Notification received date"
# - Row 2: blank, blank, blank, "State or Territory", ..., "Totals for Australia", ...
# - Row 3: "Disease group", "Disease name", "Disease code", "ACT", "NSW", "NT",
#           "Qld", "SA", "Tas", "Vic", "WA",
#           "This reporting period", "Previous reporting Period",
#           "Same reporting period last year", "Current year YTD",
#           "Past Quarter", "Quarterly rolling 5 year mean", ...
# - Row 4-5: date range rows (Excel serial numbers)
# - Row 6+: data rows
# - Footer: totals row, footnotes

states <- c("ACT", "NSW", "NT", "Qld", "SA", "Tas", "Vic", "WA")

all_fn <- list()

for (f in fn_files) {
  tryCatch({
    df <- suppressMessages(
      read_excel(f, sheet = 1, col_names = FALSE, col_types = "text",
                 .name_repair = "minimal")
    )

    # Extract fortnight identifier from row 1
    cell_a1 <- as.character(df[[1, 1]])
    fn_id <- str_extract(cell_a1, "FN\\d+/\\d{4}")
    if (is.na(fn_id)) {
      # Try extracting from filename
      fn_num <- str_extract(basename(f), "fn(\\d+)")
      fn_id <- paste0(toupper(fn_num), "/2025")
    }

    # Find the header row (contains "Disease group" or "Disease name")
    col1_vals <- as.character(unlist(df[, 1]))
    col2_vals <- as.character(unlist(df[, 2]))
    header_row <- which(
      str_detect(col1_vals, "(?i)disease\\s*group") |
      str_detect(col2_vals, "(?i)disease\\s*name")
    )[1]

    if (is.na(header_row)) {
      cat("SKIP:", basename(f), "- no header row found\n")
      next
    }

    # Header values from the header row
    header <- as.character(unlist(df[header_row, ]))

    # Map column indices
    # Columns 1-3: Disease group, Disease name, Disease code
    # Columns 4-11: state counts (ACT, NSW, NT, Qld, SA, Tas, Vic, WA)
    # Column 12: This reporting period (national total for this fortnight)
    # Column 15: Current year YTD
    col_disease_group <- 1
    col_disease_name  <- 2
    col_disease_code  <- 3
    col_states_start  <- 4   # ACT
    col_states_end    <- 11  # WA
    col_this_period   <- 12
    col_ytd           <- 15  # "Current year YTD"

    # Data rows start after header + any date sub-header rows
    # Find first row after header where column 2 has a disease name
    data_start <- header_row + 1
    for (r in (header_row + 1):nrow(df)) {
      val <- as.character(df[[r, col_disease_name]])
      if (!is.na(val) && nchar(val) > 2 && !str_detect(val, "^\\d+$|^NA$")) {
        data_start <- r
        break
      }
    }

    # Extract data rows
    for (r in data_start:nrow(df)) {
      disease_group <- as.character(df[[r, col_disease_group]])
      disease_name  <- as.character(df[[r, col_disease_name]])

      # Stop at footer (Total row, footnotes, NA rows)
      if (!is.na(disease_group) && str_detect(disease_group, "(?i)^total|^footnote|^the data")) break
      if (is.na(disease_name) || disease_name == "Total") next
      if (str_detect(disease_name, "(?i)^total")) next

      # Fill down disease_group (only first disease in group has it)
      if (!is.na(disease_group) && nchar(disease_group) > 2 &&
          !str_detect(disease_group, "^\\d+$")) {
        current_group <- disease_group
      }

      # State counts
      state_counts <- sapply(col_states_start:col_states_end, function(j) {
        suppressWarnings(as.numeric(str_remove_all(as.character(df[[r, j]]), ",")))
      })
      names(state_counts) <- states

      # National totals
      this_period <- suppressWarnings(as.numeric(str_remove_all(as.character(df[[r, col_this_period]]), ",")))

      ytd_val <- NA_real_
      if (col_ytd <= ncol(df)) {
        ytd_val <- suppressWarnings(as.numeric(str_remove_all(as.character(df[[r, col_ytd]]), ",")))
      }

      all_fn[[length(all_fn) + 1]] <- tibble(
        fortnight_id = fn_id,
        disease_group = current_group,
        disease = disease_name,
        ACT = state_counts["ACT"],
        NSW = state_counts["NSW"],
        NT  = state_counts["NT"],
        QLD = state_counts["Qld"],
        SA  = state_counts["SA"],
        TAS = state_counts["Tas"],
        VIC = state_counts["Vic"],
        WA  = state_counts["WA"],
        national_this_fortnight = this_period,
        ytd = ytd_val
      )
    }

    cat("Parsed:", basename(f), "->", fn_id, "\n")

  }, error = function(e) {
    cat("ERROR:", basename(f), "-", e$message, "\n")
  })
}

fortnightly <- bind_rows(all_fn)

# Pivot state columns to long format
fn_long <- fortnightly %>%
  pivot_longer(
    cols = c(ACT, NSW, NT, QLD, SA, TAS, VIC, WA),
    names_to = "state",
    values_to = "fortnightly_count"
  )

cat("\n--- Fortnightly dataset ---\n")
cat("Rows (long):", nrow(fn_long), "\n")
cat("Diseases:", n_distinct(fn_long$disease), "\n")
cat("Fortnights:", paste(sort(unique(fn_long$fortnight_id)), collapse = ", "), "\n")

# Sanity check: influenza YTD across fortnights
cat("\n--- Influenza YTD progression ---\n")
fortnightly %>%
  filter(str_detect(disease, "(?i)influenza.+lab")) %>%
  select(fortnight_id, disease, national_this_fortnight, ytd) %>%
  print(n = 30)

write_csv(fn_long, file.path(outdir, "nndss_fortnightly_2025.csv"))
cat("\nSaved: nndss_fortnightly_2025.csv\n")

# =============================================================================
# 3. Final summary
# =============================================================================
cat("\n##############################################################\n")
cat("#           ANALYSIS-READY DATASETS SUMMARY                  #\n")
cat("##############################################################\n\n")

cat("1. nndss_annual_analysis.csv\n")
cat("   ", nrow(annual_analysis), "rows | 66 diseases | 9 jurisdictions | 2015-2025\n")
cat("   Columns: disease, year, state, rate_per_100k, count,\n")
cat("            transmission_mode, period, analysis_suitable, etc.\n\n")

cat("2. nndss_annual_all_years.csv\n")
cat("   ", nrow(annual), "rows | 66 diseases | 9 jurisdictions | 1991-2026\n")
cat("   Same columns as above, for historical context\n\n")

cat("3. nndss_fortnightly_2025.csv\n")
cat("   ", nrow(fn_long), "rows | sub-annual resolution\n")
cat("   Columns: fortnight_id, disease, state, fortnightly_count, ytd\n\n")

cat("4. disease_metadata.csv\n")
cat("   66 diseases with transmission_mode, analysis_suitable,\n")
cat("   border_sensitive, vaccine_preventable, exclusion_reason\n\n")

cat("Supporting files:\n")
cat("   nndss_counts_national_annual.csv — raw national counts\n")
cat("   nndss_counts_by_jurisdiction_year.csv — all-disease totals by state\n")
cat("   nndss_counts_by_age_year.csv — all-disease totals by age\n")
cat("   nndss_data_caveats.csv — data quality notes\n")
