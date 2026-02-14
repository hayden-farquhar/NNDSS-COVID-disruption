#!/usr/bin/env Rscript
# =============================================================================
# 02_parse_rate_exports.R
#
# Parse the 66 per-disease Power BI rate exports into a single tidy CSV.
# Input:  data/raw/powerbi/rates/*.xlsx  (one file per disease)
# Output: data/processed/nndss_rates_by_state_year.csv
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

rates_dir <- "data/raw/powerbi/rates"
files <- list.files(rates_dir, pattern = "[.]xlsx$", full.names = TRUE)
cat("Processing", length(files), "files...\n\n")

all_data <- list()

for (f in files) {
  tryCatch({
    df <- suppressMessages(
      read_excel(f, sheet = 1, col_names = FALSE, col_types = "text",
                 .name_repair = "minimal")
    )

    # Disease name from "Applied filters" row at the bottom
    all_col1 <- as.character(unlist(df[, 1]))
    filter_idx <- which(str_detect(all_col1, "Applied filters"))
    if (length(filter_idx) == 0) next

    disease_name <- str_extract(all_col1[filter_idx[1]],
                                "DISEASE NAME is ([^\r\n]+)")
    disease_name <- str_remove(disease_name, "DISEASE NAME is ")
    disease_name <- str_trim(disease_name)
    if (is.na(disease_name)) next

    # Header from row 1: "Diagnosis Year", "ACT", "NSW", "NT", ...
    header_vals <- as.character(unlist(df[1, ]))

    # Build data frame with proper column names
    data_df <- as.data.frame(df[2:nrow(df), ], stringsAsFactors = FALSE)
    names(data_df) <- header_vals

    # Keep only rows where first column is a 4-digit year
    year_rows <- grepl("^[0-9]{4}$", data_df[["Diagnosis Year"]])
    data_df <- data_df[year_rows, ]

    if (nrow(data_df) == 0) next

    # Pivot to long format
    state_cols <- setdiff(names(data_df), "Diagnosis Year")

    long <- data_df %>%
      pivot_longer(
        cols = all_of(state_cols),
        names_to = "state",
        values_to = "rate_raw"
      ) %>%
      transmute(
        disease = disease_name,
        year = as.integer(.data[["Diagnosis Year"]]),
        state = state,
        rate_per_100k = suppressWarnings(as.numeric(str_remove_all(rate_raw, ",")))
      )

    all_data[[disease_name]] <- long

  }, error = function(e) {
    cat("ERROR:", basename(f), "-", e$message, "\n")
  })
}

combined <- bind_rows(all_data)

# --- Summary -----------------------------------------------------------------

cat("=== SUMMARY ===\n")
cat("Diseases:", length(unique(combined$disease)), "\n")
cat("Years:", paste(range(combined$year, na.rm = TRUE), collapse = "-"), "\n")
cat("States:", paste(sort(unique(combined$state)), collapse = ", "), "\n")
cat("Total rows:", nrow(combined), "\n")
cat("Non-NA rates:", sum(!is.na(combined$rate_per_100k)), "\n\n")

# Influenza sanity check
cat("--- Influenza (lab confirmed), AUS, 2015-2025 ---\n")
combined %>%
  filter(disease == "Influenza (laboratory confirmed)", state == "AUS") %>%
  filter(year >= 2015) %>%
  arrange(year) %>%
  print(n = 15)

# Top diseases by 2019 rate
cat("\n--- Top 10 diseases by 2019 national rate ---\n")
combined %>%
  filter(year == 2019, state == "AUS", !is.na(rate_per_100k)) %>%
  arrange(desc(rate_per_100k)) %>%
  select(disease, rate_per_100k) %>%
  print(n = 10)

# COVID disruption snapshot
cat("\n--- COVID disruption: 2019 vs 2020 (AUS, >20% change, rate > 0.5) ---\n")
combined %>%
  filter(year %in% c(2019, 2020), state == "AUS", !is.na(rate_per_100k)) %>%
  pivot_wider(names_from = year, values_from = rate_per_100k, names_prefix = "y") %>%
  mutate(pct_change = round((y2020 - y2019) / y2019 * 100, 1)) %>%
  filter(abs(pct_change) > 20, y2019 > 0.5) %>%
  arrange(pct_change) %>%
  print(n = 25)

# --- Save --------------------------------------------------------------------

outfile <- "data/processed/nndss_rates_by_state_year.csv"
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
write_csv(combined, outfile)
cat("\nSaved to:", outfile, "\n")
