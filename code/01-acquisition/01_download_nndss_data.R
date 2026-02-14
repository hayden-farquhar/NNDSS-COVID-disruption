#!/usr/bin/env Rscript
# =============================================================================
# 01_download_nndss_data.R
#
# Download all available NNDSS data from Australian CDC and data.gov.au.
#
# Data sources:
#   1. Record-level datasets (4 diseases, XLSX, 2009-2024) from cdc.gov.au
#   2. Fortnightly surveillance reports (XLSX, 2024-2025) from cdc.gov.au
#      — contain all 70+ diseases by state with YTD totals
#   3. Contextual data: Oxford Stringency Index, ABS population
#
# For historical aggregate data (2010-2023, all diseases by state/year),
# this script also attempts to query the NINDSS Power BI dashboard API.
# If that fails, it prints instructions for manual export.
#
# Output:
#   data/raw/record_level/     — 4 disease XLSX files
#   data/raw/fortnightly/      — 24 fortnightly report XLSX files
#   data/raw/contextual/       — stringency index, population data
#   data/raw/powerbi/          — Power BI export (if API works)
#
# Usage:
#   Rscript code/01-acquisition/01_download_nndss_data.R
# =============================================================================

library(httr)
library(rvest)
library(tidyverse)
library(cli)
library(readxl)
library(jsonlite)

# --- Configuration -----------------------------------------------------------

CDC_BASE <- "https://www.cdc.gov.au"
DELAY_SEC <- 1
TIMEOUT_SEC <- 60

# Output directories
RECORD_DIR     <- "data/raw/record_level"
FORTNIGHTLY_DIR <- "data/raw/fortnightly"
CONTEXTUAL_DIR <- "data/raw/contextual"
POWERBI_DIR    <- "data/raw/powerbi"

for (d in c(RECORD_DIR, FORTNIGHTLY_DIR, CONTEXTUAL_DIR, POWERBI_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Helper: download a file with resume support -----------------------------
# Uses system curl with --http1.1 to avoid HTTP/2 stream errors on gov.au sites

download_file <- function(url, destfile, description = basename(destfile)) {
  if (file.exists(destfile) && file.size(destfile) > 0) {
    cli_alert_info("  Skipping {description} (already exists)")
    return(TRUE)
  }

  cli_alert("  Downloading: {description}")

  # Use system curl for gov.au sites (HTTP/2 issues with httr)
  use_system_curl <- str_detect(url, "gov\\.au")

  tryCatch({
    if (use_system_curl) {
      exit_code <- system2(
        "curl",
        args = c(
          "--http1.1", "-s", "-L",
          "-o", shQuote(destfile),
          "-w", "'%{http_code}'",
          "--max-time", as.character(TIMEOUT_SEC),
          shQuote(url)
        ),
        stdout = TRUE, stderr = TRUE
      )
      http_code <- as.integer(str_extract(paste(exit_code, collapse = ""), "\\d{3}"))
      ok <- file.exists(destfile) && file.size(destfile) > 100 &&
            (is.na(http_code) || http_code == 200)
    } else {
      resp <- GET(url, write_disk(destfile, overwrite = TRUE),
                  timeout(TIMEOUT_SEC), progress())
      ok <- status_code(resp) == 200 && file.size(destfile) > 0
    }

    if (ok) {
      size_kb <- round(file.size(destfile) / 1024, 1)
      cli_alert_success("  Saved {description} ({size_kb} KB)")
      return(TRUE)
    } else {
      cli_alert_danger("  Failed to download {description}")
      unlink(destfile)
      return(FALSE)
    }
  }, error = function(e) {
    cli_alert_danger("  Error: {e$message}")
    unlink(destfile)
    return(FALSE)
  })
}


# =============================================================================
# SECTION 1: Record-level datasets (4 diseases)
# =============================================================================

download_record_level <- function() {
  cli_h1("Section 1: Record-level datasets")

  datasets <- tribble(
    ~disease,        ~filename,                                            ~url_path,
    "Influenza",     "nndss-public-dataset-influenza-laboratory-confirmed.xlsx",
                     "/system/files/2025-10/nndss-public-dataset-influenza-laboratory-confirmed.xlsx",
    "Meningococcal", "nndss-public-dataset-meningococcal-disease-invasive.xlsx",
                     "/system/files/2025-10/nndss-public-dataset-meningococcal-disease-invasive.xlsx",
    "Pneumococcal",  "nndss-public-dataset-pneumococcal-disease-invasive.xlsx",
                     "/system/files/2025-10/nndss-public-dataset-pneumococcal-disease-invasive.xlsx",
    "Salmonella",    "nndss-public-dataset-salmonella.xlsx",
                     "/system/files/2025-10/nndss-public-dataset-salmonella.xlsx"
  )

  results <- map2_lgl(
    datasets$url_path, datasets$filename,
    ~ download_file(
      url = paste0(CDC_BASE, .x),
      destfile = file.path(RECORD_DIR, .y),
      description = .y
    )
  )

  cli_alert_info("Record-level: {sum(results)}/{length(results)} downloaded")
  return(results)
}


# =============================================================================
# SECTION 2: Fortnightly surveillance report XLSXs
# =============================================================================

download_fortnightly <- function() {
  cli_h1("Section 2: Fortnightly report XLSXs")

  # Known XLSX URLs (discovered from cdc.gov.au, Dec 2024 - Nov 2025)
  xlsx_urls <- c(
    "https://www.cdc.gov.au/system/files/2025-09/nndss_fortnightly_table_-_09_december_2024_to_05_january_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/nndss_fortnightly_table_-_06_january_to_19_january_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn03_-_10feb2025-web.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn04_-_24feb2025-web.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn05_-_11mar2025-web.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn06_-_24mar2025-web_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn07_-_07apr2025-web_0_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn08_-_22apr2025-web_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/fortnight_fn09_-_05_may_2025_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/fortnight_fn10_28th_april_to_11th_may_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/fortnight_fn11_12_may_2025_to_25_may_2025_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn12_-_16june2025_-_web_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/3._fortnight_fn13_-_30june2025_-_web_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/nndss_fortnightly_reports_-_23_june_2025_to_6_july_2025_0.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/fortnight_fn15_28_july_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/fortnight_fn16_-_11_august_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/fortnight_fn17_-_25_august_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-09/fortnight_fn18_-_08_sept_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-11/fortnight_fn19_-_22_sept2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-11/fortnight_fn20_-_7_oct_2025.xlsx",
    "https://www.cdc.gov.au/system/files/2025-11/fortnight_fn21_-_20_oct_2025.xlsx",
    "https://www.cdc.gov.au/sites/default/files/2025-12/fortnight_fn22_-_03_nov_2025.xlsx",
    "https://www.cdc.gov.au/sites/default/files/2025-12/nndss_fortnightly_table_-_27_october_2025_to_9_november_2025.xlsx",
    "https://www.cdc.gov.au/sites/default/files/2025-12/fortnight_fn24_-_10th_november_2025_to_23rd_november_2025.xlsx"
  )

  results <- map_lgl(xlsx_urls, function(url) {
    fname <- basename(url)
    ok <- download_file(url, file.path(FORTNIGHTLY_DIR, fname), fname)
    Sys.sleep(DELAY_SEC)
    ok
  })

  cli_alert_info("Fortnightly: {sum(results)}/{length(results)} downloaded")
  return(results)
}


# =============================================================================
# SECTION 3: Contextual data
# =============================================================================

download_contextual <- function() {
  cli_h1("Section 3: Contextual data")

  # --- Oxford COVID-19 Government Response Tracker (Stringency Index) ---
  cli_h2("Oxford Stringency Index")
  download_file(
    url = "https://raw.githubusercontent.com/OxCGRT/covid-policy-dataset/main/data/OxCGRT_compact_subnational_v1.csv",
    destfile = file.path(CONTEXTUAL_DIR, "oxcgrt_subnational.csv"),
    description = "Oxford Stringency Index (subnational)"
  )

  # Also get national-level as backup
  download_file(
    url = "https://raw.githubusercontent.com/OxCGRT/covid-policy-dataset/main/data/OxCGRT_compact_national_v1.csv",
    destfile = file.path(CONTEXTUAL_DIR, "oxcgrt_national.csv"),
    description = "Oxford Stringency Index (national)"
  )

  # --- ABS Population ---
  cli_h2("ABS Population estimates")
  # ABS Estimated Resident Population by state, quarterly
  # Table 04 from 3101.0 National, state and territory population
  download_file(
    url = "https://www.abs.gov.au/statistics/people/population/national-state-and-territory-population/latest-release/310104.xlsx",
    destfile = file.path(CONTEXTUAL_DIR, "abs_population_3101_table4.xlsx"),
    description = "ABS Population by state (3101.0 Table 4)"
  )

  # --- Google Mobility (archived) ---
  cli_h2("Google Mobility data")
  download_file(
    url = "https://www.gstatic.com/covid19/mobility/Region_Mobility_Report_CSVs.zip",
    destfile = file.path(CONTEXTUAL_DIR, "google_mobility.zip"),
    description = "Google Mobility Reports (ZIP)"
  )
}


# =============================================================================
# SECTION 4: Power BI dashboard data (historical aggregate)
# =============================================================================

attempt_powerbi_export <- function() {
  cli_h1("Section 4: Power BI dashboard (historical aggregate data)")

  pbi_url <- "https://nindss.health.gov.au/pbi-dashboard/"

  cli_alert("Fetching Power BI dashboard page for embed token...")

  tryCatch({
    resp <- GET(pbi_url, timeout(TIMEOUT_SEC))
    if (status_code(resp) != 200) {
      cli_alert_warning("Dashboard returned HTTP {status_code(resp)}")
      return(FALSE)
    }

    page_text <- content(resp, as = "text", encoding = "UTF-8")

    # Extract the embedConfig base64 blob
    config_match <- str_extract(
      page_text,
      'embedConfig="([A-Za-z0-9+/=]+)"'
    )

    if (is.na(config_match)) {
      cli_alert_warning("Could not find embedConfig on dashboard page")
      return(FALSE)
    }

    config_b64 <- str_extract(config_match, '(?<=embedConfig=")[^"]+')
    config_json <- rawToChar(base64enc::base64decode(config_b64))

    # Parse the JSON (may be truncated, try to extract key fields)
    # Extract report ID, group ID, and embed token
    report_id <- str_extract(config_json, '"Id":"([^"]+)"') %>%
      str_extract('[0-9a-f-]{36}')
    embed_url <- str_extract(config_json, '"EmbedUrl":"([^"]+)"') %>%
      str_remove('^"EmbedUrl":"') %>% str_remove('"$')

    # Extract the token
    token_match <- str_extract(config_json, '"token":"([^"]+)"')
    embed_token <- str_remove(token_match, '"token":"') %>% str_remove('"$')

    if (is.na(report_id) || is.na(embed_token)) {
      cli_alert_warning("Could not parse Power BI config")
      return(FALSE)
    }

    cli_alert_success("Report ID: {report_id}")
    cli_alert_info("Token length: {nchar(embed_token)}")

    # Save the config for reference
    writeLines(config_json, file.path(POWERBI_DIR, "embed_config.json"))

    # Attempt to query the Power BI API
    # The cluster URL for Australia Southeast
    cluster_url <- "https://wabi-australia-southeast-redirect.analysis.windows.net"

    # Try the public report querydata endpoint
    query_url <- paste0(
      cluster_url,
      "/public/reports/querydata?synchronous=true"
    )

    cli_alert("Attempting Power BI API query...")

    # DAX query to get all notification data
    query_body <- list(
      version = "1.0.0",
      queries = list(
        list(
          Query = list(
            Commands = list(
              list(
                SemanticQueryDataShapeCommand = list(
                  Query = list(
                    Version = 2L,
                    From = list(
                      list(Name = "n", Entity = "Notifications", Type = 0L)
                    ),
                    Select = list(
                      list(
                        Column = list(
                          Expression = list(SourceRef = list(Source = "n")),
                          Property = "Disease"
                        ),
                        Name = "Disease"
                      )
                    )
                  )
                )
              )
            )
          )
        )
      ),
      cancelQueries = list(),
      modelId = 0L
    )

    pbi_resp <- POST(
      query_url,
      add_headers(
        Authorization = paste("EmbedToken", embed_token),
        `Content-Type` = "application/json"
      ),
      body = toJSON(query_body, auto_unbox = TRUE),
      timeout(TIMEOUT_SEC)
    )

    if (status_code(pbi_resp) == 200) {
      result_text <- content(pbi_resp, as = "text", encoding = "UTF-8")
      writeLines(result_text, file.path(POWERBI_DIR, "query_result.json"))
      cli_alert_success("Power BI query succeeded! Saved to powerbi/query_result.json")
      cli_alert_info("Run code/01-acquisition/02_parse_powerbi.R to extract the data")
      return(TRUE)
    } else {
      cli_alert_warning(
        "Power BI API returned HTTP {status_code(pbi_resp)}"
      )
      result_text <- content(pbi_resp, as = "text", encoding = "UTF-8")
      cli_alert_info("Response: {str_trunc(result_text, 200)}")
      return(FALSE)
    }

  }, error = function(e) {
    cli_alert_warning("Power BI attempt failed: {e$message}")
    return(FALSE)
  })
}


print_manual_export_instructions <- function() {
  cli_h2("Manual export instructions for historical aggregate data")
  cat("
  The Power BI API query did not succeed. To get historical aggregate data
  for all 70+ diseases (2010-2025), please export manually from the NINDSS
  Power BI dashboard:

  1. Open: https://nindss.health.gov.au/pbi-dashboard/
  2. In the dashboard, look for a table or chart showing notification counts
  3. Click the '...' menu (three dots) on the visual

  4. Select 'Export data' > 'Underlying data' or 'Summarized data'
  5. Choose CSV format
  6. Save as: data/raw/powerbi/nndss_powerbi_export.csv

  Alternative: Use the 'Analyze in Excel' option if available.

  The exported file should contain columns like:
  Disease, State, Year, Month, Notification Count

  Once saved, the cleaning script will automatically detect and process it.
  ")
}


# =============================================================================
# SECTION 5: Quick-inspect downloaded fortnightly data
# =============================================================================

inspect_fortnightly <- function() {
  cli_h2("Inspecting fortnightly data")

  xlsx_files <- list.files(FORTNIGHTLY_DIR, pattern = "\\.xlsx$",
                           full.names = TRUE)

  if (length(xlsx_files) == 0) {
    cli_alert_warning("No fortnightly XLSX files found")
    return(NULL)
  }

  # Read the latest one to show structure
  latest <- sort(xlsx_files, decreasing = TRUE)[1]
  cli_alert_info("Inspecting: {basename(latest)}")

  tryCatch({
    df <- read_excel(latest, sheet = 1, col_names = FALSE)
    cli_alert_info("  Dimensions: {nrow(df)} rows x {ncol(df)} cols")

    # Find disease name column (usually column 2)
    diseases <- df %>%
      pull(2) %>%
      na.omit() %>%
      unique() %>%
      discard(~ str_detect(.x, "(?i)disease|name|total|footnote|data"))

    cli_alert_success("  Found {length(diseases)} diseases")
    cat("  First 10: ", paste(head(diseases, 10), collapse = ", "), "\n")

  }, error = function(e) {
    cli_alert_warning("  Could not inspect: {e$message}")
  })
}


# =============================================================================
# Main execution
# =============================================================================

main <- function() {
  cli_h1("NNDSS Data Acquisition")
  cli_alert_info("Timestamp: {Sys.time()}")
  cli_rule()

  # Section 1: Record-level datasets
  download_record_level()

  # Section 2: Fortnightly reports
  download_fortnightly()

  # Section 3: Contextual data
  download_contextual()

  # Section 4: Power BI (historical aggregate)
  pbi_ok <- attempt_powerbi_export()
  if (!pbi_ok) {
    print_manual_export_instructions()
  }

  # Section 5: Quick inspection
  inspect_fortnightly()

  # Summary
  cli_h1("Download Summary")

  record_files <- list.files(RECORD_DIR, pattern = "\\.xlsx$")
  fort_files   <- list.files(FORTNIGHTLY_DIR, pattern = "\\.xlsx$")
  ctx_files    <- list.files(CONTEXTUAL_DIR)
  pbi_files    <- list.files(POWERBI_DIR)

  cli_alert_info("Record-level datasets: {length(record_files)} files")
  cli_alert_info("Fortnightly reports:   {length(fort_files)} files")
  cli_alert_info("Contextual data:       {length(ctx_files)} files")
  cli_alert_info("Power BI exports:      {length(pbi_files)} files")

  # Save metadata
  metadata <- list(
    download_date = as.character(Sys.time()),
    record_level = list(
      source = "https://www.cdc.gov.au/resources/collections/nndss-public-datasets",
      diseases = c("influenza", "meningococcal", "pneumococcal", "salmonella"),
      coverage = "2009-2024",
      files = record_files
    ),
    fortnightly = list(
      source = "https://www.cdc.gov.au/resources/collections/nndss-fortnightly-reports",
      coverage = "Dec 2024 - Nov 2025",
      files = fort_files
    ),
    contextual = list(
      oxcgrt = "https://github.com/OxCGRT/covid-policy-dataset",
      abs_pop = "https://www.abs.gov.au/statistics/people/population",
      files = ctx_files
    ),
    powerbi = list(
      source = "https://nindss.health.gov.au/pbi-dashboard/",
      note = "Contains all 71 diseases, all states, historical data updated daily",
      api_success = pbi_ok,
      files = pbi_files
    )
  )

  write_json(metadata, "data/metadata.json", pretty = TRUE, auto_unbox = TRUE)
  cli_alert_success("Metadata saved to data/metadata.json")

  cli_rule()
  cli_alert_info("Next step: Run code/02-cleaning/01_clean_fortnightly.R")
}

# Run
main()
