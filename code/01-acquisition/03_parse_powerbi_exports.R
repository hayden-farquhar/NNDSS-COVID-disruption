#!/usr/bin/env Rscript
# =============================================================================
# 03_parse_powerbi_exports.R
#
# Parse the remaining Power BI exports into tidy CSVs:
#   1. National annual counts by disease  → nndss_counts_national_annual.csv
#   2. Notifications by jurisdiction      → nndss_counts_by_jurisdiction_year.csv
#   3. Notifications by age group         → nndss_counts_by_age_year.csv
#   4. Data caveats                       → nndss_data_caveats.csv
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

outdir <- "data/processed"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Helper: strip commas from numbers
parse_count <- function(x) suppressWarnings(as.numeric(str_remove_all(x, ",")))

# =============================================================================
# 1. National annual counts by disease
# =============================================================================
cat("=== 1. National annual counts by disease ===\n")

df <- suppressMessages(
  read_excel("data/raw/powerbi/nndss_all_diseases_annual.xlsx",
             sheet = 1, col_names = FALSE, col_types = "text",
             .name_repair = "minimal")
)

# Row 1: "YearForMatrixDisplay", NA, 2026, 2025, ..., 1991
# Row 2: "DISEASE GROUP", "Disease Name", "Count_Notification", ...
# Rows 3+: data (groups + individual diseases)
header_years <- as.character(unlist(df[1, ]))
year_cols <- which(!is.na(header_years) & grepl("^[0-9]{4}$", header_years))
years <- as.integer(header_years[year_cols])

data_rows <- df[3:nrow(df), ]
# Remove footer rows (NA disease name, "Total", "Applied filters")
col1 <- as.character(unlist(data_rows[, 1]))
col2 <- as.character(unlist(data_rows[, 2]))
keep <- !is.na(col1) & !is.na(col2) &
  col1 != "Total" & !str_detect(col1, "Applied filters")
data_rows <- data_rows[keep, ]

# Filter to individual diseases (not group totals)
col2_clean <- as.character(unlist(data_rows[, 2]))
data_rows <- data_rows[col2_clean != "Total", ]

all_counts <- list()
for (i in seq_len(nrow(data_rows))) {
  group <- as.character(data_rows[[i, 1]])
  disease <- as.character(data_rows[[i, 2]])
  counts <- sapply(year_cols, function(j) as.character(data_rows[[i, j]]))
  all_counts[[i]] <- tibble(
    disease_group = group,
    disease = disease,
    year = years,
    count = parse_count(counts)
  )
}

counts_national <- bind_rows(all_counts)

cat("Diseases:", n_distinct(counts_national$disease), "\n")
cat("Years:", paste(range(counts_national$year), collapse = "-"), "\n")
cat("Rows:", nrow(counts_national), "\n")
cat("Non-NA counts:", sum(!is.na(counts_national$count)), "\n\n")

write_csv(counts_national, file.path(outdir, "nndss_counts_national_annual.csv"))
cat("Saved: nndss_counts_national_annual.csv\n\n")

# =============================================================================
# 2. Notifications by jurisdiction (all diseases combined)
# =============================================================================
cat("=== 2. Notifications by jurisdiction ===\n")

df2 <- suppressMessages(
  read_excel("data/raw/powerbi/nndss_notifications_by_jurisdiction.xlsx",
             sheet = 1, col_names = FALSE, col_types = "text",
             .name_repair = "minimal")
)

# Row 1: "Diagnosis Year", 1991, 1992, ..., 2026
# Row 2: "State", "Count_Notification", ...
# Rows 3+: ACT, NSW, NT, QLD, SA, TAS, VIC, WA, Total
header_years2 <- as.character(unlist(df2[1, ]))
year_cols2 <- which(!is.na(header_years2) & grepl("^[0-9]{4}$", header_years2))
years2 <- as.integer(header_years2[year_cols2])

data_rows2 <- df2[3:nrow(df2), ]
col1_2 <- as.character(unlist(data_rows2[, 1]))
keep2 <- !is.na(col1_2) & col1_2 != "Total" & !str_detect(col1_2, "Applied filters")
data_rows2 <- data_rows2[keep2, ]

juris_list <- list()
for (i in seq_len(nrow(data_rows2))) {
  state <- as.character(data_rows2[[i, 1]])
  counts <- sapply(year_cols2, function(j) as.character(data_rows2[[i, j]]))
  juris_list[[i]] <- tibble(
    state = state,
    year = years2,
    total_notifications = parse_count(counts)
  )
}

counts_jurisdiction <- bind_rows(juris_list)

cat("States:", paste(sort(unique(counts_jurisdiction$state)), collapse = ", "), "\n")
cat("Years:", paste(range(counts_jurisdiction$year), collapse = "-"), "\n")
cat("Rows:", nrow(counts_jurisdiction), "\n\n")

write_csv(counts_jurisdiction, file.path(outdir, "nndss_counts_by_jurisdiction_year.csv"))
cat("Saved: nndss_counts_by_jurisdiction_year.csv\n\n")

# =============================================================================
# 3. Notifications by age group (all diseases combined)
# =============================================================================
cat("=== 3. Notifications by age group ===\n")

df3 <- suppressMessages(
  read_excel("data/raw/powerbi/nndss_notifications_by_age_group.xlsx",
             sheet = 1, col_names = FALSE, col_types = "text",
             .name_repair = "minimal")
)

# Row 1: "Age Group", "00-04", "05-09", ..., "85+"
# Row 2: "Diagnosis Year", "Count_Notification", ...
# Rows 3+: year rows (2026, 2025, ..., 1991)
header_ages <- as.character(unlist(df3[1, ]))
age_cols <- which(!is.na(header_ages) & header_ages != "Age Group")
age_groups <- header_ages[age_cols]

data_rows3 <- df3[3:nrow(df3), ]
col1_3 <- as.character(unlist(data_rows3[, 1]))
keep3 <- !is.na(col1_3) & grepl("^[0-9]{4}$", col1_3)
data_rows3 <- data_rows3[keep3, ]

age_list <- list()
for (i in seq_len(nrow(data_rows3))) {
  yr <- as.integer(as.character(data_rows3[[i, 1]]))
  counts <- sapply(age_cols, function(j) as.character(data_rows3[[i, j]]))
  age_list[[i]] <- tibble(
    year = yr,
    age_group = age_groups,
    count = parse_count(counts)
  )
}

counts_age <- bind_rows(age_list)

cat("Age groups:", paste(sort(unique(counts_age$age_group)), collapse = ", "), "\n")
cat("Years:", paste(range(counts_age$year), collapse = "-"), "\n")
cat("Rows:", nrow(counts_age), "\n\n")

write_csv(counts_age, file.path(outdir, "nndss_counts_by_age_year.csv"))
cat("Saved: nndss_counts_by_age_year.csv\n\n")

# =============================================================================
# 4. Data caveats
# =============================================================================
cat("=== 4. Data caveats ===\n")

df4 <- suppressMessages(
  read_excel("data/raw/powerbi/nndss_data_caveats.xlsx",
             sheet = 1, col_names = FALSE, col_types = "text",
             .name_repair = "minimal")
)

# Row 1: header (DiseaseName, Jurisdiction, Notes)
# Row 2+: data
caveats <- df4[2:nrow(df4), ]
names(caveats) <- c("disease", "jurisdiction", "notes")
caveats <- caveats %>%
  filter(!is.na(disease), !str_detect(disease, "Applied filters"))

cat("Caveat entries:", nrow(caveats), "\n")
cat("Diseases with caveats:", n_distinct(caveats$disease), "\n\n")

# Show "nationally notifiable from" entries
nat_from <- caveats %>%
  filter(jurisdiction == "AUS", str_detect(notes, "(?i)nationally notifiable from")) %>%
  mutate(year_from = as.integer(str_extract(notes, "[0-9]{4}"))) %>%
  arrange(year_from)

cat("--- Diseases with late national notification start ---\n")
print(nat_from %>% select(disease, year_from, notes), n = 30)

write_csv(caveats, file.path(outdir, "nndss_data_caveats.csv"))
cat("\nSaved: nndss_data_caveats.csv\n")
