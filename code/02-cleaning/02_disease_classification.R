#!/usr/bin/env Rscript
# =============================================================================
# 02_disease_classification.R
#
# Create disease metadata CSV mapping each disease to:
#   - transmission_mode (respiratory, enteric, STI, vector-borne, etc.)
#   - analysis_suitable (has full baseline + COVID + post-COVID)
#   - year_notifiable_from
#   - key caveats
#
# Output: data/processed/disease_metadata.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ---- Build classification table ----

disease_meta <- tribble(
  ~disease, ~transmission_mode, ~subcategory,

  # --- Respiratory ---
  "Influenza (laboratory confirmed)", "respiratory", "viral_respiratory",
  "Pertussis", "respiratory", "vaccine_preventable",
  "Pneumococcal disease (invasive)", "respiratory", "vaccine_preventable",
  "Meningococcal disease (invasive)", "respiratory", "vaccine_preventable",
  "Measles", "respiratory", "vaccine_preventable",
  "Mumps", "respiratory", "vaccine_preventable",
  "Rubella", "respiratory", "vaccine_preventable",
  "Rubella congenital", "respiratory", "vaccine_preventable",
  "Tuberculosis", "respiratory", "mycobacterial",
  "Legionellosis", "respiratory", "environmental_respiratory",
  "Diphtheria", "respiratory", "vaccine_preventable",
  "Haemophilus influenzae type b", "respiratory", "vaccine_preventable",
  "Varicella zoster (chickenpox)", "respiratory", "vaccine_preventable",
  "Varicella zoster (shingles)", "respiratory", "viral_respiratory",
  "Varicella zoster (unspecified)", "respiratory", "viral_respiratory",
  "Respiratory syncytial virus (RSV)", "respiratory", "viral_respiratory",
  "Invasive Group A Streptococcal disease (iGAS)", "respiratory", "bacterial_respiratory",
  "Ornithosis", "respiratory", "zoonotic_respiratory",

  # --- Enteric ---
  "Salmonellosis", "enteric", "bacterial_enteric",
  "Campylobacteriosis", "enteric", "bacterial_enteric",
  "Cryptosporidiosis", "enteric", "parasitic_enteric",
  "Shigellosis", "enteric", "bacterial_enteric",
  "Listeriosis", "enteric", "bacterial_enteric",
  "Hepatitis A", "enteric", "viral_enteric",
  "Hepatitis E", "enteric", "viral_enteric",
  "Typhoid Fever", "enteric", "bacterial_enteric",
  "Paratyphoid", "enteric", "bacterial_enteric",
  "STEC", "enteric", "bacterial_enteric",
  "Rotavirus", "enteric", "viral_enteric",
  "Cholera", "enteric", "bacterial_enteric",
  "Botulism", "enteric", "bacterial_enteric",
  "Haemolytic uraemic syndrome (HUS)", "enteric", "bacterial_enteric",
  "Vibrio parahaemolyticus infection", "enteric", "bacterial_enteric",

  # --- STI ---
  "Chlamydial infection", "STI", "bacterial_STI",
  "Gonococcal infection", "STI", "bacterial_STI",
  "Syphilis < 2 years", "STI", "bacterial_STI",
  "Syphilis > 2 years or unspecified duration", "STI", "bacterial_STI",
  "Syphilis congenital", "STI", "bacterial_STI",
  "Donovanosis", "STI", "bacterial_STI",

  # --- Blood-borne ---
  "Hepatitis B (newly acquired)", "blood_borne", "viral_bloodborne",
  "Hepatitis B (unspecified)", "blood_borne", "viral_bloodborne",
  "Hepatitis C (newly acquired)", "blood_borne", "viral_bloodborne",
  "Hepatitis C (unspecified)", "blood_borne", "viral_bloodborne",
  "Hepatitis D", "blood_borne", "viral_bloodborne",
  "Hepatitis (NEC)", "blood_borne", "viral_bloodborne",

  # --- Vector-borne ---
  "Ross River virus infection", "vector_borne", "mosquito_borne",
  "Barmah Forest virus infection", "vector_borne", "mosquito_borne",
  "Dengue virus infection", "vector_borne", "mosquito_borne",
  "Malaria", "vector_borne", "mosquito_borne",
  "Murray Valley encephalitis virus infection", "vector_borne", "mosquito_borne",
  "Japanese encephalitis virus infection", "vector_borne", "mosquito_borne",
  "West Nile/Kunjin virus infection", "vector_borne", "mosquito_borne",
  "Chikungunya virus infection", "vector_borne", "mosquito_borne",
  "Flavivirus infection (unspecified)", "vector_borne", "mosquito_borne",

  # --- Zoonotic / environmental ---
  "Q fever", "zoonotic", "zoonotic_bacterial",
  "Leptospirosis", "zoonotic", "zoonotic_bacterial",
  "Brucellosis", "zoonotic", "zoonotic_bacterial",
  "Anthrax", "zoonotic", "zoonotic_bacterial",
  "Australian bat lyssavirus infection", "zoonotic", "zoonotic_viral",
  "Tularaemia", "zoonotic", "zoonotic_bacterial",
  "Avian influenza in humans (AIH)", "zoonotic", "zoonotic_viral",
  "Leprosy", "zoonotic", "mycobacterial",
  "Tetanus", "zoonotic", "environmental",
  "Poliovirus infection", "other", "vaccine_preventable",

  # --- Other / pandemic ---
  "COVID-19", "other", "pandemic",
  "Mpox", "other", "emerging"
)

# ---- Add analysis suitability from QC results ----

rates <- read_csv("data/processed/nndss_rates_by_state_year.csv", show_col_types = FALSE)

suitability <- rates %>%
  filter(state == "AUS") %>%
  group_by(disease) %>%
  summarise(
    year_data_from = min(year[!is.na(rate_per_100k)], na.rm = TRUE),
    year_data_to = max(year[!is.na(rate_per_100k)], na.rm = TRUE),
    baseline_years_available = sum(year >= 2015 & year <= 2019 & !is.na(rate_per_100k)),
    has_full_baseline = baseline_years_available == 5,
    has_covid_data = sum(year >= 2020 & year <= 2021 & !is.na(rate_per_100k)) == 2,
    has_post_data = sum(year >= 2022 & year <= 2025 & !is.na(rate_per_100k)) >= 2,
    analysis_suitable = has_full_baseline & has_covid_data & has_post_data,
    rate_2019 = rate_per_100k[year == 2019][1],
    rate_2020 = rate_per_100k[year == 2020][1],
    .groups = "drop"
  ) %>%
  mutate(
    year_data_from = ifelse(is.infinite(year_data_from), NA_real_, year_data_from),
    year_data_to = ifelse(is.infinite(year_data_to), NA_real_, year_data_to)
  )

# ---- Add caveats ----

caveats <- read_csv("data/processed/nndss_data_caveats.csv", show_col_types = FALSE)

caveat_summary <- caveats %>%
  filter(jurisdiction == "AUS") %>%
  group_by(disease) %>%
  summarise(national_caveats = paste(notes, collapse = "; "), .groups = "drop")

# ---- Merge everything ----

result <- disease_meta %>%
  left_join(suitability, by = "disease") %>%
  left_join(caveat_summary, by = "disease") %>%
  mutate(
    # Border-sensitive: acquired overseas vs locally
    border_sensitive = disease %in% c(
      "Dengue virus infection", "Malaria", "Chikungunya virus infection",
      "Measles", "Typhoid Fever", "Paratyphoid", "Hepatitis A",
      "Cholera", "Shigellosis", "Hepatitis E", "Leprosy"
    ),
    # Vaccine-preventable
    vaccine_preventable = subcategory == "vaccine_preventable" |
      disease %in% c("Rotavirus", "Hepatitis A", "Hepatitis B (newly acquired)",
                      "Influenza (laboratory confirmed)"),
    # Exclusion reason for non-suitable diseases
    exclusion_reason = case_when(
      analysis_suitable ~ NA_character_,
      disease == "COVID-19" ~ "Pandemic pathogen (not a disrupted disease)",
      disease == "Mpox" ~ "Emerged 2022 (no pre-COVID baseline)",
      disease == "Respiratory syncytial virus (RSV)" ~ "Notifiable from 2020 only",
      disease == "Invasive Group A Streptococcal disease (iGAS)" ~ "Notifiable from 2021 only",
      disease == "Vibrio parahaemolyticus infection" ~ "Notifiable from 2025 only",
      !has_full_baseline ~ "Insufficient pre-COVID baseline data",
      !has_covid_data ~ "Missing COVID-period data",
      !has_post_data ~ "Missing post-COVID data",
      TRUE ~ "Other"
    )
  )

# ---- Check for unmatched diseases ----

rate_diseases <- unique(rates$disease)
unmatched <- setdiff(rate_diseases, result$disease)
if (length(unmatched) > 0) {
  cat("WARNING: Diseases in rate data but NOT classified:\n")
  for (d in unmatched) cat("  -", d, "\n")
}

# ---- Summary ----

cat("\n=== DISEASE CLASSIFICATION SUMMARY ===\n\n")
cat("Total diseases classified:", nrow(result), "\n")
cat("Analysis-suitable:", sum(result$analysis_suitable, na.rm = TRUE), "\n\n")

cat("By transmission mode:\n")
result %>%
  group_by(transmission_mode) %>%
  summarise(
    total = n(),
    suitable = sum(analysis_suitable, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(suitable)) %>%
  print(n = 10)

cat("\nBorder-sensitive diseases:", sum(result$border_sensitive), "\n")
cat("Vaccine-preventable diseases:", sum(result$vaccine_preventable), "\n")

cat("\n--- Excluded diseases and reasons ---\n")
result %>%
  filter(!analysis_suitable | is.na(analysis_suitable)) %>%
  select(disease, transmission_mode, exclusion_reason) %>%
  print(n = 25)

# ---- Save ----

outfile <- "data/processed/disease_metadata.csv"
write_csv(result, outfile)
cat("\nSaved to:", outfile, "\n")
