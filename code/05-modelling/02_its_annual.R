#!/usr/bin/env Rscript
# =============================================================================
# 02_its_annual.R
#
# Annual Interrupted Time Series for 47 diseases:
#   Part A: Individual disease ITS (44, excl. 3 ultra-low)
#           Simplified 2-segment model with Newey-West HAC SEs (n=11)
#   Part B: Pooled mixed-effects model (primary manuscript analysis)
#           Random intercept + random COVID effect by disease
#   Part C: 5 sensitivity analyses
#
# Inputs:
#   data/processed/annual_with_baselines.csv
#   data/processed/baselines_annual_national.csv
#   data/processed/disease_metadata.csv
#
# Outputs:
#   data/processed/its_annual_individual.csv  — 44 diseases × coefficients
#   data/processed/its_annual_pooled.csv      — pooled model coefficients + RE
#   data/processed/its_sensitivity.csv        — sensitivity analysis comparison
#   output/tables/its_results_table.csv       — manuscript Table 3
#
# Reference: Bernal et al. 2017, Int J Epidemiol (ITS with short series)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lme4)
  library(sandwich)
  library(lmtest)
  library(car)
  library(broom)
  library(broom.mixed)
})

outdir <- "data/processed"
tabdir <- "output/tables"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Phase 5, Script 2: Annual ITS Models ===\n\n")

# =============================================================================
# 1. Load and prepare data
# =============================================================================

abl       <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                      show_col_types = FALSE)
baselines <- read_csv(file.path(outdir, "baselines_annual_national.csv"),
                      show_col_types = FALSE)
metadata  <- read_csv(file.path(outdir, "disease_metadata.csv"),
                      show_col_types = FALSE)

# National level, analysis-suitable, 2015-2025
national <- abl %>%
  filter(level == "national", analysis_suitable == TRUE,
         year >= 2015, year <= 2025) %>%
  arrange(disease, year)

# Merge transmission mode info from baselines
disease_info <- baselines %>%
  select(disease, count_tier) %>%
  left_join(
    metadata %>% select(disease, transmission_mode, border_sensitive,
                        vaccine_preventable),
    by = "disease"
  )

# count_tier, transmission_mode, border_sensitive, vaccine_preventable
# already present in national (from abl)

n_diseases <- n_distinct(national$disease)
cat("Diseases in dataset:", n_diseases, "\n")

# Identify ultra-low count diseases (baseline <50/yr) for quasi-Poisson
ultra_low <- baselines %>%
  filter(baseline_mean_count < 50) %>%
  pull(disease)

# Identify true ultra-low (tetanus, congenital syphilis, diphtheria) to exclude
exclude_diseases <- c("Tetanus", "Congenital syphilis", "Diphtheria")
analysis_diseases <- setdiff(unique(national$disease), exclude_diseases)
low_count_diseases <- setdiff(
  baselines %>% filter(baseline_mean_count < 50) %>% pull(disease),
  exclude_diseases
)

cat("Analysis diseases:", length(analysis_diseases),
    "\nLow-count (quasi-Poisson):", length(low_count_diseases),
    "\nExcluded ultra-low:", length(exclude_diseases), "\n\n")

# Create ITS variables
national <- national %>%
  mutate(
    year_centered     = year - 2017,
    covid             = as.integer(year >= 2020),
    time_since_covid  = pmax(0L, as.integer(year - 2019))
  )

# =============================================================================
# Part A: Individual disease ITS
# =============================================================================

cat("--- Part A: Individual Disease ITS ---\n\n")

individual_results <- list()

for (dis in analysis_diseases) {
  d <- national %>% filter(disease == dis)

  if (nrow(d) < 5) {
    cat("  Skipping", dis, "— insufficient data (n =", nrow(d), ")\n")
    next
  }

  use_poisson <- dis %in% low_count_diseases

  if (use_poisson) {
    # Quasi-Poisson GLM on counts with offset
    # Use count column and a population offset
    d_model <- d %>%
      mutate(
        count_val = as.numeric(count),
        log_pop   = log(1e5)  # rates already per 100k; offset normalises
      )

    tryCatch({
      fit <- glm(count_val ~ year_centered + covid + time_since_covid,
                 family = quasipoisson(link = "log"),
                 data = d_model)

      coefs <- coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))
      tidy_coefs <- tidy(coefs, conf.int = TRUE)

      tidy_coefs <- tidy_coefs %>%
        mutate(
          disease    = dis,
          model_type = "quasipoisson",
          pct_change = (exp(estimate) - 1) * 100,
          pct_change_lo = (exp(conf.low) - 1) * 100,
          pct_change_hi = (exp(conf.high) - 1) * 100,
          n_obs      = nrow(d_model),
          r_squared  = NA_real_,
          dispersion = summary(fit)$dispersion
        )

      individual_results[[dis]] <- tidy_coefs
    }, error = function(e) {
      cat("  ERROR for", dis, ":", conditionMessage(e), "\n")
    })

  } else {
    # OLS on rate_per_100k with Newey-West HAC SEs (lag=1)
    tryCatch({
      fit <- lm(rate_per_100k ~ year_centered + covid + time_since_covid,
                data = d)

      # Newey-West heteroscedasticity and autocorrelation consistent SEs
      nw_vcov <- NeweyWest(fit, lag = 1, prewhite = FALSE)
      coefs   <- coeftest(fit, vcov. = nw_vcov)
      tidy_coefs <- tidy(coefs, conf.int = TRUE)

      tidy_coefs <- tidy_coefs %>%
        mutate(
          disease    = dis,
          model_type = "ols_neweywest",
          # For OLS on rates, pct change relative to baseline mean
          baseline_mean = baselines %>%
            filter(disease == dis) %>%
            pull(baseline_mean_rate),
          pct_change = estimate / baseline_mean * 100,
          pct_change_lo = conf.low / baseline_mean * 100,
          pct_change_hi = conf.high / baseline_mean * 100,
          n_obs      = nrow(d),
          r_squared  = summary(fit)$r.squared,
          dispersion = NA_real_
        )

      individual_results[[dis]] <- tidy_coefs
    }, error = function(e) {
      cat("  ERROR for", dis, ":", conditionMessage(e), "\n")
    })
  }
}

its_individual <- bind_rows(individual_results) %>%
  select(disease, term, estimate, std.error, statistic, p.value,
         conf.low, conf.high, pct_change, pct_change_lo, pct_change_hi,
         model_type, n_obs, r_squared, dispersion)

write_csv(its_individual, file.path(outdir, "its_annual_individual.csv"))
cat("Saved individual ITS:", nrow(its_individual), "rows for",
    n_distinct(its_individual$disease), "diseases\n\n")

# Print COVID coefficient summary
cat("--- COVID coefficient summary (top 10 most disrupted) ---\n")
its_individual %>%
  filter(term == "covid") %>%
  arrange(estimate) %>%
  head(10) %>%
  select(disease, estimate, pct_change, p.value, model_type) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  as.data.frame() %>%
  print()

# Count significant results
n_sig <- its_individual %>%
  filter(term == "covid", p.value < 0.05) %>%
  nrow()
cat("\nDiseases with significant COVID coefficient (p<0.05):", n_sig,
    "of", n_distinct(its_individual$disease), "\n\n")

# =============================================================================
# Part B: Pooled Mixed-Effects Model
# =============================================================================

cat("--- Part B: Pooled Mixed-Effects Model ---\n\n")

# Prepare data: standardise rates within each disease (z-score vs baseline)
# Note: baseline_mean_rate, baseline_sd_rate, transmission_mode,
# border_sensitive, vaccine_preventable already in national (from abl)
pooled_data <- national %>%
  filter(disease %in% analysis_diseases) %>%
  mutate(
    # Z-score relative to baseline
    rate_z = (rate_per_100k - baseline_mean_rate) /
             ifelse(baseline_sd_rate > 0, baseline_sd_rate, 1),
    transmission_mode = factor(transmission_mode,
                               levels = c("respiratory", "enteric",
                                           "vector_borne", "zoonotic",
                                           "STI", "blood_borne")),
    border_sensitive   = factor(border_sensitive),
    vaccine_preventable = factor(vaccine_preventable)
  ) %>%
  filter(!is.na(rate_z), is.finite(rate_z))

cat("Pooled data:", nrow(pooled_data), "rows,",
    n_distinct(pooled_data$disease), "diseases\n")

# --- Model 1: Transmission mode interaction ---
cat("\nFitting pooled model with transmission_mode interaction...\n")
pooled_tm <- tryCatch({
  lmer(rate_z ~ year_centered + covid * transmission_mode +
         time_since_covid * transmission_mode +
         (1 + covid | disease),
       data = pooled_data, REML = FALSE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}, error = function(e) {
  cat("  Full model failed:", conditionMessage(e), "\n")
  cat("  Trying simplified random effects...\n")
  lmer(rate_z ~ year_centered + covid * transmission_mode +
         time_since_covid * transmission_mode +
         (1 | disease),
       data = pooled_data, REML = FALSE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
})

cat("  AIC:", round(AIC(pooled_tm), 1), ", BIC:", round(BIC(pooled_tm), 1), "\n")

# --- Model 2: Border-sensitive interaction ---
cat("Fitting pooled model with border_sensitive interaction...\n")
pooled_bs <- tryCatch({
  lmer(rate_z ~ year_centered + covid * border_sensitive +
         time_since_covid * border_sensitive +
         (1 + covid | disease),
       data = pooled_data, REML = FALSE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}, error = function(e) {
  lmer(rate_z ~ year_centered + covid * border_sensitive +
         time_since_covid * border_sensitive +
         (1 | disease),
       data = pooled_data, REML = FALSE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
})

cat("  AIC:", round(AIC(pooled_bs), 1), ", BIC:", round(BIC(pooled_bs), 1), "\n")

# --- Model 3: Vaccine-preventable interaction ---
cat("Fitting pooled model with vaccine_preventable interaction...\n")
pooled_vp <- tryCatch({
  lmer(rate_z ~ year_centered + covid * vaccine_preventable +
         time_since_covid * vaccine_preventable +
         (1 + covid | disease),
       data = pooled_data, REML = FALSE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}, error = function(e) {
  lmer(rate_z ~ year_centered + covid * vaccine_preventable +
         time_since_covid * vaccine_preventable +
         (1 | disease),
       data = pooled_data, REML = FALSE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
})

cat("  AIC:", round(AIC(pooled_vp), 1), ", BIC:", round(BIC(pooled_vp), 1), "\n\n")

# Collect pooled results
extract_pooled <- function(model, model_name) {
  fe <- tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      effect_type = "fixed",
      # Approximate p-value from z-distribution (conservative for lmer)
      p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE)
    )
  re <- tidy(model, effects = "ran_pars") %>%
    mutate(effect_type = "random", conf.low = NA, conf.high = NA,
           p.value = NA_real_, statistic = NA_real_)
  bind_rows(fe, re) %>%
    mutate(model_name = model_name,
           aic = AIC(model), bic = BIC(model),
           n_obs = nrow(model@frame),
           n_groups = ngrps(model))
}

pooled_results <- bind_rows(
  extract_pooled(pooled_tm, "transmission_mode"),
  extract_pooled(pooled_bs, "border_sensitive"),
  extract_pooled(pooled_vp, "vaccine_preventable")
)

# Random effects by disease (for primary model)
re_by_disease <- ranef(pooled_tm)$disease %>%
  tibble::rownames_to_column("disease") %>%
  rename(re_intercept = `(Intercept)`)

# Check if covid random slope exists
if ("covid" %in% names(ranef(pooled_tm)$disease)) {
  re_by_disease <- re_by_disease %>% rename(re_covid = covid)
} else {
  re_by_disease <- re_by_disease %>% mutate(re_covid = NA_real_)
}

pooled_all <- list(
  coefficients = pooled_results,
  random_effects = re_by_disease
)

# Save pooled results
write_csv(pooled_results, file.path(outdir, "its_annual_pooled.csv"))
write_csv(re_by_disease, file.path(outdir, "its_annual_random_effects.csv"))
cat("Saved pooled model results\n")

# Print key interaction terms
cat("\n--- Key interaction terms ---\n")
pooled_results %>%
  filter(grepl("covid:", term) | grepl(":covid", term),
         effect_type == "fixed") %>%
  select(model_name, term, estimate, std.error, p.value) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  as.data.frame() %>%
  print()

# =============================================================================
# Part C: Sensitivity Analyses
# =============================================================================

cat("\n\n--- Part C: Sensitivity Analyses ---\n\n")

sensitivity_results <- list()

# Helper: fit individual ITS and return COVID coefficient
fit_its_covid <- function(data, label) {
  results <- list()
  for (dis in unique(data$disease)) {
    d <- data %>% filter(disease == dis)
    if (nrow(d) < 5 || dis %in% exclude_diseases) next

    tryCatch({
      fit <- lm(rate_per_100k ~ year_centered + covid + time_since_covid,
                data = d)
      nw_vcov <- NeweyWest(fit, lag = 1, prewhite = FALSE)
      coefs   <- coeftest(fit, vcov. = nw_vcov)
      tidy_c  <- tidy(coefs, conf.int = TRUE) %>%
        filter(term == "covid") %>%
        mutate(disease = dis, sensitivity = label)
      results[[dis]] <- tidy_c
    }, error = function(e) NULL)
  }
  bind_rows(results)
}

# S1: Trend-adjusted baselines (use O/E ratio as outcome for 13 trended diseases)
cat("S1: Trend-adjusted baselines...\n")
trended <- baselines %>% filter(trend_significant == TRUE) %>% pull(disease)

s1_data <- national %>%
  filter(disease %in% analysis_diseases) %>%
  mutate(
    oe_outcome = ifelse(disease %in% trended & !is.na(oe_ratio),
                        oe_ratio, rate_per_100k)
  )

s1_results <- list()
for (dis in analysis_diseases) {
  d <- s1_data %>% filter(disease == dis)
  if (nrow(d) < 5) next
  outcome_col <- if (dis %in% trended) "oe_ratio" else "rate_per_100k"
  tryCatch({
    d$outcome <- d[[outcome_col]]
    if (all(is.na(d$outcome))) next
    fit <- lm(outcome ~ year_centered + covid + time_since_covid, data = d)
    nw_vcov <- NeweyWest(fit, lag = 1, prewhite = FALSE)
    coefs <- coeftest(fit, vcov. = nw_vcov)
    s1_results[[dis]] <- tidy(coefs, conf.int = TRUE) %>%
      filter(term == "covid") %>%
      mutate(disease = dis, sensitivity = "S1_trend_adjusted")
  }, error = function(e) NULL)
}
sensitivity_results[["S1"]] <- bind_rows(s1_results)
cat("  Fitted", length(s1_results), "models\n")

# S2: Exclude 2025 (use 2015-2024 only)
cat("S2: Exclude 2025...\n")
s2_data <- national %>%
  filter(year <= 2024, disease %in% analysis_diseases)
sensitivity_results[["S2"]] <- fit_its_covid(s2_data, "S2_exclude_2025")
cat("  Fitted", n_distinct(sensitivity_results[["S2"]]$disease), "models\n")

# S3: Alternative period — COVID ends 2022 (covid = 1 if year >= 2020)
cat("S3: COVID includes 2022...\n")
s3_data <- national %>%
  filter(disease %in% analysis_diseases) %>%
  mutate(
    covid             = as.integer(year >= 2020),
    time_since_covid  = pmax(0L, as.integer(year - 2019))
  )
# This is actually the same coding since we already code covid = year >= 2020
# What changes is the interpretation: post-COVID starts 2023
# The model structure is the same, so S3 = primary for this coding
# Instead, recode: covid = 1 only for 2020-2022, post starts 2023
s3_data <- national %>%
  filter(disease %in% analysis_diseases) %>%
  mutate(
    covid             = as.integer(year >= 2020 & year <= 2022),
    time_since_covid  = ifelse(year >= 2020, pmin(as.integer(year - 2019), 3L), 0L)
  )
sensitivity_results[["S3"]] <- fit_its_covid(s3_data, "S3_covid_includes_2022")
cat("  Fitted", n_distinct(sensitivity_results[["S3"]]$disease), "models\n")

# S4: Poisson vs OLS comparison (for high-count diseases only)
cat("S4: Poisson vs OLS comparison...\n")
high_count <- setdiff(analysis_diseases, low_count_diseases)
s4_results <- list()
for (dis in high_count) {
  d <- national %>% filter(disease == dis)
  if (nrow(d) < 5) next
  tryCatch({
    # Poisson on counts
    d$count_val <- as.numeric(d$count)
    fit_pois <- glm(count_val ~ year_centered + covid + time_since_covid,
                    family = quasipoisson(link = "log"), data = d)
    coefs_pois <- coeftest(fit_pois, vcov. = vcovHC(fit_pois, type = "HC1"))
    s4_results[[dis]] <- tidy(coefs_pois, conf.int = TRUE) %>%
      filter(term == "covid") %>%
      mutate(disease = dis, sensitivity = "S4_poisson",
             pct_change_pois = (exp(estimate) - 1) * 100)
  }, error = function(e) NULL)
}
sensitivity_results[["S4"]] <- bind_rows(s4_results)
cat("  Fitted", length(s4_results), "models\n")

# S5: Exclude high-CV diseases (CV > 0.5)
cat("S5: Exclude high-CV diseases...\n")
high_cv <- baselines %>% filter(baseline_cv_rate > 0.5) %>% pull(disease)
cat("  High-CV diseases excluded:", paste(high_cv, collapse = ", "), "\n")
s5_diseases <- setdiff(analysis_diseases, high_cv)
s5_data <- national %>% filter(disease %in% s5_diseases)
sensitivity_results[["S5"]] <- fit_its_covid(s5_data, "S5_exclude_high_cv")
cat("  Fitted", n_distinct(sensitivity_results[["S5"]]$disease), "models\n")

# Combine sensitivity results
its_sensitivity <- bind_rows(sensitivity_results) %>%
  select(disease, sensitivity, estimate, std.error, p.value, conf.low, conf.high)

# Compare with primary analysis
primary_covid <- its_individual %>%
  filter(term == "covid") %>%
  select(disease, estimate, p.value) %>%
  rename(primary_estimate = estimate, primary_p = p.value)

sensitivity_comparison <- its_sensitivity %>%
  left_join(primary_covid, by = "disease") %>%
  mutate(
    same_direction = sign(estimate) == sign(primary_estimate),
    both_significant = (p.value < 0.05) & (primary_p < 0.05)
  )

write_csv(its_sensitivity, file.path(outdir, "its_sensitivity.csv"))

cat("\n--- Sensitivity Concordance ---\n")
sensitivity_comparison %>%
  group_by(sensitivity) %>%
  summarise(
    n_diseases = n(),
    same_direction_pct = round(mean(same_direction, na.rm = TRUE) * 100, 1),
    both_sig_pct = round(mean(both_significant, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print()

# =============================================================================
# Manuscript Table 3: ITS Results Summary
# =============================================================================

cat("\nGenerating manuscript Table 3...\n")

# Select representative diseases for table
representative <- c(
  "Influenza (laboratory confirmed)",
  "Pertussis",
  "Invasive pneumococcal disease",
  "Measles",
  "Meningococcal disease (invasive)",
  "Salmonellosis",
  "Campylobacteriosis",
  "Cryptosporidiosis",
  "Dengue virus infection",
  "Ross River virus infection",
  "Chlamydia",
  "Gonorrhoea",
  "Hepatitis C (newly acquired)"
)

table3 <- its_individual %>%
  filter(term == "covid", disease %in% representative) %>%
  left_join(
    disease_info %>% select(disease, transmission_mode, border_sensitive,
                            vaccine_preventable),
    by = "disease"
  ) %>%
  select(disease, transmission_mode, border_sensitive, vaccine_preventable,
         estimate, conf.low, conf.high, p.value, pct_change,
         model_type, r_squared) %>%
  arrange(estimate) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

write_csv(table3, file.path(tabdir, "its_results_table.csv"))
cat("Saved Table 3:", nrow(table3), "representative diseases\n")

cat("\n=== Script 2 complete ===\n")
