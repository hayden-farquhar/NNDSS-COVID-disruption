#!/usr/bin/env Rscript
# =============================================================================
# 03_predictors_stringency.R
#
# Part A: Nested predictor regression hierarchy (AIC/BIC comparison)
#         Tests what predicts disruption magnitude: transmission mode,
#         border sensitivity, vaccine-preventability, baseline rate
# Part B: State-level stringency correlation (OxCGRT data)
#         Spearman correlation + mixed-effects model
# Part C: Recovery projections for 12 still-suppressed diseases
#
# Inputs:
#   data/processed/disruption_metrics_national.csv
#   data/processed/disruption_by_state.csv
#   data/processed/baselines_annual_national.csv
#   data/processed/recovery_trajectories.csv
#   data/processed/annual_with_baselines.csv
#   data/raw/contextual/oxcgrt_subnational.csv
#
# Outputs:
#   data/processed/stringency_by_state.csv           — parsed stringency data
#   data/processed/predictor_regression_results.csv   — nested model comparison
#   data/processed/stringency_correlation.csv         — state-level results
#   data/processed/recovery_projections.csv           — projected recovery years
#   output/tables/predictor_regression_table.csv      — manuscript table
#   output/tables/stringency_table.csv                — stringency × disruption
#   output/figures/fig7_stringency_scatter.png/.pdf   — scatter with Spearman rho
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(lme4)
  library(car)
  library(broom)
  library(broom.mixed)
})

outdir <- "data/processed"
tabdir <- "output/tables"
figdir <- "output/figures"
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Phase 5, Script 3: Predictors, Stringency & Recovery ===\n\n")

# =============================================================================
# Load data
# =============================================================================

disruption  <- read_csv(file.path(outdir, "disruption_metrics_national.csv"),
                        show_col_types = FALSE)
state_disr  <- read_csv(file.path(outdir, "disruption_by_state.csv"),
                        show_col_types = FALSE)
baselines   <- read_csv(file.path(outdir, "baselines_annual_national.csv"),
                        show_col_types = FALSE)
recovery    <- read_csv(file.path(outdir, "recovery_trajectories.csv"),
                        show_col_types = FALSE)
abl         <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                        show_col_types = FALSE)

cat("Loaded disruption metrics:", nrow(disruption), "diseases\n")
cat("Loaded state disruption:", nrow(state_disr), "rows\n")
cat("Loaded recovery trajectories:", nrow(recovery), "diseases\n\n")

# =============================================================================
# Part A: Predictor Regression Hierarchy
# =============================================================================

cat("--- Part A: Predictor Regression ---\n\n")

# Prepare predictor data
# Note: baseline_mean_rate, transmission_mode, etc. already in disruption
pred_data <- disruption %>%
  filter(!is.na(covid_oe_rate), is.finite(covid_oe_rate)) %>%
  mutate(
    transmission_mode   = factor(transmission_mode,
                                 levels = c("respiratory", "enteric",
                                            "vector_borne", "zoonotic",
                                            "STI", "blood_borne")),
    border_sensitive    = factor(border_sensitive),
    vaccine_preventable = factor(vaccine_preventable),
    log_baseline_rate   = log1p(baseline_mean_rate)
  ) %>%
  filter(!is.na(transmission_mode))

cat("Predictor dataset:", nrow(pred_data), "diseases\n\n")

# --- Nested model hierarchy ---
M1 <- lm(covid_oe_rate ~ transmission_mode, data = pred_data)
M2 <- lm(covid_oe_rate ~ transmission_mode + border_sensitive, data = pred_data)
M3 <- lm(covid_oe_rate ~ transmission_mode + border_sensitive + vaccine_preventable,
         data = pred_data)
M4 <- lm(covid_oe_rate ~ transmission_mode + border_sensitive + vaccine_preventable +
           log_baseline_rate, data = pred_data)
M5 <- lm(covid_oe_rate ~ border_sensitive + vaccine_preventable + log_baseline_rate,
         data = pred_data)

models <- list(M1 = M1, M2 = M2, M3 = M3, M4 = M4, M5 = M5)

# Model comparison table
model_comparison <- tibble(
  model  = names(models),
  formula = sapply(models, function(m) deparse(formula(m), width.cutoff = 200)),
  r_squared     = sapply(models, function(m) summary(m)$r.squared),
  adj_r_squared = sapply(models, function(m) summary(m)$adj.r.squared),
  aic    = sapply(models, AIC),
  bic    = sapply(models, BIC),
  df     = sapply(models, function(m) m$rank),
  f_stat = sapply(models, function(m) {
    s <- summary(m)
    if (!is.null(s$fstatistic)) s$fstatistic[1] else NA
  }),
  f_pvalue = sapply(models, function(m) {
    s <- summary(m)
    if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA
  })
)

cat("--- Nested Model Comparison ---\n")
model_comparison %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  select(model, r_squared, adj_r_squared, aic, bic) %>%
  as.data.frame() %>%
  print()

# VIF check on M4
cat("\n--- VIF for M4 (full model) ---\n")
vif_m4 <- tryCatch({
  vif(M4)
}, error = function(e) {
  cat("  VIF computation failed:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(vif_m4)) {
  # For factors, GVIF^(1/2*Df) is the relevant statistic
  if (is.matrix(vif_m4)) {
    print(round(vif_m4, 3))
    max_gvif <- max(vif_m4[, "GVIF^(1/(2*Df))"])
    cat("\nMax GVIF^(1/2Df):", round(max_gvif, 3),
        ifelse(max_gvif > 2.236, " — CONCERN (>sqrt(5))", " — OK"), "\n")
  } else {
    print(round(vif_m4, 3))
    cat("\nMax VIF:", round(max(vif_m4), 3),
        ifelse(max(vif_m4) > 5, " — CONCERN (>5)", " — OK"), "\n")
  }
}

# Best model coefficients
best_model_name <- model_comparison %>%
  filter(bic == min(bic)) %>%
  pull(model)
best_model <- models[[best_model_name]]
cat("\nBest model by BIC:", best_model_name, "\n\n")

best_coefs <- tidy(best_model, conf.int = TRUE) %>%
  mutate(model = best_model_name)

cat("--- Best model coefficients ---\n")
best_coefs %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  as.data.frame() %>%
  print()

# --- Recovery regression: more disrupted → greater overshoot? ---
cat("\n--- Recovery regression (immunity debt hypothesis) ---\n")
recovery_reg_data <- disruption %>%
  filter(!is.na(post_oe_rate), is.finite(post_oe_rate),
         !is.na(covid_oe_rate), is.finite(covid_oe_rate)) %>%
  mutate(
    border_sensitive    = factor(border_sensitive),
    vaccine_preventable = factor(vaccine_preventable)
  )

recovery_model <- lm(post_oe_rate ~ covid_oe_rate + border_sensitive + vaccine_preventable,
                     data = recovery_reg_data)

cat("Recovery model R²:", round(summary(recovery_model)$r.squared, 4), "\n")
tidy(recovery_model, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  as.data.frame() %>%
  print()

# Save Part A results
all_model_coefs <- bind_rows(
  lapply(names(models), function(nm) {
    tidy(models[[nm]], conf.int = TRUE) %>%
      mutate(model = nm)
  })
)

write_csv(model_comparison, file.path(outdir, "predictor_regression_results.csv"))
write_csv(all_model_coefs, file.path(tabdir, "predictor_regression_table.csv"))
cat("\nSaved predictor regression results\n\n")

# =============================================================================
# Part B: State Stringency Correlation
# =============================================================================

cat("--- Part B: State Stringency Correlation ---\n\n")

# Parse OxCGRT subnational data
cat("Parsing OxCGRT subnational data...\n")
oxcgrt <- read_csv("data/raw/contextual/oxcgrt_subnational.csv",
                   show_col_types = FALSE)

# Filter to Australian states
aus_states <- c("Australian Capital Territory", "New South Wales",
                "Northern Territory", "Queensland", "South Australia",
                "Tasmania", "Victoria", "Western Australia")

oxcgrt_aus <- oxcgrt %>%
  filter(CountryName == "Australia",
         RegionName %in% aus_states,
         Jurisdiction == "STATE_TOTAL") %>%
  mutate(
    date = as.Date(as.character(Date), format = "%Y%m%d"),
    state = case_when(
      RegionName == "Australian Capital Territory" ~ "ACT",
      RegionName == "New South Wales" ~ "NSW",
      RegionName == "Northern Territory" ~ "NT",
      RegionName == "Queensland" ~ "QLD",
      RegionName == "South Australia" ~ "SA",
      RegionName == "Tasmania" ~ "TAS",
      RegionName == "Victoria" ~ "VIC",
      RegionName == "Western Australia" ~ "WA",
      TRUE ~ RegionName
    )
  ) %>%
  select(state, date, StringencyIndex_Average)

cat("OxCGRT Australian records:", nrow(oxcgrt_aus), "\n")

# Filter to COVID-acute period (Mar 2020 – Dec 2021)
covid_period <- oxcgrt_aus %>%
  filter(date >= as.Date("2020-03-01"), date <= as.Date("2021-12-31")) %>%
  filter(!is.na(StringencyIndex_Average))

cat("COVID-acute period records:", nrow(covid_period), "\n")

# Compute per-state stringency metrics
stringency_by_state <- covid_period %>%
  group_by(state) %>%
  summarise(
    mean_stringency      = mean(StringencyIndex_Average, na.rm = TRUE),
    median_stringency    = median(StringencyIndex_Average, na.rm = TRUE),
    max_stringency       = max(StringencyIndex_Average, na.rm = TRUE),
    cumulative_stringency = sum(StringencyIndex_Average, na.rm = TRUE),
    days_above_70        = sum(StringencyIndex_Average > 70, na.rm = TRUE),
    days_total           = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_stringency))

cat("\n--- State Stringency Summary ---\n")
stringency_by_state %>%
  mutate(across(where(is.numeric), ~ round(., 1))) %>%
  as.data.frame() %>%
  print()

write_csv(stringency_by_state, file.path(outdir, "stringency_by_state.csv"))

# --- Merge stringency with state-level disruption ---
# Compute median COVID O/E per state across diseases
state_median_oe <- state_disr %>%
  filter(!is.na(covid_oe), is.finite(covid_oe)) %>%
  group_by(state) %>%
  summarise(
    median_covid_oe = median(covid_oe, na.rm = TRUE),
    mean_covid_oe   = mean(covid_oe, na.rm = TRUE),
    n_diseases      = n(),
    .groups = "drop"
  )

state_merged <- stringency_by_state %>%
  inner_join(state_median_oe, by = "state")

cat("\n--- Merged state data ---\n")
state_merged %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  as.data.frame() %>%
  print()

# Primary: Spearman rank correlation (n=8 states)
spearman_test <- cor.test(state_merged$mean_stringency,
                          state_merged$median_covid_oe,
                          method = "spearman", exact = FALSE)

cat("\nSpearman correlation: rho =", round(spearman_test$estimate, 4),
    ", p =", round(spearman_test$p.value, 4), "\n")

# Also test with days_above_70
spearman_days70 <- cor.test(state_merged$days_above_70,
                            state_merged$median_covid_oe,
                            method = "spearman", exact = FALSE)

cat("Spearman (days >70): rho =", round(spearman_days70$estimate, 4),
    ", p =", round(spearman_days70$p.value, 4), "\n")

# --- Secondary: Mixed-effects model on disease-state level ---
cat("\nFitting mixed-effects stringency model (disease-state level)...\n")

state_disease_data <- state_disr %>%
  filter(!is.na(covid_oe), is.finite(covid_oe)) %>%
  inner_join(stringency_by_state %>% select(state, mean_stringency, days_above_70),
             by = "state") %>%
  mutate(
    mean_stringency_z = scale(mean_stringency)[, 1]
  )

cat("Disease-state observations:", nrow(state_disease_data), "\n")

mixed_stringency <- tryCatch({
  lmer(covid_oe ~ mean_stringency_z + (1 | disease),
       data = state_disease_data, REML = FALSE)
}, error = function(e) {
  cat("  Mixed model failed:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(mixed_stringency)) {
  cat("Mixed model results:\n")
  tidy(mixed_stringency, effects = "fixed", conf.int = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(., 4))) %>%
    as.data.frame() %>%
    print()
  cat("  AIC:", round(AIC(mixed_stringency), 1), "\n")
}

# Save stringency results
stringency_results <- tibble(
  analysis = c("spearman_mean_stringency", "spearman_days_above_70"),
  rho      = c(spearman_test$estimate, spearman_days70$estimate),
  p_value  = c(spearman_test$p.value, spearman_days70$p.value),
  n        = c(nrow(state_merged), nrow(state_merged))
)

if (!is.null(mixed_stringency)) {
  mixed_coefs <- tidy(mixed_stringency, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "mean_stringency_z") %>%
    mutate(
      p_value  = 2 * pnorm(abs(statistic), lower.tail = FALSE),
      analysis = "mixed_effects_stringency_z",
      rho = estimate, n = nrow(state_disease_data)
    ) %>%
    select(analysis, rho, p_value, n)
  stringency_results <- bind_rows(stringency_results, mixed_coefs)
}

write_csv(stringency_results, file.path(outdir, "stringency_correlation.csv"))
write_csv(state_merged, file.path(tabdir, "stringency_table.csv"))
cat("\nSaved stringency results\n")

# --- Figure 7: Stringency scatter ---
cat("\nGenerating Figure 7...\n")

fig7 <- ggplot(state_merged, aes(x = mean_stringency, y = median_covid_oe)) +
  geom_point(aes(size = n_diseases), colour = "#2166AC", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#B2182B",
              fill = "#FDDBC7", linewidth = 0.8) +
  geom_text(aes(label = state), nudge_y = 0.03, size = 3.5, fontface = "bold") +
  scale_size_continuous(name = "N diseases", range = c(3, 8)) +
  annotate("text", x = max(state_merged$mean_stringency) - 2,
           y = max(state_merged$median_covid_oe) - 0.02,
           label = paste0("Spearman \u03C1 = ",
                          round(spearman_test$estimate, 3),
                          "\np = ", round(spearman_test$p.value, 3)),
           hjust = 1, vjust = 1, size = 4, fontface = "italic") +
  labs(
    title    = "Figure 7. State-Level Stringency and Disease Disruption",
    subtitle = "Mean Oxford Stringency Index (Mar 2020 \u2013 Dec 2021) vs Median O/E Ratio",
    x        = "Mean Stringency Index",
    y        = "Median COVID-period O/E Ratio",
    caption  = "Each point represents one Australian state/territory.\nHigher stringency is associated with lower O/E ratios (greater disruption)."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    plot.caption  = element_text(size = 8, hjust = 0),
    legend.position = "right"
  )

ggsave(file.path(figdir, "fig7_stringency_scatter.png"), fig7,
       width = 9, height = 7, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig7_stringency_scatter.pdf"), fig7,
       width = 9, height = 7, bg = "white")

cat("Saved Figure 7\n\n")

# =============================================================================
# Part C: Recovery Projections
# =============================================================================

cat("--- Part C: Recovery Projections ---\n\n")

# Identify still-suppressed diseases
still_suppressed <- recovery %>%
  filter(still_suppressed == TRUE | trajectory_refined %in%
           c("sustained_suppression", "partial_slow", "partial_nearing"))

cat("Still-suppressed diseases:", nrow(still_suppressed), "\n")
cat(paste(" ", still_suppressed$disease, collapse = "\n"), "\n\n")

# Get post-COVID O/E ratios by year
post_oe <- abl %>%
  filter(level == "national",
         disease %in% still_suppressed$disease,
         year >= 2022, year <= 2025) %>%
  select(disease, year, oe_ratio) %>%
  filter(!is.na(oe_ratio), is.finite(oe_ratio))

# Linear extrapolation: project when O/E reaches 0.9
projections <- list()

for (dis in unique(post_oe$disease)) {
  d <- post_oe %>% filter(disease == dis)

  if (nrow(d) < 2) {
    projections[[dis]] <- tibble(
      disease = dis,
      n_points = nrow(d),
      slope = NA_real_,
      current_oe = d$oe_ratio[d$year == max(d$year)],
      projected_recovery_year = NA_real_,
      method = "insufficient_data"
    )
    next
  }

  fit <- lm(oe_ratio ~ year, data = d)
  slope <- coef(fit)["year"]
  intercept <- coef(fit)["(Intercept)"]
  current_oe <- tail(d$oe_ratio, 1)

  # Project year when O/E = 0.9
  if (slope > 0.001) {  # meaningful positive slope
    recovery_year <- (0.9 - intercept) / slope
    if (recovery_year < max(d$year)) {
      recovery_year <- NA  # already past
      method <- "already_recovered"
    } else if (recovery_year > 2040) {
      method <- "slow_recovery"
    } else {
      method <- "projected"
    }
  } else if (current_oe >= 0.9) {
    recovery_year <- NA
    method <- "already_recovered"
  } else {
    recovery_year <- NA
    method <- "not_recovering"
  }

  projections[[dis]] <- tibble(
    disease = dis,
    n_points = nrow(d),
    slope = round(slope, 4),
    current_oe = round(current_oe, 3),
    projected_recovery_year = if (!is.na(recovery_year)) ceiling(recovery_year) else NA_integer_,
    method = method
  )
}

recovery_projections <- bind_rows(projections)

cat("--- Recovery Projections ---\n")
recovery_projections %>%
  as.data.frame() %>%
  print()

write_csv(recovery_projections, file.path(outdir, "recovery_projections.csv"))
cat("\nSaved recovery projections:", nrow(recovery_projections), "diseases\n")

cat("\n=== Script 3 complete ===\n")
