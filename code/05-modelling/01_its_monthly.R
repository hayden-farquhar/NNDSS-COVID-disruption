#!/usr/bin/env Rscript
# =============================================================================
# 01_its_monthly.R
#
# Interrupted Time Series (ITS) models for 3 Tier 1 diseases with monthly data:
#   - Influenza (laboratory confirmed)
#   - Meningococcal disease (invasive)
#   - Salmonellosis
#
# Model:
#   log(count+1) ~ time + covid + time_since_covid + post + time_since_post
#                + sin(2pi*month/12) + cos(2pi*month/12)
#                + sin(4pi*month/12) + cos(4pi*month/12)
#
# Fitted with nlme::gls(), correlation = corAR1(form = ~time), method = "ML"
# Salmonellosis: try corARMA(p=2) if AR(1) residuals still autocorrelated
# Fallback: OLS if GLS fails to converge
#
# Inputs:
#   data/processed/tier1_monthly_national.csv
#
# Outputs:
#   data/processed/its_monthly_results.csv    — coefficients, CIs, p-values, diagnostics
#   data/processed/its_monthly_fitted.csv     — fitted values + counterfactual
#   output/figures/fig6_its_monthly.png/.pdf  — observed + ITS fitted + counterfactual
#
# Reference: Bernal et al. 2017, Int J Epidemiol (ITS tutorial)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(nlme)
  library(broom)
  library(patchwork)
})

outdir <- "data/processed"
figdir <- "output/figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Phase 5, Script 1: Monthly ITS Models ===\n\n")

# =============================================================================
# 1. Load and prepare data
# =============================================================================

monthly <- read_csv(file.path(outdir, "tier1_monthly_national.csv"),
                    show_col_types = FALSE)

cat("Loaded monthly data:", nrow(monthly), "rows,",
    n_distinct(monthly$disease), "diseases\n")

# Use baseline period 2015-2019, full data through Dec 2024
# Truncate: drop months with zero counts that represent missing data (2025 partial)
monthly <- monthly %>%
  filter(year >= 2015, !(year >= 2025 & month > 1)) %>%
  # For meningococcal, data ends Dec 2024; for flu/salm, Jan 2025 exists but truncate
  filter(year <= 2024)

cat("After filtering to 2015-2024:", nrow(monthly), "rows\n\n")

# Create ITS variables
monthly <- monthly %>%
  group_by(disease) %>%
  arrange(disease, year, month) %>%
  mutate(
    # Sequential time index
    time = row_number(),
    # Date for plotting
    date = as.Date(paste(year, month, "01", sep = "-")),
    # COVID indicator: 1 from March 2020
    covid = as.integer(date >= as.Date("2020-03-01")),
    # Time since COVID onset
    time_since_covid = cumsum(covid),
    # Post-COVID indicator: 1 from January 2023
    post = as.integer(date >= as.Date("2023-01-01")),
    # Time since post-COVID onset
    time_since_post = cumsum(post),
    # Harmonic seasonal terms (first and second order)
    sin1 = sin(2 * pi * month / 12),
    cos1 = cos(2 * pi * month / 12),
    sin2 = sin(4 * pi * month / 12),
    cos2 = cos(4 * pi * month / 12),
    # Log-transformed outcome
    log_count = log(count + 1)
  ) %>%
  ungroup()

diseases <- unique(monthly$disease)

# =============================================================================
# 2. Fit ITS models
# =============================================================================

its_formula <- log_count ~ time + covid + time_since_covid + post + time_since_post +
  sin1 + cos1 + sin2 + cos2

results_list <- list()
fitted_list  <- list()

for (dis in diseases) {
  cat("Fitting ITS for:", dis, "\n")
  d <- monthly %>% filter(disease == dis)

  # --- Try GLS with AR(1) ---
  model      <- NULL
  model_type <- "gls_ar1"
  ar_corr    <- NULL

  tryCatch({
    model <- gls(its_formula, data = d,
                 correlation = corAR1(form = ~ time),
                 method = "ML")
    ar_corr <- coef(model$modelStruct$corStruct, unconstrained = FALSE)
    cat("  GLS AR(1) converged. Phi =", round(ar_corr, 3), "\n")
  }, error = function(e) {
    cat("  GLS AR(1) failed:", conditionMessage(e), "\n")
  })

  # --- Salmonellosis: always try corARMA(p=2) (known residual autocorrelation) ---
  if (!is.null(model) && grepl("Salmonella", dis)) {
    cat("  Trying ARMA(2,0) for Salmonella (known autocorrelation issue)...\n")
    tryCatch({
      model_arma <- gls(its_formula, data = d,
                        correlation = corARMA(form = ~ time, p = 2, q = 0),
                        method = "ML")
      if (AIC(model_arma) < AIC(model)) {
        model      <- model_arma
        model_type <- "gls_arma2"
        ar_corr    <- coef(model$modelStruct$corStruct, unconstrained = FALSE)
        cat("  ARMA(2,0) preferred (lower AIC). Phi =",
            round(ar_corr, 3), "\n")
      } else {
        cat("  AR(1) preferred (lower AIC), keeping original\n")
      }
    }, error = function(e) {
      cat("  ARMA(2,0) failed:", conditionMessage(e), ", keeping AR(1)\n")
    })
  }

  # --- High AR(1): if phi > 0.85, note near-unit-root concern ---
  if (!is.null(model) && inherits(model, "gls")) {
    phi_check <- coef(model$modelStruct$corStruct, unconstrained = FALSE)
    if (length(phi_check) >= 1 && phi_check[1] > 0.85) {
      cat("  NOTE: AR(1) phi =", round(phi_check[1], 3),
          "> 0.85 — near unit root. COVID level-shift test has reduced power.\n")
      cat("  Also fitting OLS (no autocorrelation) for comparison...\n")
      model_ols  <- lm(its_formula, data = d)
      ols_covid  <- coef(model_ols)["covid"]
      ols_p      <- summary(model_ols)$coefficients["covid", "Pr(>|t|)"]
      cat("  OLS COVID coef =", round(ols_covid, 3),
          ", p =", round(ols_p, 4), "(for reference)\n")
    }
  }

  # --- Fallback: OLS ---
  if (is.null(model)) {
    cat("  Falling back to OLS\n")
    model      <- lm(its_formula, data = d)
    model_type <- "ols"
    ar_corr    <- NA
  }

  # --- Extract coefficients ---
  if (inherits(model, "gls")) {
    tidy_coefs <- broom.mixed::tidy(model, conf.int = TRUE)
  } else {
    tidy_coefs <- broom::tidy(model, conf.int = TRUE)
  }

  # Convert to percent change: (exp(beta) - 1) * 100
  tidy_coefs <- tidy_coefs %>%
    mutate(
      pct_change   = (exp(estimate) - 1) * 100,
      pct_change_lo = (exp(conf.low) - 1) * 100,
      pct_change_hi = (exp(conf.high) - 1) * 100,
      disease      = dis,
      model_type   = model_type,
      ar_phi       = if (length(ar_corr) >= 1) ar_corr[1] else NA_real_
    )

  # --- Diagnostics ---
  if (inherits(model, "gls")) {
    resid_norm <- residuals(model, type = "normalized")
  } else {
    resid_norm <- residuals(model)
  }

  lb12     <- Box.test(resid_norm, lag = 12, type = "Ljung-Box")
  lb24     <- Box.test(resid_norm, lag = 24, type = "Ljung-Box")
  sw_test  <- shapiro.test(resid_norm)
  dw_stat  <- sum(diff(resid_norm)^2) / sum(resid_norm^2)  # manual DW

  tidy_coefs <- tidy_coefs %>%
    mutate(
      ljung_box_12_p = lb12$p.value,
      ljung_box_24_p = lb24$p.value,
      shapiro_p      = sw_test$p.value,
      durbin_watson  = dw_stat,
      n_obs          = nrow(d),
      aic            = if (inherits(model, "gls")) AIC(model) else AIC(model),
      bic            = if (inherits(model, "gls")) BIC(model) else BIC(model)
    )

  results_list[[dis]] <- tidy_coefs

  cat("  Diagnostics: Ljung-Box(12) p =", round(lb12$p.value, 4),
      ", Shapiro p =", round(sw_test$p.value, 4),
      ", DW =", round(dw_stat, 3), "\n")

  # --- Fitted values and counterfactual ---
  d_fitted <- d %>%
    mutate(
      fitted = as.numeric(fitted(model)),
      fitted_count = exp(fitted) - 1
    )

  # Counterfactual: what would have happened without COVID/post
  d_counter <- d %>%
    mutate(covid = 0L, time_since_covid = 0L, post = 0L, time_since_post = 0L)

  if (inherits(model, "gls")) {
    d_fitted$counterfactual <- as.numeric(predict(model, newdata = d_counter))
  } else {
    d_fitted$counterfactual <- predict(model, newdata = d_counter)
  }

  d_fitted <- d_fitted %>%
    mutate(counterfactual_count = exp(counterfactual) - 1)

  fitted_list[[dis]] <- d_fitted

  cat("\n")
}

# =============================================================================
# 3. Combine and save results
# =============================================================================

its_results <- bind_rows(results_list) %>%
  select(disease, term, estimate, std.error, statistic, p.value,
         conf.low, conf.high, pct_change, pct_change_lo, pct_change_hi,
         model_type, ar_phi, ljung_box_12_p, ljung_box_24_p,
         shapiro_p, durbin_watson, n_obs, aic, bic)

its_fitted <- bind_rows(fitted_list) %>%
  select(disease, year, month, date, time, count, log_count,
         covid, post, fitted, fitted_count,
         counterfactual, counterfactual_count)

write_csv(its_results, file.path(outdir, "its_monthly_results.csv"))
write_csv(its_fitted, file.path(outdir, "its_monthly_fitted.csv"))

cat("Saved:", nrow(its_results), "coefficient rows to its_monthly_results.csv\n")
cat("Saved:", nrow(its_fitted), "fitted value rows to its_monthly_fitted.csv\n\n")

# Print key results summary
cat("=== Key Results Summary ===\n\n")
its_results %>%
  filter(term == "covid") %>%
  select(disease, estimate, pct_change, p.value, model_type, ar_phi,
         ljung_box_12_p) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  as.data.frame() %>%
  print()

# =============================================================================
# 4. Figure 6: Monthly ITS panels
# =============================================================================

cat("\n\nGenerating Figure 6...\n")

# Colour scheme
col_obs   <- "grey30"
col_fit   <- "#2166AC"
col_cf    <- "#B2182B"
col_shade <- "#FFEDA0"

disease_labels <- c(

  "Influenza (laboratory confirmed)" = "Influenza",
  "Meningococcal disease (invasive)" = "Meningococcal disease",
  "Salmonellosis"                    = "Salmonellosis"
)

plot_list <- list()

for (dis in diseases) {
  d <- its_fitted %>% filter(disease == dis)

  # COVID and post-COVID shading regions
  covid_start <- as.Date("2020-03-01")
  post_start  <- as.Date("2023-01-01")
  end_date    <- max(d$date)

  p <- ggplot(d, aes(x = date)) +
    # Shading for periods
    annotate("rect", xmin = covid_start, xmax = post_start,
             ymin = -Inf, ymax = Inf, fill = "#FEE0D2", alpha = 0.4) +
    annotate("rect", xmin = post_start, xmax = end_date,
             ymin = -Inf, ymax = Inf, fill = "#DEEBF7", alpha = 0.4) +
    # Observed counts
    geom_line(aes(y = count, colour = "Observed"), linewidth = 0.3, alpha = 0.7) +
    geom_point(aes(y = count, colour = "Observed"), size = 0.4, alpha = 0.5) +
    # ITS fitted
    geom_line(aes(y = fitted_count, colour = "ITS fitted"), linewidth = 0.7) +
    # Counterfactual
    geom_line(aes(y = counterfactual_count, colour = "Counterfactual"),
              linewidth = 0.7, linetype = "dashed") +
    # Intervention lines
    geom_vline(xintercept = covid_start, linetype = "dotted", colour = "red", linewidth = 0.5) +
    geom_vline(xintercept = post_start, linetype = "dotted", colour = "blue", linewidth = 0.5) +
    scale_colour_manual(
      values = c("Observed" = col_obs, "ITS fitted" = col_fit,
                 "Counterfactual" = col_cf),
      name = NULL
    ) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title = disease_labels[dis],
      x = NULL,
      y = "Monthly notifications"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11)
    )

  # Log scale for influenza (very large dynamic range)
  if (grepl("Influenza", dis)) {
    p <- p + scale_y_continuous(
      labels = scales::comma,
      trans  = "pseudo_log"
    )
  } else {
    p <- p + scale_y_continuous(labels = scales::comma)
  }

  plot_list[[dis]] <- p
}

# Combine panels
fig6 <- wrap_plots(plot_list, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

fig6 <- fig6 + plot_annotation(
  title   = "Interrupted Time Series Analysis of Monthly Notifications",
  caption = "Red shading: COVID-acute period (Mar 2020 \u2013 Dec 2022). Blue shading: Post-COVID (Jan 2023+).\nDotted lines mark intervention points. Counterfactual shows expected trajectory without COVID.",
  theme   = theme(
    plot.title   = element_text(face = "bold", size = 13),
    plot.caption = element_text(size = 8, hjust = 0)
  )
)

ggsave(file.path(figdir, "fig6_its_monthly.png"), fig6,
       width = 10, height = 12, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig6_its_monthly.pdf"), fig6,
       width = 10, height = 12, bg = "white")

cat("Saved Figure 6 to", figdir, "\n")
cat("\n=== Script 1 complete ===\n")
