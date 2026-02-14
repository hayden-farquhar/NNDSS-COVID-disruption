#!/usr/bin/env Rscript
# =============================================================================
# 05_survival_gam.R
#
# Phase 5b Script 2: Survival Analysis + GAM Recovery + Bootstrap Immunity Debt
#
# Part A: Time-to-recovery survival analysis
#   - Kaplan-Meier curves stratified by border_sensitive, VPD, cluster
#   - Cox proportional hazards model for recovery predictors
#
# Part B: GAM recovery curves for 3 monthly Tier 1 diseases
#   - gamm() with AR(1) for fair comparison to GLS ITS
#   - Visual comparison: GAM smooth vs ITS piecewise linear
#
# Part C: Bootstrap immunity debt uncertainty
#   - 999 bootstrap resamples of baseline years
#   - 95% CIs on cumulative deficit, excess, and net debt
#
# Inputs:
#   data/processed/recovery_trajectories.csv
#   data/processed/trajectory_clusters.csv
#   data/processed/disruption_metrics_national.csv
#   data/processed/tier1_monthly_national.csv
#   data/processed/its_monthly_fitted.csv
#   data/processed/annual_with_baselines.csv
#
# Outputs:
#   data/processed/survival_recovery.csv
#   data/processed/gam_monthly_results.csv
#   data/processed/gam_monthly_fitted.csv
#   data/processed/immunity_debt_bootstrap.csv
#   output/figures/fig10_km_survival.png/.pdf
#   output/figures/fig11_gam_recovery.png/.pdf
#   output/figures/fig12_debt_forest.png/.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(mgcv)
  library(nlme)
  library(patchwork)
})

outdir <- "data/processed"
figdir <- "output/figures"

cat("=== Phase 5b, Script 5: Survival + GAM + Bootstrap ===\n\n")

# =============================================================================
# 0. Load data
# =============================================================================

recovery   <- read_csv(file.path(outdir, "recovery_trajectories.csv"), show_col_types = FALSE)
clusters   <- read_csv(file.path(outdir, "trajectory_clusters.csv"), show_col_types = FALSE)
disruption <- read_csv(file.path(outdir, "disruption_metrics_national.csv"), show_col_types = FALSE)
monthly    <- read_csv(file.path(outdir, "tier1_monthly_national.csv"), show_col_types = FALSE)
its_fitted <- read_csv(file.path(outdir, "its_monthly_fitted.csv"), show_col_types = FALSE)
abl        <- read_csv(file.path(outdir, "annual_with_baselines.csv"), show_col_types = FALSE)

# =============================================================================
# PART A: Time-to-Recovery Survival Analysis
# =============================================================================

cat("--- Part A: Time-to-Recovery Survival Analysis ---\n\n")

# Build survival dataset: diseases suppressed during COVID (O/E < 0.9)
surv_data <- recovery %>%
  filter(covid_oe_rate < 0.9) %>%
  left_join(disruption %>% select(disease, baseline_mean_rate), by = "disease") %>%
  left_join(clusters %>% select(disease, cluster), by = "disease") %>%
  mutate(
    event = as.integer(!is.na(recovery_year)),
    time  = if_else(event == 1L, as.integer(recovery_year - 2019), 6L),
    time  = pmax(time, 1L),
    log_baseline   = log1p(baseline_mean_rate),
    disruption_mag = covid_oe_rate,
    cluster_f      = factor(cluster)
  )

n_total  <- nrow(surv_data)
n_events <- sum(surv_data$event)
cat("Survival dataset:", n_total, "diseases suppressed (O/E < 0.9)\n")
cat("Events (recovered):", n_events, "; Censored:", n_total - n_events, "\n\n")

# --- Kaplan-Meier: Overall ---
km_overall <- survfit(Surv(time, event) ~ 1, data = surv_data)
cat("Median time-to-recovery:", summary(km_overall)$table["median"], "years\n")

# --- KM by border_sensitive ---
km_border <- survfit(Surv(time, event) ~ border_sensitive, data = surv_data)
lr_border <- survdiff(Surv(time, event) ~ border_sensitive, data = surv_data)
cat("Log-rank test (border_sensitive): chi2 =",
    round(lr_border$chisq, 2), ", p =",
    round(pchisq(lr_border$chisq, 1, lower.tail = FALSE), 4), "\n")

# --- KM by vaccine_preventable ---
km_vpd <- survfit(Surv(time, event) ~ vaccine_preventable, data = surv_data)
lr_vpd <- survdiff(Surv(time, event) ~ vaccine_preventable, data = surv_data)
cat("Log-rank test (VPD): chi2 =",
    round(lr_vpd$chisq, 2), ", p =",
    round(pchisq(lr_vpd$chisq, 1, lower.tail = FALSE), 4), "\n")

# --- KM by cluster ---
if (n_distinct(surv_data$cluster_f) > 1) {
  km_cluster <- survfit(Surv(time, event) ~ cluster_f, data = surv_data)
  lr_cluster <- survdiff(Surv(time, event) ~ cluster_f, data = surv_data)
  cat("Log-rank test (cluster): chi2 =",
      round(lr_cluster$chisq, 2), ", p =",
      round(pchisq(lr_cluster$chisq, df = n_distinct(surv_data$cluster_f) - 1,
                    lower.tail = FALSE), 4), "\n")
}

# --- Cox PH: Univariate screening ---
cat("\nCox PH univariate screening:\n")
covariates <- c("border_sensitive", "vaccine_preventable", "disruption_mag", "log_baseline")
univ_results <- list()
for (cv in covariates) {
  f <- as.formula(paste0("Surv(time, event) ~ ", cv))
  cx <- coxph(f, data = surv_data)
  s  <- summary(cx)
  univ_results[[cv]] <- tibble(
    covariate = cv,
    HR = s$conf.int[1, "exp(coef)"],
    HR_lo = s$conf.int[1, "lower .95"],
    HR_hi = s$conf.int[1, "upper .95"],
    p = s$coefficients[1, "Pr(>|z|)"]
  )
  cat("  ", cv, ": HR =", round(s$conf.int[1, 1], 3),
      "(", round(s$conf.int[1, 3], 2), "-", round(s$conf.int[1, 4], 2), ")",
      "p =", round(s$coefficients[1, 5], 4), "\n")
}

# --- Cox PH: Multivariate (parsimonious — 2-3 covariates for n events) ---
cat("\nCox PH multivariate model:\n")
# Include covariates based on sample size (rule: ~10 events per covariate)
max_covars <- floor(n_events / 10)
cat("  Max covariates (rule of 10):", max_covars, "\n")

# Always include disruption_mag (our key hypothesis covariate)
# Add border_sensitive if room
if (max_covars >= 3) {
  cox_mv <- coxph(Surv(time, event) ~ disruption_mag + border_sensitive + vaccine_preventable,
                  data = surv_data)
} else if (max_covars >= 2) {
  cox_mv <- coxph(Surv(time, event) ~ disruption_mag + border_sensitive,
                  data = surv_data)
} else {
  cox_mv <- coxph(Surv(time, event) ~ disruption_mag, data = surv_data)
}

cox_summary <- summary(cox_mv)
cat("\n")
print(cox_summary)

# PH assumption test
ph_test <- cox.zph(cox_mv)
cat("\nProportional hazards assumption test:\n")
print(ph_test)
cat("\n")

# --- Save survival results ---
surv_results <- surv_data %>%
  select(disease, transmission_mode, border_sensitive, vaccine_preventable,
         covid_oe_rate, post_oe_rate, cluster, time, event,
         disruption_mag, log_baseline)

cox_coefs <- broom::tidy(cox_mv, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(model = "cox_multivariate")

univ_df <- bind_rows(univ_results)

write_csv(surv_results, file.path(outdir, "survival_recovery.csv"))
write_csv(bind_rows(cox_coefs, univ_df %>% mutate(model = "cox_univariate")),
          file.path(outdir, "survival_cox_results.csv"))
cat("Saved survival results\n\n")

# --- Figure 10: Kaplan-Meier curves ---
cat("Generating Figure 10...\n")

# Panel A: by border_sensitive
p10a <- ggsurvplot(
  km_border, data = surv_data,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, risk.table.height = 0.25,
  palette = c("#4575B4", "#D73027"),
  legend.labs = c("Not border-sensitive", "Border-sensitive"),
  legend.title = "",
  xlab = "Years since COVID onset (2020)",
  ylab = "Proportion still suppressed",
  title = "A. Recovery by border sensitivity",
  ggtheme = theme_minimal(base_size = 10)
)

# Panel B: by VPD
p10b <- ggsurvplot(
  km_vpd, data = surv_data,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, risk.table.height = 0.25,
  palette = c("#4575B4", "#D73027"),
  legend.labs = c("Not VPD", "VPD"),
  legend.title = "",
  xlab = "Years since COVID onset (2020)",
  ylab = "Proportion still suppressed",
  title = "B. Recovery by vaccine-preventable status",
  ggtheme = theme_minimal(base_size = 10)
)

# Combine using arrange_ggsurvplots
fig10_list <- list(p10a, p10b)
fig10 <- arrange_ggsurvplots(fig10_list, ncol = 2, nrow = 1,
                              title = "Figure 10. Time-to-Recovery Survival Analysis")

ggsave(file.path(figdir, "fig10_km_survival.png"), fig10,
       width = 14, height = 7, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig10_km_survival.pdf"), fig10,
       width = 14, height = 7, bg = "white")
cat("Saved Figure 10\n\n")

# =============================================================================
# PART B: GAM Recovery Curves for Monthly Tier 1 Diseases
# =============================================================================

cat("--- Part B: GAM Recovery Curves ---\n\n")

# Prepare monthly data (same filter as Script 1: 2015-2024)
monthly_prep <- monthly %>%
  filter(year >= 2015, year <= 2024) %>%
  group_by(disease) %>%
  arrange(disease, year, month) %>%
  mutate(
    time      = row_number(),
    date      = as.Date(paste(year, month, "01", sep = "-")),
    log_count = log(count + 1),
    covid     = as.integer(date >= as.Date("2020-03-01")),
    post      = as.integer(date >= as.Date("2023-01-01"))
  ) %>%
  ungroup()

diseases <- unique(monthly_prep$disease)
gam_results_list <- list()
gam_fitted_list  <- list()

for (dis in diseases) {
  cat("Fitting GAM for:", dis, "\n")
  d <- monthly_prep %>% filter(disease == dis)

  # GAM with AR(1) for fair comparison to GLS ITS
  gam_fit <- NULL
  gam_type <- "gamm_ar1"

  tryCatch({
    gam_fit <- gamm(
      log_count ~ s(time, k = 20) + s(month, bs = "cc", k = 12),
      correlation = corAR1(form = ~ time),
      data = d, method = "ML"
    )
    cat("  gamm AR(1) converged\n")
  }, error = function(e) {
    cat("  gamm AR(1) failed:", conditionMessage(e), "\n")
  })

  # Fallback: simple gam without AR(1)
  if (is.null(gam_fit)) {
    cat("  Falling back to gam (no AR1)\n")
    gam_fit_simple <- gam(
      log_count ~ s(time, k = 20) + s(month, bs = "cc", k = 12),
      data = d, method = "ML"
    )
    # Wrap in list to match gamm structure
    gam_fit <- list(gam = gam_fit_simple, lme = NULL)
    gam_type <- "gam_simple"
  }

  # Extract diagnostics
  gam_aic <- if (!is.null(gam_fit$lme)) AIC(gam_fit$lme) else AIC(gam_fit$gam)
  gam_edf <- sum(gam_fit$gam$edf)
  gam_r2  <- summary(gam_fit$gam)$r.sq

  # Fitted values
  gam_fitted_vals <- as.numeric(predict(gam_fit$gam))

  # Compare to ITS
  its_d <- its_fitted %>% filter(disease == dis)
  if (nrow(its_d) == nrow(d)) {
    rmse_its <- sqrt(mean((d$log_count - its_d$fitted)^2))
    rmse_gam <- sqrt(mean((d$log_count - gam_fitted_vals)^2))
    cat("  RMSE: ITS =", round(rmse_its, 4), ", GAM =", round(rmse_gam, 4),
        "(", if_else(rmse_gam < rmse_its, "GAM better", "ITS better"), ")\n")
  } else {
    rmse_its <- NA
    rmse_gam <- sqrt(mean((d$log_count - gam_fitted_vals)^2))
    cat("  RMSE: GAM =", round(rmse_gam, 4), "(ITS row count mismatch)\n")
  }

  gam_results_list[[dis]] <- tibble(
    disease  = dis,
    gam_type = gam_type,
    gam_aic  = gam_aic,
    gam_edf  = gam_edf,
    gam_r2   = gam_r2,
    rmse_gam = rmse_gam,
    rmse_its = rmse_its
  )

  gam_fitted_list[[dis]] <- d %>%
    select(disease, year, month, date, time, count, log_count) %>%
    mutate(gam_fitted_log = gam_fitted_vals,
           gam_fitted_count = exp(gam_fitted_vals) - 1)

  cat("\n")
}

gam_results <- bind_rows(gam_results_list)
gam_fitted_all <- bind_rows(gam_fitted_list)

write_csv(gam_results, file.path(outdir, "gam_monthly_results.csv"))
write_csv(gam_fitted_all, file.path(outdir, "gam_monthly_fitted.csv"))

cat("GAM vs ITS comparison:\n")
print(as.data.frame(gam_results %>% mutate(across(where(is.numeric), ~ round(., 4)))))
cat("\n")

# --- Figure 11: GAM vs ITS overlay ---
cat("Generating Figure 11...\n")

disease_labels <- c(
  "Influenza (laboratory confirmed)" = "Influenza",
  "Meningococcal disease (invasive)" = "Meningococcal disease",
  "Salmonellosis"                    = "Salmonellosis"
)

plot_list_11 <- list()
for (dis in diseases) {
  gam_d <- gam_fitted_all %>% filter(disease == dis)
  its_d <- its_fitted %>% filter(disease == dis)

  # Merge ITS fitted values if row counts match
  if (nrow(its_d) == nrow(gam_d)) {
    plot_d <- gam_d %>%
      mutate(its_fitted_count = its_d$fitted_count,
             its_cf_count = its_d$counterfactual_count)
  } else {
    plot_d <- gam_d %>%
      mutate(its_fitted_count = NA_real_, its_cf_count = NA_real_)
  }

  covid_start <- as.Date("2020-03-01")
  post_start  <- as.Date("2023-01-01")

  p <- ggplot(plot_d, aes(x = date)) +
    annotate("rect", xmin = covid_start, xmax = post_start,
             ymin = -Inf, ymax = Inf, fill = "#FEE0D2", alpha = 0.4) +
    annotate("rect", xmin = post_start, xmax = max(plot_d$date),
             ymin = -Inf, ymax = Inf, fill = "#DEEBF7", alpha = 0.4) +
    geom_point(aes(y = count), colour = "grey50", size = 0.4, alpha = 0.5) +
    geom_line(aes(y = its_fitted_count, colour = "ITS (piecewise linear)"),
              linewidth = 0.6, na.rm = TRUE) +
    geom_line(aes(y = its_cf_count, colour = "ITS counterfactual"),
              linewidth = 0.6, linetype = "dashed", na.rm = TRUE) +
    geom_line(aes(y = gam_fitted_count, colour = "GAM (smooth)"),
              linewidth = 0.8) +
    geom_vline(xintercept = covid_start, linetype = "dotted", colour = "red",
               linewidth = 0.4) +
    geom_vline(xintercept = post_start, linetype = "dotted", colour = "blue",
               linewidth = 0.4) +
    scale_colour_manual(
      values = c("ITS (piecewise linear)" = "#2166AC",
                 "ITS counterfactual" = "#B2182B",
                 "GAM (smooth)" = "#1B7837"),
      name = NULL
    ) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = disease_labels[dis], x = NULL, y = "Monthly notifications") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))

  if (grepl("Influenza", dis)) {
    p <- p + scale_y_continuous(labels = scales::comma, trans = "pseudo_log")
  } else {
    p <- p + scale_y_continuous(labels = scales::comma)
  }

  plot_list_11[[dis]] <- p
}

fig11 <- wrap_plots(plot_list_11, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

fig11 <- fig11 + plot_annotation(
  title = "Figure 11. GAM vs ITS Recovery Curves",
  caption = paste0("Green: GAM smooth (captures non-linear recovery). ",
                   "Blue: ITS piecewise linear. Red dashed: ITS counterfactual."),
  theme = theme(
    plot.title   = element_text(face = "bold", size = 13),
    plot.caption = element_text(size = 8, hjust = 0)
  )
)

ggsave(file.path(figdir, "fig11_gam_recovery.png"), fig11,
       width = 10, height = 12, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig11_gam_recovery.pdf"), fig11,
       width = 10, height = 12, bg = "white")
cat("Saved Figure 11\n\n")

# =============================================================================
# PART C: Bootstrap Immunity Debt Uncertainty
# =============================================================================

cat("--- Part C: Bootstrap Immunity Debt CIs ---\n\n")

national <- abl %>%
  filter(level == "national", disease %in% recovery$disease,
         year >= 2015, year <= 2025)

analysis_diseases <- unique(recovery$disease)
B <- 999
set.seed(42)

cat("Running", B, "bootstrap iterations for", length(analysis_diseases), "diseases...\n")

boot_results <- list()

for (dis in analysis_diseases) {
  d <- national %>% filter(disease == dis) %>% arrange(year)

  baseline_counts <- d %>% filter(period == "baseline") %>% pull(count)
  post_data       <- d %>% filter(year >= 2020) %>% arrange(year)
  post_counts     <- post_data$count

  if (length(baseline_counts) < 3 || length(post_counts) == 0) {
    cat("  Skipping", dis, "(insufficient data)\n")
    next
  }

  # Original point estimate
  orig_expected <- mean(baseline_counts)
  orig_deficit  <- sum(pmax(0, orig_expected - post_counts))
  orig_excess   <- sum(pmax(0, post_counts - orig_expected))
  orig_net      <- orig_excess - orig_deficit

  # Bootstrap
  net_boot     <- numeric(B)
  deficit_boot <- numeric(B)
  excess_boot  <- numeric(B)

  for (b in seq_len(B)) {
    boot_baseline <- sample(baseline_counts, length(baseline_counts), replace = TRUE)
    boot_expected <- mean(boot_baseline)
    deficit_boot[b] <- sum(pmax(0, boot_expected - post_counts))
    excess_boot[b]  <- sum(pmax(0, post_counts - boot_expected))
    net_boot[b]     <- excess_boot[b] - deficit_boot[b]
  }

  boot_results[[dis]] <- tibble(
    disease       = dis,
    net_point     = orig_net,
    net_ci_lo     = quantile(net_boot, 0.025),
    net_ci_hi     = quantile(net_boot, 0.975),
    deficit_point = orig_deficit,
    deficit_ci_lo = quantile(deficit_boot, 0.025),
    deficit_ci_hi = quantile(deficit_boot, 0.975),
    excess_point  = orig_excess,
    excess_ci_lo  = quantile(excess_boot, 0.025),
    excess_ci_hi  = quantile(excess_boot, 0.975),
    ci_excludes_zero = (quantile(net_boot, 0.025) > 0) |
                       (quantile(net_boot, 0.975) < 0)
  )
}

boot_df <- bind_rows(boot_results) %>%
  left_join(recovery %>% select(disease, transmission_mode, immunity_debt_signal),
            by = "disease")

write_csv(boot_df, file.path(outdir, "immunity_debt_bootstrap.csv"))

n_sig <- sum(boot_df$ci_excludes_zero)
cat("Diseases with CI excluding zero:", n_sig, "/", nrow(boot_df), "\n")
cat("Diseases with positive net debt (immunity debt):",
    sum(boot_df$net_point > 0), "\n")
cat("Diseases with negative net (still in deficit):",
    sum(boot_df$net_point < 0), "\n\n")

# Print top 10 by absolute net debt
cat("Top 10 diseases by |net debt|:\n")
boot_df %>%
  arrange(desc(abs(net_point))) %>%
  slice_head(n = 10) %>%
  mutate(across(where(is.numeric), ~ round(., 0))) %>%
  select(disease, net_point, net_ci_lo, net_ci_hi, ci_excludes_zero,
         immunity_debt_signal) %>%
  as.data.frame() %>%
  print()
cat("\n")

# --- Figure 12: Immunity debt forest plot ---
cat("Generating Figure 12...\n")

# Order by net debt — keep full names to avoid duplicates
plot_boot <- boot_df %>%
  mutate(
    disease_short = disease,
    disease_short = factor(disease_short, levels = unique(disease_short[order(net_point)])),
    sig_group = case_when(
      net_ci_lo > 0 ~ "Confirmed excess (immunity debt)",
      net_ci_hi < 0 ~ "Confirmed deficit (still suppressed)",
      TRUE          ~ "Inconclusive"
    )
  )

fig12 <- ggplot(plot_boot, aes(x = net_point, y = disease_short, colour = sig_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = net_ci_lo, xmax = net_ci_hi),
                 height = 0.3, linewidth = 0.4) +
  geom_point(size = 1.5) +
  scale_colour_manual(
    values = c("Confirmed excess (immunity debt)" = "#D73027",
               "Confirmed deficit (still suppressed)" = "#4575B4",
               "Inconclusive" = "grey50"),
    name = NULL
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "Figure 12. Cumulative Net Immunity Debt with 95% Bootstrap CIs",
    subtitle = paste0("Positive = more excess cases than deficit (immunity debt). ",
                      "n = ", nrow(plot_boot), " diseases, ", B, " bootstrap iterations."),
    x = "Net cases (cumulative excess \u2212 cumulative deficit, 2020\u20132025)",
    y = NULL
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9)
  )

ggsave(file.path(figdir, "fig12_debt_forest.png"), fig12,
       width = 10, height = 12, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig12_debt_forest.pdf"), fig12,
       width = 10, height = 12, bg = "white")
cat("Saved Figure 12\n")

cat("\n=== Script 5 complete ===\n")
