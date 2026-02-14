#!/usr/bin/env Rscript
# =============================================================================
# 04_did_clustering.R
#
# Phase 5b Script 1: Difference-in-Differences + Trajectory Clustering
#
# Part A: Difference-in-Differences (DiD)
#   Treatment: NPI-sensitive diseases (respiratory OR border-sensitive)
#   Control:   All other diseases (enteric, STI, domestic vector-borne,
#              zoonotic, blood-borne)
#   Three models: simple DiD, three-period DiD (+ recovery), full trend DiD
#   Validation: pre-trend test (parallel trends), placebo test (fake 2018)
#
# Part B: Trajectory Clustering
#   k-medoids (PAM) on 47 × 11 z-scored rate matrix
#   Optimal k by average silhouette width
#   Cluster characterisation by transmission mode and recovery status
#
# Inputs:
#   data/processed/annual_with_baselines.csv
#   data/processed/recovery_trajectories.csv
#
# Outputs:
#   data/processed/did_results.csv
#   data/processed/trajectory_clusters.csv
#   output/figures/fig8_did_panel.png/.pdf
#   output/figures/fig9_cluster_profiles.png/.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(lme4)
  library(broom.mixed)
  library(cluster)
  library(factoextra)
  library(patchwork)
})

outdir <- "data/processed"
figdir <- "output/figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Phase 5b, Script 4: DiD + Trajectory Clustering ===\n\n")

# =============================================================================
# 0. Load data
# =============================================================================

abl <- read_csv(file.path(outdir, "annual_with_baselines.csv"), show_col_types = FALSE)
recovery <- read_csv(file.path(outdir, "recovery_trajectories.csv"), show_col_types = FALSE)

analysis_diseases <- recovery$disease

national <- abl %>%
  filter(level == "national",
         disease %in% analysis_diseases,
         year >= 2015, year <= 2025)

cat("Loaded:", n_distinct(national$disease), "diseases,",
    nrow(national), "observations\n")
cat("Years:", paste(range(national$year), collapse = " to "), "\n\n")

# =============================================================================
# PART A: Difference-in-Differences
# =============================================================================

cat("--- Part A: Difference-in-Differences ---\n\n")

# Z-score rates within each disease
did_data <- national %>%
  group_by(disease) %>%
  mutate(
    rate_z = (rate_per_100k - mean(rate_per_100k, na.rm = TRUE)) /
              sd(rate_per_100k, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(rate_z = if_else(is.nan(rate_z) | is.infinite(rate_z), 0, rate_z))

# ITS variables
did_data <- did_data %>%
  mutate(
    year_c = year - 2017,
    covid = as.integer(year >= 2020),
    time_since_covid = pmax(0L, as.integer(year - 2019)),
    post = as.integer(year >= 2023),
    time_since_post = pmax(0L, as.integer(year - 2022)),
    # NPI-sensitive: respiratory OR border-sensitive
    npi_sensitive = as.integer(transmission_mode == "respiratory" | border_sensitive)
  )

n_treatment <- n_distinct(did_data$disease[did_data$npi_sensitive == 1])
n_control   <- n_distinct(did_data$disease[did_data$npi_sensitive == 0])
cat("Treatment (NPI-sensitive):", n_treatment, "diseases\n")
cat("Control:", n_control, "diseases\n\n")

# Helper: extract coefficients with p-values from lmer
extract_lmer <- function(model, model_name) {
  broom.mixed::tidy(model, conf.int = TRUE, effects = "fixed") %>%
    mutate(
      p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE),
      model   = model_name
    )
}

# --- Model 1: Simple two-period DiD ---
cat("Fitting Model 1: Simple DiD (pre vs COVID)...\n")
m1 <- lmer(rate_z ~ year_c + npi_sensitive + covid + covid:npi_sensitive +
             (1 | disease), data = did_data, REML = FALSE)
m1_tidy <- extract_lmer(m1, "M1_simple_DiD")

did_coef <- filter(m1_tidy, term == "covid:npi_sensitive")
cat("  DiD coef (covid:npi_sensitive):", round(did_coef$estimate, 3),
    ", p =", round(did_coef$p.value, 4), "\n")

# --- Model 2: Three-period DiD (baseline vs COVID vs recovery) ---
cat("Fitting Model 2: Three-period DiD (+ recovery)...\n")
m2 <- lmer(rate_z ~ year_c + npi_sensitive +
             covid + covid:npi_sensitive +
             post + post:npi_sensitive +
             (1 | disease), data = did_data, REML = FALSE)
m2_tidy <- extract_lmer(m2, "M2_three_period")

recov_coef <- filter(m2_tidy, term == "npi_sensitive:post")
if (nrow(recov_coef) == 0)
  recov_coef <- filter(m2_tidy, grepl("post.*npi|npi.*post", term))
cat("  Recovery DiD coef:", round(recov_coef$estimate[1], 3),
    ", p =", round(recov_coef$p.value[1], 4), "\n")

# --- Model 3: Full trend DiD ---
cat("Fitting Model 3: Full trend DiD...\n")
m3 <- lmer(rate_z ~ year_c + npi_sensitive +
             covid + covid:npi_sensitive +
             time_since_covid + time_since_covid:npi_sensitive +
             post + post:npi_sensitive +
             time_since_post + time_since_post:npi_sensitive +
             (1 | disease), data = did_data, REML = FALSE)
m3_tidy <- extract_lmer(m3, "M3_full_trend")

# --- Pre-trend test (2015–2019 only) ---
cat("\nPre-trend test (2015-2019)...\n")
baseline_data <- did_data %>% filter(year <= 2019)
m_pre <- lmer(rate_z ~ year_c + npi_sensitive + year_c:npi_sensitive +
                (1 | disease), data = baseline_data, REML = FALSE)
pre_tidy <- extract_lmer(m_pre, "pretrend_test")
pretrend_p <- filter(pre_tidy, grepl("year_c.*npi|npi.*year_c", term))$p.value[1]
cat("  year_c:npi_sensitive p =", round(pretrend_p, 4),
    if_else(pretrend_p > 0.05, " (PASS: parallel trends)", " (FAIL)"), "\n")

# --- Placebo test (fake intervention at 2018) ---
cat("Placebo test (fake COVID at 2018)...\n")
placebo_data <- baseline_data %>% mutate(fake_covid = as.integer(year >= 2018))
m_plac <- lmer(rate_z ~ year_c + npi_sensitive +
                 fake_covid + fake_covid:npi_sensitive +
                 (1 | disease), data = placebo_data, REML = FALSE)
plac_tidy <- extract_lmer(m_plac, "placebo_test")
placebo_p <- filter(plac_tidy, grepl("fake_covid.*npi|npi.*fake_covid", term))$p.value[1]
cat("  fake_covid:npi_sensitive p =", round(placebo_p, 4),
    if_else(placebo_p > 0.05, " (PASS: no false positive)", " (FAIL)"), "\n\n")

# --- Sensitivity: alternative definitions ---
cat("Sensitivity: border_sensitive only...\n")
did_data_sens <- did_data %>%
  mutate(npi_sensitive_border = as.integer(border_sensitive))
m_sens_border <- lmer(rate_z ~ year_c + npi_sensitive_border +
                        covid + covid:npi_sensitive_border +
                        (1 | disease),
                      data = did_data_sens, REML = FALSE)
sens_border_tidy <- extract_lmer(m_sens_border, "sensitivity_border_only")

cat("Sensitivity: respiratory only...\n")
did_data_sens <- did_data_sens %>%
  mutate(npi_sensitive_resp = as.integer(transmission_mode == "respiratory"))
m_sens_resp <- lmer(rate_z ~ year_c + npi_sensitive_resp +
                      covid + covid:npi_sensitive_resp +
                      (1 | disease),
                    data = did_data_sens, REML = FALSE)
sens_resp_tidy <- extract_lmer(m_sens_resp, "sensitivity_respiratory_only")

# --- Combine and save DiD results ---
did_results <- bind_rows(m1_tidy, m2_tidy, m3_tidy, pre_tidy, plac_tidy,
                         sens_border_tidy, sens_resp_tidy) %>%
  select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high)

write_csv(did_results, file.path(outdir, "did_results.csv"))
cat("Saved DiD results:", nrow(did_results), "rows\n\n")

cat("=== Key DiD Coefficients ===\n")
did_results %>%
  filter(grepl("npi_sensitive", term)) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  as.data.frame() %>%
  print()
cat("\n")

# --- Figure 8: DiD panel ---
cat("Generating Figure 8...\n")

group_means <- did_data %>%
  mutate(group = if_else(npi_sensitive == 1, "NPI-sensitive", "Control")) %>%
  group_by(group, year) %>%
  summarise(
    mean_z  = mean(rate_z, na.rm = TRUE),
    se_z    = sd(rate_z, na.rm = TRUE) / sqrt(n()),
    mean_oe = mean(oe_ratio, na.rm = TRUE),
    se_oe   = sd(oe_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

covid_start <- 2020
post_start  <- 2023

p8a <- ggplot(group_means, aes(x = year, y = mean_z, colour = group, fill = group)) +
  annotate("rect", xmin = covid_start - 0.5, xmax = post_start - 0.5,
           ymin = -Inf, ymax = Inf, fill = "#FEE0D2", alpha = 0.4) +
  annotate("rect", xmin = post_start - 0.5, xmax = 2025.5,
           ymin = -Inf, ymax = Inf, fill = "#DEEBF7", alpha = 0.4) +
  geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = covid_start - 0.5, linetype = "dotted", colour = "red") +
  geom_vline(xintercept = post_start - 0.5, linetype = "dotted", colour = "blue") +
  scale_colour_manual(values = c("NPI-sensitive" = "#D73027", "Control" = "#4575B4"),
                      name = NULL) +
  scale_fill_manual(values = c("NPI-sensitive" = "#D73027", "Control" = "#4575B4"),
                    name = NULL) +
  scale_x_continuous(breaks = 2015:2025) +
  labs(title = "A. Difference-in-Differences: Z-scored rates",
       subtitle = paste0("Treatment: ", n_treatment,
                         " NPI-sensitive; Control: ", n_control, " diseases"),
       x = NULL, y = "Mean z-scored rate (\u00b1 SE)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p8b <- ggplot(group_means, aes(x = year, y = mean_oe, colour = group, fill = group)) +
  annotate("rect", xmin = covid_start - 0.5, xmax = post_start - 0.5,
           ymin = -Inf, ymax = Inf, fill = "#FEE0D2", alpha = 0.4) +
  annotate("rect", xmin = post_start - 0.5, xmax = 2025.5,
           ymin = -Inf, ymax = Inf, fill = "#DEEBF7", alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = mean_oe - se_oe, ymax = mean_oe + se_oe),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = covid_start - 0.5, linetype = "dotted", colour = "red") +
  geom_vline(xintercept = post_start - 0.5, linetype = "dotted", colour = "blue") +
  scale_colour_manual(values = c("NPI-sensitive" = "#D73027", "Control" = "#4575B4"),
                      name = NULL) +
  scale_fill_manual(values = c("NPI-sensitive" = "#D73027", "Control" = "#4575B4"),
                    name = NULL) +
  scale_x_continuous(breaks = 2015:2025) +
  labs(title = "B. Observed/Expected ratios by group",
       x = "Year", y = "Mean O/E ratio (\u00b1 SE)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

fig8 <- p8a / p8b + plot_annotation(
  title = "Figure 8. Difference-in-Differences: NPI-Sensitive vs Control Diseases",
  caption = paste0(
    "DiD coefficient (covid \u00d7 NPI-sensitive): ",
    round(did_coef$estimate, 3), " (p = ",
    round(did_coef$p.value, 4), "). Pre-trend p = ",
    round(pretrend_p, 3), "; Placebo p = ", round(placebo_p, 3), "."
  ),
  theme = theme(
    plot.title   = element_text(face = "bold", size = 13),
    plot.caption = element_text(size = 8, hjust = 0)
  )
)

ggsave(file.path(figdir, "fig8_did_panel.png"), fig8,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig8_did_panel.pdf"), fig8,
       width = 10, height = 10, bg = "white")
cat("Saved Figure 8\n\n")

# =============================================================================
# PART B: Trajectory Clustering
# =============================================================================

cat("--- Part B: Trajectory Clustering ---\n\n")

# Build rate matrix: 47 diseases × 11 years
rate_wide <- national %>%
  select(disease, year, rate_per_100k) %>%
  pivot_wider(names_from = year, values_from = rate_per_100k) %>%
  as.data.frame()

rownames(rate_wide) <- rate_wide$disease
rate_matrix <- as.matrix(rate_wide[, -1])

cat("Rate matrix:", nrow(rate_matrix), "diseases x", ncol(rate_matrix), "years\n")

# Impute rare NAs with row mean
n_missing <- sum(is.na(rate_matrix))
if (n_missing > 0) {
  cat("Imputing", n_missing, "missing values with row means\n")
  for (i in seq_len(nrow(rate_matrix))) {
    na_cols <- is.na(rate_matrix[i, ])
    if (any(na_cols)) {
      rate_matrix[i, na_cols] <- mean(rate_matrix[i, !na_cols], na.rm = TRUE)
    }
  }
}

# Z-score within each disease (row)
rate_z_matrix <- t(apply(rate_matrix, 1, function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}))
colnames(rate_z_matrix) <- colnames(rate_matrix)

# --- Determine optimal k ---
cat("Testing k = 2 to 6...\n")
sil_widths <- numeric(5)
for (k in 2:6) {
  pam_fit <- pam(rate_z_matrix, k = k, metric = "euclidean")
  sil_widths[k - 1] <- pam_fit$silinfo$avg.width
  cat("  k =", k, "-> avg silhouette =", round(sil_widths[k - 1], 3), "\n")
}

# Gap statistic
set.seed(42)
gap_stat <- clusGap(rate_z_matrix, FUNcluster = pam, K.max = 6, B = 100,
                    verbose = FALSE)
gap_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"],
               method = "Tibs2001SEmax")
cat("\nGap statistic optimal k:", gap_k, "\n")

sil_k <- which.max(sil_widths) + 1
cat("Max silhouette k:", sil_k, "(avg sil =", round(max(sil_widths), 3), ")\n")

# Use silhouette-based k (standard in applied clustering)
optimal_k <- sil_k
cat("Selected k:", optimal_k, "\n\n")

# --- Final clustering ---
set.seed(42)
final_pam <- pam(rate_z_matrix, k = optimal_k, metric = "euclidean")

cluster_df <- tibble(
  disease    = rownames(rate_z_matrix),
  cluster    = final_pam$clustering,
  silhouette = final_pam$silinfo$widths[
    match(rownames(rate_z_matrix), rownames(final_pam$silinfo$widths)),
    "sil_width"
  ]
) %>%
  left_join(
    recovery %>% select(disease, transmission_mode, border_sensitive,
                        vaccine_preventable, covid_oe_rate, post_oe_rate,
                        trajectory_refined, immunity_debt_signal,
                        recovery_year, still_suppressed),
    by = "disease"
  )

# --- Label clusters by their mean COVID O/E (most disrupted = cluster 1) ---
cluster_order <- cluster_df %>%
  group_by(cluster) %>%
  summarise(mean_covid_oe = mean(covid_oe_rate, na.rm = TRUE)) %>%
  arrange(mean_covid_oe) %>%
  mutate(new_label = row_number())

cluster_df <- cluster_df %>%
  left_join(cluster_order %>% select(cluster, new_label), by = "cluster") %>%
  mutate(cluster = new_label) %>%
  select(-new_label)

# Update PAM clustering vector to match new labels
cluster_remap <- setNames(cluster_order$new_label, cluster_order$cluster)
final_pam$clustering <- cluster_remap[as.character(final_pam$clustering)]

# --- Characterise clusters ---
cat("=== Cluster Characteristics ===\n\n")
for (cl in sort(unique(cluster_df$cluster))) {
  cl_d <- cluster_df %>% filter(cluster == cl)
  cat("Cluster", cl, "(n =", nrow(cl_d), "):\n")
  cat("  Diseases:", paste(head(cl_d$disease, 5), collapse = ", "),
      if (nrow(cl_d) > 5) paste0("... +", nrow(cl_d) - 5, " more") else "", "\n")
  tm_tab <- sort(table(cl_d$transmission_mode), decreasing = TRUE)
  cat("  Transmission:", paste(paste0(names(tm_tab), "(", tm_tab, ")"),
                               collapse = ", "), "\n")
  cat("  Mean COVID O/E:", round(mean(cl_d$covid_oe_rate), 3), "\n")
  cat("  Mean post O/E:", round(mean(cl_d$post_oe_rate), 3), "\n")
  cat("  Border-sensitive:", sum(cl_d$border_sensitive), "/", nrow(cl_d), "\n")
  cat("  VPD:", sum(cl_d$vaccine_preventable), "/", nrow(cl_d), "\n")
  cat("  Still suppressed:", sum(cl_d$still_suppressed, na.rm = TRUE), "/",
      nrow(cl_d), "\n\n")
}

# Cross-tabulate with threshold-based trajectories
cat("Cluster vs Threshold Trajectory:\n")
print(table(Cluster = cluster_df$cluster,
            Trajectory = cluster_df$trajectory_refined))
cat("\n")

write_csv(cluster_df, file.path(outdir, "trajectory_clusters.csv"))
cat("Saved cluster assignments:", nrow(cluster_df), "diseases\n\n")

# --- Figure 9: Cluster profiles ---
cat("Generating Figure 9...\n")

profile_long <- as.data.frame(rate_z_matrix) %>%
  tibble::rownames_to_column("disease") %>%
  pivot_longer(-disease, names_to = "year", values_to = "rate_z") %>%
  mutate(year = as.integer(year)) %>%
  left_join(cluster_df %>% select(disease, cluster), by = "disease")

# Cluster labels with n
cluster_labels <- cluster_df %>%
  count(cluster) %>%
  mutate(label = paste0("Cluster ", cluster, " (n=", n, ")"))

profile_long <- profile_long %>%
  left_join(cluster_labels %>% select(cluster, label), by = "cluster")

# Cluster means ± SE
cluster_profiles <- profile_long %>%
  group_by(cluster, label, year) %>%
  summarise(mean_z = mean(rate_z), se_z = sd(rate_z) / sqrt(n()),
            .groups = "drop")

p9a <- ggplot() +
  annotate("rect", xmin = covid_start - 0.5, xmax = post_start - 0.5,
           ymin = -Inf, ymax = Inf, fill = "#FEE0D2", alpha = 0.3) +
  annotate("rect", xmin = post_start - 0.5, xmax = 2025.5,
           ymin = -Inf, ymax = Inf, fill = "#DEEBF7", alpha = 0.3) +
  geom_line(data = profile_long,
            aes(x = year, y = rate_z, group = disease),
            colour = "grey70", alpha = 0.3, linewidth = 0.3) +
  geom_ribbon(data = cluster_profiles,
              aes(x = year, ymin = mean_z - se_z, ymax = mean_z + se_z),
              fill = "#2166AC", alpha = 0.3) +
  geom_line(data = cluster_profiles,
            aes(x = year, y = mean_z), colour = "#2166AC", linewidth = 1.2) +
  geom_vline(xintercept = covid_start - 0.5, linetype = "dotted", colour = "red") +
  geom_vline(xintercept = post_start - 0.5, linetype = "dotted", colour = "blue") +
  facet_wrap(~ label, ncol = min(optimal_k, 3)) +
  scale_x_continuous(breaks = seq(2015, 2025, 2)) +
  labs(title = "A. Cluster trajectory profiles (z-scored rates)",
       x = "Year", y = "Z-scored rate") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

# Panel B: cluster composition by transmission mode
comp_data <- cluster_df %>%
  count(cluster, transmission_mode) %>%
  mutate(cluster_label = paste0("Cluster ", cluster))

p9b <- ggplot(comp_data, aes(x = cluster_label, y = n, fill = transmission_mode)) +
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Set2", name = "Transmission mode") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "B. Cluster composition by transmission mode",
       x = NULL, y = "Proportion") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right")

fig9 <- p9a / p9b + plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = paste0("Figure 9. Trajectory Clustering (k=", optimal_k,
                   ", avg silhouette=", round(final_pam$silinfo$avg.width, 2), ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(figdir, "fig9_cluster_profiles.png"), fig9,
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig9_cluster_profiles.pdf"), fig9,
       width = 12, height = 10, bg = "white")
cat("Saved Figure 9\n")

cat("\n=== Script 4 complete ===\n")
