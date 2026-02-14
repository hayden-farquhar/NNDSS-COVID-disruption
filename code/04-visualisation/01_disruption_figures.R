#!/usr/bin/env Rscript
# =============================================================================
# 01_disruption_figures.R
#
# Publication-quality figures for the disruption analysis:
#   Figure 1: Heatmap of O/E ratios by disease and year
#   Figure 2: Forest plot of COVID disruption by transmission mode
#   Figure 3: Time series panels for selected diseases (Tier 1 monthly)
#   Figure 4: State-level disruption comparison
#   Figure 5: Recovery trajectory summary
#
# Inputs:
#   data/processed/disruption_metrics_national.csv
#   data/processed/disruption_year_by_year.csv
#   data/processed/annual_with_baselines.csv
#   data/processed/disruption_by_state.csv
#   data/processed/tier1_arima_forecasts.csv
#   data/processed/disruption_monthly_tier1.csv
#
# Outputs:
#   output/figures/fig1_heatmap_disruption.png/.pdf
#   output/figures/fig2_forest_plot.png/.pdf
#   output/figures/fig3_time_series_panels.png/.pdf
#   output/figures/fig4_state_comparison.png/.pdf
#   output/figures/fig5_recovery_trajectories.png/.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(patchwork)
  library(scales)
})

outdir <- "data/processed"
figdir <- "output/figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

cat("=== Generating Disruption Figures ===\n\n")

# Consistent colour palette for transmission modes
tm_colours <- c(
  "respiratory"   = "#E41A1C",
  "enteric"       = "#377EB8",
  "vector_borne"  = "#4DAF4A",
  "STI"           = "#984EA3",
  "blood_borne"   = "#FF7F00",
  "zoonotic"      = "#A65628"
)

# Period colours
period_colours <- c(
  "baseline"    = "#999999",
  "covid_acute" = "#E41A1C",
  "transition"  = "#FF7F00",
  "post_covid"  = "#377EB8"
)

# Trajectory colours
traj_colours <- c(
  "overshoot"             = "#E41A1C",
  "returned"              = "#4DAF4A",
  "partial_recovery"      = "#FF7F00",
  "sustained_suppression" = "#984EA3",
  "insufficient_data"     = "#999999"
)

# Nice labels
tm_labels <- c(
  "respiratory"  = "Respiratory",
  "enteric"      = "Enteric",
  "vector_borne" = "Vector-borne",
  "STI"          = "STI",
  "blood_borne"  = "Blood-borne",
  "zoonotic"     = "Zoonotic"
)

# Load data
metrics   <- read_csv(file.path(outdir, "disruption_metrics_national.csv"),
                      show_col_types = FALSE)
yby       <- read_csv(file.path(outdir, "disruption_year_by_year.csv"),
                      show_col_types = FALSE)
abl       <- read_csv(file.path(outdir, "annual_with_baselines.csv"),
                      show_col_types = FALSE)
state_dis <- read_csv(file.path(outdir, "disruption_by_state.csv"),
                      show_col_types = FALSE)
arima_fc  <- read_csv(file.path(outdir, "tier1_arima_forecasts.csv"),
                      show_col_types = FALSE)

# Shorten long disease names for figures
short_name <- function(x) {
  x %>%
    str_replace("Influenza \\(laboratory confirmed\\)", "Influenza") %>%
    str_replace("Meningococcal disease \\(invasive\\)", "Meningococcal") %>%
    str_replace("Pneumococcal disease \\(invasive\\)", "Pneumococcal") %>%
    str_replace("Haemophilus influenzae type b", "Hib") %>%
    str_replace("Haemolytic uraemic syndrome \\(HUS\\)", "HUS") %>%
    str_replace("Varicella zoster \\(chickenpox\\)", "Varicella (chickenpox)") %>%
    str_replace("Varicella zoster \\(shingles\\)", "Varicella (shingles)") %>%
    str_replace("Varicella zoster \\(unspecified\\)", "Varicella (unspec.)") %>%
    str_replace("Syphilis < 2 years", "Syphilis (<2yr)") %>%
    str_replace("Syphilis > 2 years or unspecified duration", "Syphilis (>2yr)") %>%
    str_replace("Syphilis congenital", "Congenital syphilis") %>%
    str_replace("Hepatitis B \\(newly acquired\\)", "Hep B (new)") %>%
    str_replace("Hepatitis B \\(unspecified\\)", "Hep B (unspec.)") %>%
    str_replace("Hepatitis C \\(newly acquired\\)", "Hep C (new)") %>%
    str_replace("Hepatitis C \\(unspecified\\)", "Hep C (unspec.)") %>%
    str_replace("Chikungunya virus infection", "Chikungunya") %>%
    str_replace("Dengue virus infection", "Dengue") %>%
    str_replace("Ross River virus infection", "Ross River") %>%
    str_replace("Barmah Forest virus infection", "Barmah Forest") %>%
    str_replace("Flavivirus infection \\(unspecified\\)", "Flavivirus (unspec.)") %>%
    str_replace("Chlamydial infection", "Chlamydia") %>%
    str_replace("Gonococcal infection", "Gonorrhoea") %>%
    str_replace("Campylobacteriosis", "Campylobacter") %>%
    str_replace("Cryptosporidiosis", "Cryptosporidium") %>%
    str_replace("Salmonellosis", "Salmonella") %>%
    str_replace("Legionellosis", "Legionella") %>%
    str_replace("Leptospirosis", "Leptospira")
}

# =============================================================================
# Figure 1: Heatmap of O/E ratios by disease and year
# =============================================================================
cat("--- Figure 1: Heatmap ---\n")

# Use year-by-year data plus baseline (O/E = 1 by definition during baseline)
heatmap_data <- abl %>%
  filter(level == "national") %>%
  select(disease, year, period, transmission_mode, oe_ratio, count_tier) %>%
  filter(count_tier != "ultra_low") %>%
  mutate(
    disease_short = short_name(disease),
    # Cap O/E for colour scale
    oe_capped = pmin(pmax(oe_ratio, 0), 3)
  )

# Order diseases by transmission mode, then by COVID O/E
disease_order <- metrics %>%
  filter(count_tier != "ultra_low") %>%
  arrange(transmission_mode, covid_oe_rate) %>%
  mutate(disease_short = short_name(disease)) %>%
  pull(disease_short)

heatmap_data <- heatmap_data %>%
  mutate(disease_short = factor(disease_short, levels = rev(disease_order)))

fig1 <- ggplot(heatmap_data, aes(x = factor(year), y = disease_short,
                                  fill = oe_capped)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 1, limits = c(0, 3),
    name = "O/E Ratio",
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3),
    labels = c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "\u22653.0")
  ) +
  geom_vline(xintercept = 5.5, linetype = "dashed", colour = "black",
             linewidth = 0.5) +
  annotate("text", x = 3, y = 0.5, label = "Baseline", size = 3,
           colour = "grey40") +
  annotate("text", x = 7, y = 0.5, label = "COVID/Post", size = 3,
           colour = "grey40") +
  labs(
    x = "Year",
    y = NULL,
    title = "Observed/Expected Ratios for Notifiable Diseases, 2015\u20132025"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold")
  )

ggsave(file.path(figdir, "fig1_heatmap_disruption.png"), fig1,
       width = 10, height = 14, dpi = 300)
ggsave(file.path(figdir, "fig1_heatmap_disruption.pdf"), fig1,
       width = 10, height = 14)
cat("Saved fig1_heatmap_disruption.png/.pdf\n")

# =============================================================================
# Figure 2: Forest plot of COVID disruption by transmission mode
# =============================================================================
cat("--- Figure 2: Forest plot ---\n")

forest_data <- metrics %>%
  filter(count_tier != "ultra_low") %>%
  mutate(
    disease_short = short_name(disease),
    tm_label = tm_labels[transmission_mode]
  ) %>%
  arrange(transmission_mode, covid_oe_rate) %>%
  mutate(disease_short = factor(disease_short,
                                levels = rev(disease_short)))

fig2 <- ggplot(forest_data,
               aes(x = covid_oe_rate, y = disease_short,
                   colour = transmission_mode)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5) +
  # Add post-COVID O/E as open circle
  geom_point(aes(x = post_oe_rate), shape = 1, size = 2.5) +
  # Connecting line
  geom_segment(aes(x = covid_oe_rate, xend = post_oe_rate,
                   y = disease_short, yend = disease_short),
               linewidth = 0.3, alpha = 0.5) +
  scale_colour_manual(values = tm_colours, labels = tm_labels,
                      name = "Transmission Mode") +
  scale_x_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3.2)) +
  labs(
    x = "Observed/Expected Ratio",
    y = NULL,
    title = "COVID-19 Disruption and Recovery by Disease",
    subtitle = "Filled = COVID-acute (2020\u20132021); Open = Post-COVID (2023\u20132025)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold")
  )

ggsave(file.path(figdir, "fig2_forest_plot.png"), fig2,
       width = 10, height = 12, dpi = 300)
ggsave(file.path(figdir, "fig2_forest_plot.pdf"), fig2,
       width = 10, height = 12)
cat("Saved fig2_forest_plot.png/.pdf\n")

# =============================================================================
# Figure 3: Time series panels for Tier 1 diseases (monthly ARIMA)
# =============================================================================
cat("--- Figure 3: Time series panels ---\n")

ts_data <- arima_fc %>%
  mutate(
    date = as.Date(paste(year, month, 1, sep = "-")),
    disease_short = short_name(disease)
  )

fig3 <- ggplot(ts_data, aes(x = date)) +
  # Forecast ribbon
  geom_ribbon(aes(ymin = expected_lo95, ymax = expected_hi95),
              fill = "steelblue", alpha = 0.15, na.rm = TRUE) +
  # Expected line
  geom_line(aes(y = expected), colour = "steelblue", linewidth = 0.5,
            linetype = "dashed", na.rm = TRUE) +
  # Observed line
  geom_line(aes(y = observed), colour = "black", linewidth = 0.4,
            na.rm = TRUE) +
  # COVID onset marker
  geom_vline(xintercept = as.Date("2020-03-01"),
             linetype = "dotted", colour = "red", linewidth = 0.5) +
  # Border reopening marker
  geom_vline(xintercept = as.Date("2022-03-01"),
             linetype = "dotted", colour = "blue", linewidth = 0.5) +
  facet_wrap(~ disease_short, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Date",
    y = "Monthly Count",
    title = "Monthly Time Series with ARIMA Baseline Forecasts",
    subtitle = "Black = observed; Blue dashed = expected; Shading = 95% PI\nRed line = COVID onset (Mar 2020); Blue line = borders reopen (Mar 2022)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold")
  )

ggsave(file.path(figdir, "fig3_time_series_panels.png"), fig3,
       width = 10, height = 10, dpi = 300)
ggsave(file.path(figdir, "fig3_time_series_panels.pdf"), fig3,
       width = 10, height = 10)
cat("Saved fig3_time_series_panels.png/.pdf\n")

# =============================================================================
# Figure 4: State-level disruption comparison
# =============================================================================
cat("--- Figure 4: State comparison ---\n")

# Aggregate state-level COVID O/E by transmission mode
state_tm <- state_dis %>%
  filter(!is.na(covid_oe)) %>%
  group_by(state, transmission_mode) %>%
  summarise(
    mean_covid_oe = mean(covid_oe, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(tm_label = tm_labels[transmission_mode])

# Order states by overall disruption (most disrupted first)
state_order <- state_dis %>%
  group_by(state) %>%
  summarise(mean_oe = mean(covid_oe, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_oe) %>%
  pull(state)

state_tm <- state_tm %>%
  mutate(state = factor(state, levels = state_order))

fig4 <- ggplot(state_tm,
               aes(x = state, y = mean_covid_oe, fill = transmission_mode)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = tm_colours, labels = tm_labels,
                    name = "Transmission Mode") +
  scale_y_continuous(breaks = seq(0, 3, 0.5), limits = c(0, NA)) +
  labs(
    x = "State/Territory",
    y = "Mean COVID-Acute O/E Ratio",
    title = "State-Level Disease Disruption During COVID-19 (2020\u20132021)",
    subtitle = "By transmission mode; dashed line = baseline (O/E = 1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold")
  )

ggsave(file.path(figdir, "fig4_state_comparison.png"), fig4,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(figdir, "fig4_state_comparison.pdf"), fig4,
       width = 10, height = 6)
cat("Saved fig4_state_comparison.png/.pdf\n")

# =============================================================================
# Figure 5: Recovery trajectory summary
# =============================================================================
cat("--- Figure 5: Recovery trajectories ---\n")

# Stacked bar by transmission mode
traj_data <- metrics %>%
  filter(count_tier != "ultra_low") %>%
  mutate(
    tm_label = tm_labels[transmission_mode],
    trajectory = factor(trajectory,
                        levels = c("overshoot", "returned",
                                   "partial_recovery", "sustained_suppression"))
  ) %>%
  count(tm_label, trajectory)

fig5a <- ggplot(traj_data, aes(x = tm_label, y = n, fill = trajectory)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(
    values = traj_colours,
    labels = c("Overshoot (>120%)", "Returned (90\u2013120%)",
               "Partial (50\u201390%)", "Sustained suppression (<50%)"),
    name = "Recovery Status"
  ) +
  labs(
    x = "Transmission Mode",
    y = "Number of Diseases",
    title = "Post-COVID Recovery Trajectory by Transmission Mode"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 12, face = "bold")
  )

# Scatter: COVID O/E vs Post-COVID O/E
scatter_data <- metrics %>%
  filter(count_tier != "ultra_low") %>%
  mutate(disease_short = short_name(disease))

fig5b <- ggplot(scatter_data,
                aes(x = covid_oe_rate, y = post_oe_rate,
                    colour = transmission_mode)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  # Quadrant shading
  annotate("rect", xmin = 0, xmax = 1, ymin = 1.2, ymax = Inf,
           fill = "#E41A1C", alpha = 0.05) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 0.5,
           fill = "#984EA3", alpha = 0.05) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(
    data = scatter_data %>%
      filter(covid_oe_rate < 0.2 | post_oe_rate > 2 |
               post_oe_rate < 0.4 | covid_oe_rate > 2),
    aes(label = disease_short),
    size = 2.5, max.overlaps = 15, show.legend = FALSE
  ) +
  scale_colour_manual(values = tm_colours, labels = tm_labels,
                      name = "Transmission Mode") +
  scale_x_continuous(limits = c(0, 3)) +
  scale_y_continuous(limits = c(0, 3.5)) +
  labs(
    x = "COVID-Acute O/E Ratio (2020\u20132021)",
    y = "Post-COVID O/E Ratio (2023\u20132025)",
    title = "Disruption vs Recovery"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold")
  )

fig5 <- fig5a + fig5b + plot_layout(widths = c(1, 1.3))

ggsave(file.path(figdir, "fig5_recovery_trajectories.png"), fig5,
       width = 14, height = 7, dpi = 300)
ggsave(file.path(figdir, "fig5_recovery_trajectories.pdf"), fig5,
       width = 14, height = 7)
cat("Saved fig5_recovery_trajectories.png/.pdf\n")

# =============================================================================
# Additional: Annual time series for all diseases (supplementary)
# =============================================================================
cat("--- Supplementary: Annual time series grid ---\n")

annual_ts <- abl %>%
  filter(level == "national", count_tier != "ultra_low") %>%
  mutate(disease_short = short_name(disease))

fig_supp <- ggplot(annual_ts, aes(x = year, y = rate_per_100k)) +
  geom_ribbon(aes(ymin = pi_lower_rate, ymax = pi_upper_rate),
              fill = "steelblue", alpha = 0.15) +
  geom_hline(aes(yintercept = baseline_mean_rate),
             linetype = "dashed", colour = "steelblue", linewidth = 0.3) +
  geom_line(linewidth = 0.4) +
  geom_point(aes(colour = period), size = 1) +
  scale_colour_manual(values = period_colours, name = "Period") +
  geom_vline(xintercept = 2019.5, linetype = "dotted",
             colour = "red", linewidth = 0.3) +
  facet_wrap(~ disease_short, scales = "free_y", ncol = 6) +
  labs(
    x = "Year",
    y = "Rate per 100,000",
    title = "Annual Notification Rates for 44 Notifiable Diseases, 2015\u20132025",
    subtitle = "Dashed line = baseline mean; Shading = 95% prediction interval; Red line = COVID onset"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    strip.text = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold")
  )

ggsave(file.path(figdir, "fig_supp_annual_timeseries.png"), fig_supp,
       width = 16, height = 18, dpi = 300)
ggsave(file.path(figdir, "fig_supp_annual_timeseries.pdf"), fig_supp,
       width = 16, height = 18)
cat("Saved fig_supp_annual_timeseries.png/.pdf\n")

cat("\n=== All figures generated ===\n")
