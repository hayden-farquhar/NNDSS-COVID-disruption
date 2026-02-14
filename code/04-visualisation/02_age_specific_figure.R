# ============================================================================
# Figure S8: Age-specific disruption and recovery for influenza and
#            pneumococcal disease
# Supports policy claim: children 0-4 bear disproportionate rebound burden
# ============================================================================

library(tidyverse)

# --- Load data ---------------------------------------------------------------
age_data <- read_csv("data/processed/age_specific_disruption.csv",
                     show_col_types = FALSE)

# Focus on influenza and pneumococcal; exclude "Other" age group
# Exclude 2025 for influenza (record-level data incomplete for age breakdown)
plot_data <- age_data %>%
  filter(
    disease %in% c("Influenza (laboratory confirmed)",
                    "Pneumococcal disease (invasive)"),
    age_broad != "Other"
  ) %>%
  # Exclude 2025 influenza age data (incomplete record-level extract)
  filter(!(disease == "Influenza (laboratory confirmed)" & year == 2025)) %>%
  mutate(
    disease_short = case_when(
      str_detect(disease, "Influenza") ~ "Influenza",
      str_detect(disease, "Pneumococcal") ~ "Pneumococcal disease",
      TRUE ~ disease
    ),
    age_broad = factor(age_broad,
                       levels = c("0-4", "5-14", "15-24", "25-49", "50-64", "65+"))
  )

# --- Colour palette (age groups) -------------------------------------------
age_colours <- c(
  "0-4"   = "#E41A1C",   # red
  "5-14"  = "#FF7F00",   # orange
  "15-24" = "#984EA3",   # purple
  "25-49" = "#377EB8",   # blue
  "50-64" = "#4DAF4A",   # green
  "65+"   = "#A65628"    # brown
)

# --- Panel A: O/E ratio time series by age group ----------------------------
p_timeseries <- ggplot(plot_data, aes(x = year, y = oe_ratio,
                                       colour = age_broad, group = age_broad)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  # COVID shading
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.3) +
  facet_wrap(~ disease_short, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = age_colours, name = "Age group") +
  scale_x_continuous(breaks = seq(2015, 2024, 2)) +
  labs(
    x = "Year",
    y = "Observed / Expected ratio",
    title = "Age-specific disruption and recovery",
    subtitle = "O/E ratios by age group relative to 2015\u20132019 baselines; grey band = COVID-acute period"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40")
  )

# --- Panel B: Post-COVID mean O/E by age group (bar chart) -----------------
post_covid_summary <- plot_data %>%
  filter(period == "post_covid") %>%
  group_by(disease_short, age_broad) %>%
  summarise(
    mean_oe = mean(oe_ratio, na.rm = TRUE),
    se_oe = sd(oe_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_bars <- ggplot(post_covid_summary, aes(x = age_broad, y = mean_oe, fill = age_broad)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_oe - se_oe, ymax = mean_oe + se_oe),
                width = 0.2, linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_wrap(~ disease_short, ncol = 2) +
  scale_fill_manual(values = age_colours, guide = "none") +
  labs(
    x = "Age group",
    y = "Mean post-COVID O/E ratio",
    title = "Post-pandemic rebound by age group",
    subtitle = "Mean O/E ratio during post-COVID period (2023\u20132024); dashed line = baseline"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40")
  )

# --- Combine panels ----------------------------------------------------------
library(patchwork)

p_combined <- p_timeseries / p_bars +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# --- Save --------------------------------------------------------------------
ggsave("output/figures/fig_supp_age_specific.png", p_combined,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave("output/figures/fig_supp_age_specific.pdf", p_combined,
       width = 10, height = 10, bg = "white")

cat("Figure S8 saved: output/figures/fig_supp_age_specific.png/.pdf\n")
