# ============================================================================
# Figure 2 (revised): Border-sensitive DiD panel
# Replaces the NPI-sensitive version as the main manuscript figure
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

outdir <- "data/processed"
figdir <- "output/figures"

# --- Load data ---------------------------------------------------------------
abl <- read_csv(file.path(outdir, "annual_with_baselines.csv"), show_col_types = FALSE)
recovery <- read_csv(file.path(outdir, "recovery_trajectories.csv"), show_col_types = FALSE)

analysis_diseases <- recovery$disease

national <- abl %>%
  filter(level == "national",
         disease %in% analysis_diseases,
         year >= 2015, year <= 2025)

# Z-score rates within each disease
did_data <- national %>%
  group_by(disease) %>%
  mutate(
    rate_z = (rate_per_100k - mean(rate_per_100k, na.rm = TRUE)) /
              sd(rate_per_100k, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(rate_z = if_else(is.nan(rate_z) | is.infinite(rate_z), 0, rate_z))

# --- Group means: border-sensitive vs non-border-sensitive -------------------
n_border <- n_distinct(did_data$disease[did_data$border_sensitive == 1])
n_control <- n_distinct(did_data$disease[did_data$border_sensitive == 0])

group_means <- did_data %>%
  mutate(group = if_else(border_sensitive == 1, "Border-sensitive", "Non-border-sensitive")) %>%
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

# --- Panel A: Z-scored rates ------------------------------------------------
p_a <- ggplot(group_means, aes(x = year, y = mean_z, colour = group, fill = group)) +
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
  scale_colour_manual(values = c("Border-sensitive" = "#D73027",
                                  "Non-border-sensitive" = "#4575B4"),
                      name = NULL) +
  scale_fill_manual(values = c("Border-sensitive" = "#D73027",
                                "Non-border-sensitive" = "#4575B4"),
                    name = NULL) +
  scale_x_continuous(breaks = 2015:2025) +
  labs(title = "A. Z-scored notification rates",
       subtitle = paste0("Border-sensitive (n=", n_border,
                         ") vs non-border-sensitive (n=", n_control, ")"),
       x = NULL, y = "Mean z-scored rate (\u00b1 SE)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# --- Panel B: O/E ratios ----------------------------------------------------
p_b <- ggplot(group_means, aes(x = year, y = mean_oe, colour = group, fill = group)) +
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
  scale_colour_manual(values = c("Border-sensitive" = "#D73027",
                                  "Non-border-sensitive" = "#4575B4"),
                      name = NULL) +
  scale_fill_manual(values = c("Border-sensitive" = "#D73027",
                                "Non-border-sensitive" = "#4575B4"),
                    name = NULL) +
  scale_x_continuous(breaks = 2015:2025) +
  labs(title = "B. Observed/Expected ratios by group",
       x = "Year", y = "Mean O/E ratio (\u00b1 SE)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# --- Combine and save --------------------------------------------------------
fig2 <- p_a / p_b + plot_annotation(
  title = "Figure 2. Difference-in-Differences: Border-Sensitive vs Non-Border-Sensitive Diseases",
  caption = "DiD coefficient (covid \u00d7 border-sensitive): \u22120.50 (95% CI \u22120.90 to \u22120.10, p = 0.016). Pre-trend test: NS.",
  theme = theme(
    plot.title   = element_text(face = "bold", size = 13),
    plot.caption = element_text(size = 8, hjust = 0)
  )
)

ggsave(file.path(figdir, "fig2_did_border_sensitive.png"), fig2,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(figdir, "fig2_did_border_sensitive.pdf"), fig2,
       width = 10, height = 10, bg = "white")

cat("Figure 2 (border-sensitive DiD) saved to output/figures/\n")
