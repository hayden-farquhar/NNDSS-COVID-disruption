# Data Dictionary

Column definitions for the key processed datasets. All files are CSV format with headers.

## Core Datasets

### `nndss_annual_analysis.csv` (6,534 rows)

Annual notification data for 66 diseases across 9 jurisdictions, 2015–2025.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | Disease name (NNDSS standard) |
| `year` | integer | Calendar year |
| `state` | character | Jurisdiction code (ACT, NSW, NT, QLD, SA, TAS, VIC, WA, National) |
| `rate_per_100k` | numeric | Notification rate per 100,000 population |
| `count` | integer | National notification count (repeated across states) |
| `transmission_mode` | character | respiratory, enteric, vector_borne, zoonotic, STI, blood_borne |
| `subcategory` | character | Disease subgroup within transmission mode |
| `analysis_suitable` | logical | TRUE if disease meets criteria for inclusion (n = 47) |
| `border_sensitive` | logical | TRUE if disease is predominantly acquired overseas |
| `vaccine_preventable` | logical | TRUE if a vaccine is available in Australia |
| `exclusion_reason` | character | Reason for exclusion if analysis_suitable = FALSE |
| `period` | character | baseline, covid_acute, transition, post_covid |

### `annual_with_baselines.csv` (4,620 rows)

Analysis-suitable diseases with baseline statistics and observed/expected ratios.

| Column | Type | Description |
|--------|------|-------------|
| `disease` … `period` | — | Same as `nndss_annual_analysis.csv` |
| `baseline_mean_rate` | numeric | Mean rate across 2015–2019 |
| `baseline_sd_rate` | numeric | Standard deviation of baseline rate |
| `baseline_cv_rate` | numeric | Coefficient of variation (SD/mean) |
| `pi_lower_rate` / `pi_upper_rate` | numeric | 95% prediction interval bounds (rate) |
| `baseline_mean_count` | numeric | Mean national count across 2015–2019 |
| `pi_lower_count` / `pi_upper_count` | numeric | 95% prediction interval bounds (count) |
| `count_tier` | character | high (≥1000/yr), medium (100–999), low (<100) |
| `trend_significant` | logical | Significant linear trend in baseline period |
| `trend_slope` | numeric | Annual change in rate (linear regression) |
| `trend_pct_per_year` | numeric | Percentage change per year |
| `oe_ratio` | numeric | Observed/expected ratio (rate-based) |
| `pct_change` | numeric | Percentage change from baseline mean |
| `below_pi` / `above_pi` | logical | Below/above 95% prediction interval |
| `oe_ratio_count` | numeric | O/E ratio (count-based) |
| `level` | character | Disruption classification |

### `baselines_annual_national.csv` (47 rows)

One row per analysis-suitable disease with baseline summary statistics.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | Disease name |
| `baseline_n_years` | integer | Number of baseline years (typically 5) |
| `baseline_mean_rate` | numeric | Mean annual rate 2015–2019 |
| `baseline_sd_rate` | numeric | Standard deviation |
| `baseline_cv_rate` | numeric | Coefficient of variation |
| `baseline_median_rate` | numeric | Median rate |
| `baseline_min_rate` / `baseline_max_rate` | numeric | Range |
| `baseline_mean_count` … `baseline_max_count` | numeric | Same statistics for counts |
| `pi_lower_rate` / `pi_upper_rate` | numeric | 95% prediction interval (rate) |
| `pi_lower_count` / `pi_upper_count` | numeric | 95% prediction interval (count) |
| `count_tier` | character | high, medium, low |
| `trend_slope` / `trend_se` / `trend_p_value` | numeric | Linear trend parameters |
| `trend_significant` | logical | p < 0.05 |
| `trend_pct_per_year` | numeric | Annual percentage change |

### `disease_metadata.csv` (69 rows)

Classification and data availability for all 66+ diseases in the NNDSS.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | Disease name |
| `transmission_mode` | character | Primary transmission route |
| `subcategory` | character | Subgroup |
| `year_data_from` / `year_data_to` | integer | Data availability range |
| `baseline_years_available` | integer | Years with data in 2015–2019 |
| `has_full_baseline` / `has_covid_data` / `has_post_data` | logical | Data completeness flags |
| `analysis_suitable` | logical | Meets inclusion criteria |
| `rate_2019` / `rate_2020` | numeric | Key reference rates |
| `national_caveats` | character | Data quality notes |
| `border_sensitive` | logical | Predominantly imported |
| `vaccine_preventable` | logical | Vaccine available |
| `exclusion_reason` | character | Reason for exclusion |

## Recovery & Trajectory Datasets

### `recovery_trajectories.csv` (47 rows)

Recovery classification and immunity debt signals per disease.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | Disease name |
| `transmission_mode` / `border_sensitive` / `vaccine_preventable` | — | Classification |
| `count_tier` | character | Volume tier |
| `covid_oe_rate` / `covid_pct_change` | numeric | O/E ratio and % change during COVID-acute |
| `post_oe_rate` / `post_pct_change` | numeric | O/E ratio and % change post-COVID |
| `trajectory` | character | Initial trajectory assignment |
| `cumulative_deficit_total` / `cumulative_excess_total` | numeric | Total deficit/excess cases |
| `net_cumulative` | numeric | Net cumulative (excess − deficit) |
| `nadir_year` / `nadir_oe` | numeric | Year and O/E of lowest point |
| `recovery_year` | integer | Year O/E first returned to ≥0.9 (NA if not recovered) |
| `years_to_recovery` | integer | Years from 2020 to recovery |
| `recovered_by_2024` / `still_suppressed` | logical | Status flags |
| `trajectory_refined` | character | overshoot, returned, partial, sustained_suppression, etc. |
| `had_rebound_peak` | logical | Exceeded baseline after suppression |
| `immunity_debt_signal` | logical | Evidence of suppression → overshoot pattern |

### `trajectory_clusters.csv` (47 rows)

k-medoids clustering results.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | Disease name |
| `cluster` | integer | Cluster assignment (1 = heavy disruption, 2 = mild) |
| `silhouette` | numeric | Silhouette width for this observation |
| `covid_oe_rate` / `post_oe_rate` | numeric | O/E ratios used in clustering |
| `trajectory_refined` | character | Refined trajectory label |
| `immunity_debt_signal` | logical | Debt signal flag |
| `recovery_year` | integer | Year recovered (NA if not) |
| `still_suppressed` | logical | Still below baseline |

### `immunity_debt_annual.csv` (282 rows)

Year-by-year cumulative deficit and excess for each disease.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | Disease name |
| `year` | integer | Calendar year |
| `period` | character | Time period classification |
| `count` | integer | Observed national count |
| `expected_count` | numeric | Expected count (baseline mean) |
| `deficit` / `excess` | numeric | Cases below/above expected |
| `cum_deficit` / `cum_excess` / `cum_net` | numeric | Running cumulative totals |
| `debt_ratio` | numeric | Cumulative excess / cumulative deficit |

### `immunity_debt_bootstrap.csv` (47 rows)

Bootstrap confidence intervals for net cumulative case balance.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | Disease name |
| `net_point` | numeric | Point estimate of net cumulative cases |
| `net_ci_lo` / `net_ci_hi` | numeric | 95% bootstrap CI |
| `deficit_point` / `deficit_ci_lo` / `deficit_ci_hi` | numeric | Cumulative deficit with CI |
| `excess_point` / `excess_ci_lo` / `excess_ci_hi` | numeric | Cumulative excess with CI |
| `ci_excludes_zero` | logical | TRUE if CI for net does not include zero |
| `transmission_mode` | character | Transmission route |
| `immunity_debt_signal` | logical | Debt signal flag |

## Modelling Datasets

### `did_results.csv` (42 rows)

Difference-in-differences regression output.

| Column | Type | Description |
|--------|------|-------------|
| `model` | character | Model identifier (M1–M4, border_sensitive_only, etc.) |
| `term` | character | Regression term |
| `estimate` | numeric | Coefficient estimate |
| `std.error` | numeric | Standard error |
| `statistic` | numeric | t-statistic |
| `p.value` | numeric | p-value |
| `conf.low` / `conf.high` | numeric | 95% confidence interval |

### `tier1_monthly_national.csv` (612 rows)

Monthly national notification counts for three Tier 1 diseases.

| Column | Type | Description |
|--------|------|-------------|
| `disease` | character | influenza, meningococcal disease (invasive), salmonellosis |
| `year` | integer | Calendar year |
| `month` | integer | Month (1–12) |
| `count` | integer | National notification count |
