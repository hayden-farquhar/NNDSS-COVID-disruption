# Replication Guide

This guide explains how to reproduce the analysis in full or in part.

## Partial Replication (Recommended)

All processed datasets are included in `data/processed/`. You can run any script from Phase 3 onward without downloading raw data.

### Steps

1. **Install dependencies**
   ```bash
   Rscript code/00-setup/install_dependencies.R
   ```

2. **Run analysis scripts** (Phase 3)
   ```bash
   Rscript code/03-analysis/01_tier2_annual_baselines.R
   Rscript code/03-analysis/05_disruption_metrics.R
   Rscript code/03-analysis/07_recovery_trajectories.R
   Rscript code/03-analysis/08_immunity_debt.R
   ```

3. **Generate figures** (Phase 4)
   ```bash
   Rscript code/04-visualisation/01_disruption_figures.R
   Rscript code/04-visualisation/02_age_specific_figure.R
   Rscript code/04-visualisation/03_did_border_sensitive_figure.R
   ```

4. **Run models** (Phase 5)
   ```bash
   Rscript code/05-modelling/01_its_monthly.R
   Rscript code/05-modelling/02_its_annual.R
   Rscript code/05-modelling/03_predictors_stringency.R
   Rscript code/05-modelling/04_did_clustering.R
   Rscript code/05-modelling/05_survival_gam.R
   ```

All scripts read from `data/processed/` and write to `output/figures/` and `output/tables/`.

## Full Replication (From Raw Data)

To replicate from scratch, you must first download the raw data sources.

### 1. Download raw data

Follow the instructions in [`data/raw/README.md`](../data/raw/README.md) to obtain:

- NNDSS Power BI exports → `data/raw/powerbi/`
- Record-level XLSX files → `data/raw/record_level/`
- Fortnightly reports → `data/raw/fortnightly/`
- OxCGRT stringency data → `data/raw/contextual/`
- ABS population data → `data/raw/contextual/`

### 2. Run the complete pipeline

Execute scripts in numeric order within each phase:

```bash
# Phase 1: Acquisition
Rscript code/01-acquisition/01_download_nndss_data.R
Rscript code/01-acquisition/02_parse_rate_exports.R
Rscript code/01-acquisition/03_parse_powerbi_exports.R

# Phase 2: Cleaning
Rscript code/02-cleaning/01_data_quality_checks.R
Rscript code/02-cleaning/02_disease_classification.R
Rscript code/02-cleaning/03_create_analysis_datasets.R

# Phase 3: Analysis
Rscript code/03-analysis/01_tier2_annual_baselines.R
Rscript code/03-analysis/02_tier1_parse_record_level.R
Rscript code/03-analysis/03_tier1_arima_baselines.R
Rscript code/03-analysis/04_baseline_validation.R
Rscript code/03-analysis/05_disruption_metrics.R
Rscript code/03-analysis/06_transmission_mode_state.R
Rscript code/03-analysis/07_recovery_trajectories.R
Rscript code/03-analysis/08_immunity_debt.R
Rscript code/03-analysis/09_age_specific_recovery.R
Rscript code/03-analysis/10_jev_case_study.R

# Phase 4: Visualisation
Rscript code/04-visualisation/01_disruption_figures.R
Rscript code/04-visualisation/02_age_specific_figure.R
Rscript code/04-visualisation/03_did_border_sensitive_figure.R

# Phase 5: Modelling
Rscript code/05-modelling/01_its_monthly.R
Rscript code/05-modelling/02_its_annual.R
Rscript code/05-modelling/03_predictors_stringency.R
Rscript code/05-modelling/04_did_clustering.R
Rscript code/05-modelling/05_survival_gam.R
```

### 3. Expected runtime

The full pipeline takes approximately 5–10 minutes on a modern machine. The ARIMA fitting (Phase 3, script 03) and bootstrap (Phase 5, script 05) are the most computationally intensive steps.

## Script Dependencies

Scripts within each phase must be run in order. Cross-phase dependencies:

```
Phase 1 (acquisition) → Phase 2 (cleaning) → Phase 3 (analysis)
                                                  ↓
                                             Phase 4 (figures)
                                                  ↓
                                             Phase 5 (models)
```

Phase 4 and Phase 5 both depend on Phase 3 outputs but are independent of each other.

## Key Input/Output Files

| Script | Key Inputs | Key Outputs |
|--------|-----------|-------------|
| `03/01_tier2_annual_baselines.R` | `nndss_annual_analysis.csv` | `baselines_annual_national.csv`, `annual_with_baselines.csv` |
| `03/05_disruption_metrics.R` | `annual_with_baselines.csv` | `disruption_metrics_national.csv` |
| `03/07_recovery_trajectories.R` | `annual_with_baselines.csv` | `recovery_trajectories.csv` |
| `03/08_immunity_debt.R` | `annual_with_baselines.csv` | `immunity_debt_annual.csv` |
| `04/01_disruption_figures.R` | Multiple processed CSVs | `fig1–fig7` (PNG + PDF) |
| `05/02_its_annual.R` | `annual_with_baselines.csv` | `its_annual_individual.csv` |
| `05/04_did_clustering.R` | `annual_with_baselines.csv` | `did_results.csv`, `trajectory_clusters.csv` |
| `05/05_survival_gam.R` | `recovery_trajectories.csv` | `survival_recovery.csv`, `survival_cox_results.csv` |

## Troubleshooting

- **HTTP/2 errors from cdc.gov.au**: Use `curl --http1.1` (handled automatically in acquisition scripts)
- **Power BI export formatting**: Commas in numbers are stripped automatically; blank cells are treated as zero
- **Singular fits in mixed models**: Expected with z-scored data — random intercept variance converges to zero. Results are interpretable as fixed-effects models.
- **ARIMA convergence**: If `auto.arima()` selects a different model specification on your system, results may differ slightly from published values
