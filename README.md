# Disruption and Recovery of Notifiable Infectious Diseases After COVID-19 in Australia, 2015–2025

**Hayden Farquhar** MBBS MPHTM
ORCID: [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)

## Abstract

COVID-19 pandemic control measures profoundly disrupted transmission of other infectious diseases worldwide. We analysed notification data for 47 nationally notifiable diseases in Australia from 2015 to 2025, quantifying disruption during the acute pandemic period (2020–2021) and characterising heterogeneous recovery trajectories through 2025. Twenty-eight diseases decreased during the acute period, with border-sensitive diseases (difference-in-differences: −0.50, p = 0.016) and vaccine-preventable diseases showing the greatest suppression. By 2025, 17 diseases had overshot pre-pandemic baselines, 12 had returned to baseline, 9 remained partially suppressed, and 6 showed sustained suppression. Five diseases exhibited immunity debt signals, with influenza accumulating a cumulative excess of 739,000 notifications above expected levels post-pandemic. Trajectory clustering identified two distinct response patterns, with cluster membership strongly predicting time-to-recovery (log-rank p = 0.004). Twelve diseases remain below baseline as of 2025, of which six show no trend toward recovery.

## Data Sources

| Source | Description | Access |
|--------|-------------|--------|
| [NNDSS](https://nndss.health.gov.au) | Annual and fortnightly aggregate notifications, all states, 1991–2025 | Public (Power BI dashboard) |
| [cdc.gov.au](https://www.cdc.gov.au) | Record-level data for influenza, meningococcal, salmonellosis | Public (XLSX downloads) |
| [OxCGRT](https://github.com/OxCGRT/covid-policy-dataset) | Government Response Tracker — stringency indices | CC BY 4.0 |
| [ABS](https://www.abs.gov.au) | Estimated Resident Population (3101.0 Table 4) | CC BY 4.0 |

See [`data/raw/README.md`](data/raw/README.md) for step-by-step download instructions.

## Requirements

- **R** ≥ 4.4.1
- **26 R packages** — install all at once:

```r
source("code/00-setup/install_dependencies.R")
```

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/hayden-farquhar/NNDSS-COVID-disruption.git
cd NNDSS-COVID-disruption

# 2. Install R dependencies
Rscript code/00-setup/install_dependencies.R

# 3. Run any analysis script (processed data is included)
Rscript code/03-analysis/01_tier2_annual_baselines.R
```

Processed data is committed to `data/processed/`, so you can run analysis, visualisation, and modelling scripts without downloading raw data. To replicate from scratch, follow the [Replication Guide](docs/REPLICATION_GUIDE.md).

## Analysis Pipeline

| Phase | Directory | Scripts | Description |
|-------|-----------|---------|-------------|
| 0 | `code/00-setup/` | 1 | Install R dependencies |
| 1 | `code/01-acquisition/` | 3 | Download and parse NNDSS data from Power BI |
| 2 | `code/02-cleaning/` | 3 | Quality checks, disease classification, analysis datasets |
| 3 | `code/03-analysis/` | 10 | Baselines, ARIMA, disruption metrics, recovery, immunity debt |
| 4 | `code/04-visualisation/` | 3 | Publication figures (13 figures, PNG + PDF) |
| 5 | `code/05-modelling/` | 5 | ITS, DiD, clustering, survival, GAM, bootstrap |

Scripts are numbered sequentially within each phase and should be run in order. See [`docs/REPLICATION_GUIDE.md`](docs/REPLICATION_GUIDE.md) for details.

## Key Results

- **28/47** diseases decreased during the acute COVID period (2020–2021)
- **Border-sensitive** (DiD: −0.50, p = 0.016) and **vaccine-preventable** diseases were most disrupted
- Most disrupted: influenza (O/E 0.065), dengue (0.076), pertussis (0.115)
- **17** diseases overshot baseline by 2025; **12** still below baseline
- **5** diseases show immunity debt signals (suppression followed by overshoot)
- Median time-to-recovery for suppressed diseases: 5 years
- Trajectory cluster membership predicts recovery time (log-rank p = 0.004)

## Repository Structure

```
├── code/
│   ├── 00-setup/              # Package installer
│   ├── 01-acquisition/        # Data download and parsing
│   ├── 02-cleaning/           # Validation and classification
│   ├── 03-analysis/           # Baselines, disruption, recovery
│   ├── 04-visualisation/      # Publication figures
│   └── 05-modelling/          # ITS, DiD, clustering, survival
├── data/
│   ├── raw/                   # Download instructions (data not included)
│   └── processed/             # 48 analysis-ready CSVs (included)
├── output/
│   ├── figures/               # 13 figures (PNG + PDF)
│   └── tables/                # 11 summary tables
└── docs/
    └── REPLICATION_GUIDE.md   # Full replication walkthrough
```

## Statistical Methods

1. **Disruption quantification** — Observed/expected ratios against 2015–2019 baselines with prediction intervals
2. **Interrupted time series** — Newey-West HAC standard errors (annual, n = 11); GLS with AR(1) correlation (monthly)
3. **Difference-in-differences** — Border-sensitive vs non-border-sensitive diseases as natural experiment
4. **Predictor regression** — Nested models identifying border sensitivity and VPD status as strongest predictors (R² = 0.41)
5. **Trajectory clustering** — k-medoids on O/E ratio profiles (k = 2, silhouette = 0.235)
6. **Survival analysis** — Kaplan-Meier and Cox proportional hazards for time-to-recovery
7. **Bootstrap immunity debt** — 999 resamples for cumulative deficit/excess confidence intervals

## Ethics and Data

This study used publicly available aggregate surveillance data from the Australian National Notifiable Diseases Surveillance System. No ethics approval was required. No individual-level data were accessed.

## License

Code: [MIT](LICENSE)

Data: Processed datasets derived from publicly available Australian Government surveillance data. Original data sources retain their respective licenses (see [data sources](#data-sources)).

## Citation

If you use this code or data, please cite:

> Farquhar H. Disruption and recovery of notifiable infectious diseases after COVID-19 in Australia, 2015–2025. 2026. Available from: https://github.com/hayden-farquhar/NNDSS-COVID-disruption

See also [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.
