# Raw Data Download Instructions

Raw data files are not included in this repository due to size and redistribution constraints. Follow the instructions below to obtain them. Each data source is free and publicly accessible.

> **Note:** If you only want to run the analysis, visualisation, and modelling scripts, you do not need raw data — all processed datasets are included in `data/processed/`.

## 1. NNDSS Power BI Dashboard

**Source:** https://nndss.health.gov.au

The NNDSS Power BI dashboard provides aggregate annual and fortnightly notification data for all nationally notifiable diseases.

### Annual data (rates by state and year)

1. Navigate to https://nndss.health.gov.au
2. Select **"Annual notifications and notification rate tables"**
3. For each disease of interest, export the table as XLSX or CSV
4. Save files to `data/raw/powerbi/`

### Fortnightly data

1. On the same dashboard, select **"Fortnightly reports"**
2. Export data for the relevant fortnightly periods
3. Save to `data/raw/fortnightly/`

**Expected files:**
- `data/raw/powerbi/*.xlsx` — One file per disease or disease group
- `data/raw/fortnightly/*.csv` — Fortnightly aggregates

## 2. Record-Level Data (cdc.gov.au)

**Source:** https://www.cdc.gov.au

Record-level notification data is available for three diseases with detailed demographic and temporal fields.

### Influenza

1. Visit https://www.cdc.gov.au/topics/communicable-diseases/influenza-flu/influenza-surveillance-reports
2. Download the surveillance XLSX file
3. Save as `data/raw/record_level/influenza_notifications.xlsx`

### Meningococcal disease

1. Visit https://www.cdc.gov.au/topics/communicable-diseases/meningococcal-disease/meningococcal-disease-surveillance-reports
2. Download the surveillance XLSX file
3. Save as `data/raw/record_level/meningococcal_notifications.xlsx`

### Salmonellosis

1. Visit https://www.cdc.gov.au/topics/communicable-diseases/salmonellosis/salmonellosis-surveillance-reports
2. Download the surveillance XLSX file
3. Save as `data/raw/record_level/salmonellosis_notifications.xlsx`

**Note:** `cdc.gov.au` may return HTTP/2 errors. If downloading programmatically, use `curl --http1.1`.

## 3. Oxford COVID-19 Government Response Tracker (OxCGRT)

**Source:** https://github.com/OxCGRT/covid-policy-dataset

Used for subnational stringency index analysis.

1. Clone or download from the OxCGRT GitHub repository
2. From the `data/` directory, obtain:
   - National time series → save as `data/raw/contextual/oxcgrt_national.csv`
   - Subnational time series → save as `data/raw/contextual/oxcgrt_subnational.csv`

These files are licensed under CC BY 4.0.

## 4. ABS Population Estimates

**Source:** https://www.abs.gov.au

Used for population denominators and rate calculations.

1. Visit ABS catalogue 3101.0 — Australian Demographic Statistics
2. Download **Table 4: Estimated Resident Population, States and Territories**
3. Save as `data/raw/contextual/abs_population_3101_table4.xlsx`

## Directory Structure After Download

```
data/raw/
├── contextual/
│   ├── oxcgrt_national.csv
│   ├── oxcgrt_subnational.csv
│   └── abs_population_3101_table4.xlsx
├── powerbi/
│   └── *.xlsx
├── record_level/
│   ├── influenza_notifications.xlsx
│   ├── meningococcal_notifications.xlsx
│   └── salmonellosis_notifications.xlsx
└── fortnightly/
    └── *.csv
```
