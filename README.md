# Reproducible analysis

Run from the workspace root:

```sh
python work/analysis/00_download.py
Rscript --vanilla work/analysis/01_prepare.R
Rscript --vanilla work/analysis/02_analyze.R
Rscript --vanilla work/analysis/03_sensitivity_tables.R
Rscript --vanilla work/analysis/05_lsm_sensitivity.R
```

Requirements: R 4.6.1 with the `survey`, `dplyr`, `ggplot2`, `splines`,
`broom`, and `haven` packages; Python 3 with `requests`.

The analysis treats NHANES 2017-March 2020 pre-pandemic and August
2021-August 2023 as separate survey periods and never stacks their person
weights. The official two-day dietary weights, strata, and PSU identifiers
are used within each period.

`02_analyze.R` also reports observed agreement, PABAK, prevalence and bias
indices, positive and negative agreement, and survey-weighted Day 1-Day 2 DKR
correlations with stratified jackknife confidence intervals.

`05_lsm_sensitivity.R` is an add-on sensitivity analysis of liver stiffness
measurement >=8 kPa using the same design, standardization, covariates, and
inference as the main analysis.

Pre-registration disclosure: exposure-distribution feasibility checks were
performed before the formal analysis. No DKR-CAP or DKR-MASLD association model
had been examined at that stage.
