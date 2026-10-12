# Forecasting Spanish Cement Consumption with SARIMA Models

Time series analysis and forecasting of monthly apparent cement consumption in
Spain (1990–2017) in **R**, following the Box–Jenkins methodology and extending
it with calendar effects, intervention analysis of the 2008 financial crisis,
automatic outlier detection, and a rolling-origin forecast evaluation.

📄 **[Read the full analysis →](analysis/cement_consumption_analysis.md)**

![Series](analysis/figures/series-plot-1.png)

## Key results

| Model | Specification | RMSE | MAPE | 95% PI width | Passes all diagnostics |
|:--|:--|--:|--:|--:|:--:|
| A | SARIMA(2,1,1)(0,1,2)₁₂ — baseline | 87.9 | 7.0 % | 417 | ✗ |
| B | SARIMA(0,1,1)(0,1,2)₁₂ + calendar + crisis | 79.5 | 6.8 % | 376 | ✗ |
| C | Model B + outlier linearisation | 77.5 | 6.6 % | 317 | ✗ (1 borderline) |
| **D** | **ARIMA(0,1,1)(0,0,1)₁₂ + monthly dummies + calendar + crisis + outliers** | **72.5** | **6.2 %** | **308** | **✓** |

*Rolling-origin evaluation: models re-estimated at Dec 2013, 2014, 2015 and 2016, each forecasting the next 12 months (48 forecasts). Errors in thousand tonnes.*

* **Calendar effects matter a lot.** Each extra working day raises monthly
  consumption by ~1.1 %, and Easter falling in March instead of April lowers
  March consumption by ~10 %. Both effects are highly significant (|t| > 12).
* **The 2008 crisis is captured twice**: as jointly significant intervention
  dummies (LR test, p ≈ 0.007) and as automatically detected outliers,
  including a permanent **−16 % level shift in December 2008**.
* **Seasonality is deterministic.** The seasonal MA coefficients of the SARIMA
  models sum to −1, i.e. a unit root that cancels the seasonal difference
  (seasonal over-differencing). Replacing the seasonal difference with monthly
  dummies gives an invertible model that passes every residual test.
* The final model cuts RMSE by **17 %** and prediction-interval width by
  **26 %** versus the baseline while keeping ~98 % coverage.

![Backtest](analysis/figures/backtest-2017-1.png)

## Methods

1. **Exploration & stationarity** – log transform for variance, regular +
   seasonal differencing, over-differencing check via variances.
2. **Identification & estimation** – ACF/PACF, candidate SARIMA models
   compared by AIC/BIC and coefficient significance.
3. **Validation** – residual plots; Shapiro–Wilk, Anderson–Darling,
   Jarque–Bera, Breusch–Pagan, Durbin–Watson and Ljung–Box tests;
   stationarity/invertibility via characteristic roots.
4. **Extensions** – trading-day and Easter regressors, step interventions
   for crisis phases (likelihood-ratio test), iterative AO/LS/TC outlier
   detection, deterministic seasonality.
5. **Evaluation** – rolling-origin out-of-sample forecasts with RMSE, MAE,
   RMSPE, MAPE, interval coverage and width.

## Repository structure

```
├── analysis/
│   ├── cement_consumption_analysis.Rmd   # full analysis (source)
│   ├── cement_consumption_analysis.md    # rendered report (GitHub-readable)
│   └── figures/
├── R/
│   ├── calendar_effects.R    # trading-day & Easter regressors (Meeus algorithm)
│   ├── outlier_detection.R   # AO / LS / TC outlier detection and linearisation
│   ├── validation.R          # residual tests and diagnostic plots
│   └── forecasting.R         # forecasts, accuracy metrics, plots
├── data/cemento.dat          # monthly series, Jan 1988 – Dec 2017
├── report/original_course_report.pdf
└── install_packages.R
```

## Reproduce

```r
source("install_packages.R")
rmarkdown::render("analysis/cement_consumption_analysis.Rmd")
```

Requires R ≥ 4.1 (native pipe). Packages: `rmarkdown`, `knitr`, `nortest`,
`tseries`, `lmtest`.

## Data

Monthly apparent cement consumption in Spain (thousand tonnes), from the
*Boletín Estadístico* of the Spanish Ministry of Public Works
(Ministerio de Fomento). The analysis uses January 1990 – December 2017.

## About this project

Developed by **Maria Sans Bosch** and **David González** as the Time Series
project of the Data Science and Engineering degree at UPC (BarcelonaTech), 2025.
The outlier-detection and calendar-effect helpers are adapted from code
provided in the course.

The repository is a revised version of the course submission
([original report](report/original_course_report.pdf)). Changes include:

* Fixed data import: `read.csv()` treated the first observation as a header,
  shifting the whole series by one month and misaligning the calendar
  regressors.
* Fixed the crisis intervention windows (year boundaries were off by
  eleven months).
* Calendar regressors for the forecast horizon are now computed instead of
  set to zero.
* Added the deterministic-seasonality model (D) and replaced the single
  hold-out year with a rolling-origin evaluation.
* Refactored helper code into documented, reusable functions; the analysis
  now renders end to end.

Because of these fixes, the numbers differ from the original report.
