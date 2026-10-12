# =============================================================================
# Forecasting utilities for ARIMA models fitted on log-transformed series
# =============================================================================

#' Point forecasts and 95% prediction intervals back on the original scale.
#'
#' @param model   fitted `arima` object (fitted on the log series)
#' @param h       forecast horizon
#' @param newxreg regressors for the forecast horizon (NULL if none)
#' @param wLS     total log-level shift removed by outlier linearisation;
#'                added back so forecasts are on the level of the observed data
#' @return data.frame with time, lower, forecast, upper (original units)
arima_forecast <- function(model, h = 12, newxreg = NULL, wLS = 0) {
  pr <- predict(model, n.ahead = h, newxreg = newxreg)
  data.frame(
    time     = as.numeric(time(pr$pred)),
    lower    = as.numeric(exp(pr$pred + wLS - 1.96 * pr$se)),
    forecast = as.numeric(exp(pr$pred + wLS)),
    upper    = as.numeric(exp(pr$pred + wLS + 1.96 * pr$se))
  )
}

#' Accuracy of a forecast against the realised values.
forecast_metrics <- function(fc, actual) {
  actual <- as.numeric(actual)
  err <- actual - fc$forecast
  c(RMSE     = sqrt(mean(err^2)),
    MAE      = mean(abs(err)),
    RMSPE    = sqrt(mean((err / actual)^2)),
    MAPE     = mean(abs(err / actual)),
    Coverage = mean(actual >= fc$lower & actual <= fc$upper),
    CI_width = mean(fc$upper - fc$lower))
}

#' Plot observed series with forecast and prediction interval.
plot_forecast <- function(series, fc, from = 2013, main = "") {
  freq <- frequency(series)
  shown <- window(series, start = from)
  lo <- ts(fc$lower, start = fc$time[1], frequency = freq)
  pt <- ts(fc$forecast, start = fc$time[1], frequency = freq)
  up <- ts(fc$upper, start = fc$time[1], frequency = freq)
  rng <- range(shown, fc$lower, fc$upper)
  yl <- c(rng[1] * 0.95, rng[2] * 1.15)            # headroom for the legend
  # ts.plot() only forwards a fixed set of graphical args; the rest go in gpars
  ts.plot(shown, lo, up, pt, lty = c(1, 2, 2, 1), col = c(1, 4, 4, 2),
          type = "o", ylim = yl, xlim = c(from, max(fc$time) + 1 / freq),
          main = main, ylab = "Thousand tonnes", xlab = "",
          gpars = list(pch = 20, cex = 0.5))
  abline(v = from:ceiling(max(fc$time) + 1), lty = 3, col = "grey60")
  legend("top", c("Observed", "Forecast", "95% PI"), col = c(1, 2, 4),
         lty = c(1, 1, 2), horiz = TRUE, bty = "n", cex = 0.8)
}
