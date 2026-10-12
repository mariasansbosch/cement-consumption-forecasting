# =============================================================================
# Residual diagnostics for (seasonal) ARIMA models
# =============================================================================

#' t-ratios of the estimated coefficients
t_ratios <- function(model) {
  round(model$coef / sqrt(diag(model$var.coef)), 2)
}

#' Diagnostic plots for the residuals of an ARIMA model.
#' @param model fitted `arima` object
#' @param s     seasonal period
plot_residual_diagnostics <- function(model, s = 12) {
  res <- residuals(model)
  op <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
  on.exit(par(op))

  plot(res, main = "Residuals", ylab = "")
  abline(h = 0)
  abline(h = c(-3, 3) * sd(res), lty = 3, col = 4)

  scatter.smooth(sqrt(abs(res)), main = "sqrt(|residuals|)  (homoscedasticity)",
                 lpars = list(col = 2), xlab = "", ylab = "")

  qqnorm(res, main = "Normal Q-Q")
  qqline(res, col = 2, lwd = 2)

  hist(res, breaks = 20, freq = FALSE, main = "Histogram", xlab = "")
  curve(dnorm(x, mean = mean(res), sd = sd(res)), col = 2, add = TRUE)

  par(mfrow = c(1, 2))
  acf(res, ylim = c(-1, 1), lag.max = 60, col = c(2, rep(1, s - 1)), main = "ACF residuals")
  pacf(res, ylim = c(-1, 1), lag.max = 60, col = c(rep(1, s - 1), 2), main = "PACF residuals")
}

#' Formal residual tests, returned as a tidy data.frame.
#' @param model fitted `arima` object
#' @param y     the series the model was fitted on (needed for Breusch-Pagan)
#' @param s     seasonal period
residual_tests <- function(model, y, s = 12) {
  res <- as.numeric(residuals(model))
  fitted_vals <- as.numeric(y) - res

  lb <- sapply(c(s, 2 * s, 4 * s), function(l)
    Box.test(res, type = "Ljung-Box", lag = l, fitdf = sum(model$arma[1:4]))$p.value)

  data.frame(
    check = c("Normality", "Normality", "Normality",
              "Homoscedasticity", "Independence",
              rep("Independence", 3)),
    test = c("Shapiro-Wilk", "Anderson-Darling", "Jarque-Bera",
             "Breusch-Pagan", "Durbin-Watson",
             paste0("Ljung-Box (lag ", c(s, 2 * s, 4 * s), ")")),
    p_value = c(
      shapiro.test(res)$p.value,
      nortest::ad.test(res)$p.value,
      tseries::jarque.bera.test(res)$p.value,
      lmtest::bptest(res ~ fitted_vals)$p.value,
      lmtest::dwtest(res ~ seq_along(res))$p.value,
      lb
    )
  ) |>
    transform(p_value = signif(p_value, 3),
              passes_5pct = ifelse(p_value > 0.05, "yes", "no"))
}

#' Moduli of the AR and MA characteristic-polynomial roots
#' (all > 1  <=>  stationary and invertible).
root_moduli <- function(model) {
  list(
    AR = if (length(model$model$phi))   Mod(polyroot(c(1, -model$model$phi)))   else NA,
    MA = if (length(model$model$theta)) Mod(polyroot(c(1,  model$model$theta))) else NA
  )
}

#' One-line check of stationarity / invertibility.
#' Roots with modulus within `tol` of 1 are flagged as lying on the unit circle.
check_roots <- function(model, tol = 0.01) {
  r <- root_moduli(model)
  verdict <- function(x, ok, bad) {
    if (all(is.na(x))) return("no roots")
    if (min(x) < 1) return(bad)
    if (min(x) < 1 + tol) return(paste(bad, "(root on the unit circle)"))
    ok
  }
  cat(sprintf("Min |AR root| = %.4f  -> %s\n", min(r$AR), verdict(r$AR, "stationary", "NOT stationary")))
  cat(sprintf("Min |MA root| = %.4f  -> %s\n", min(r$MA), verdict(r$MA, "invertible", "NOT invertible")))
}
