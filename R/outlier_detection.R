# =============================================================================
# Automatic outlier detection for ARIMA models (Chen & Liu, 1993 style)
#
# Adapted from helper code provided in the Time Series course at UPC.
# Logic is unchanged; the code has been reformatted and documented.
#
# Detects three outlier types iteratively until no statistic exceeds `crit`:
#   AO - Additive Outlier   : one-off spike
#   LS - Level Shift        : permanent change of level
#   TC - Transitory Change  : shock decaying geometrically (delta = 0.7)
# =============================================================================

#' @param object fitted `arima` model (on the log series)
#' @param dif    differencing orders, e.g. c(1, 12) for (1-B)(1-B^12)
#' @param crit   critical value for |t|-like statistic (typically 2.8-3.5)
#' @param LS     whether to search for level shifts
#' @return list(atip = data.frame of outliers, sigma2, resid)
outdetec <- function(object, dif = c(0, 0), crit, LS = TRUE) {
  stopifnot(crit > 0)
  residuals <- object$residuals
  m <- length(residuals)

  # pi-weights of the full (differenced) model
  piweight <- -ARMAtoMA(ar = -object$model$theta, ma = -object$model$phi,
                        lag.max = m + sum(dif))
  if (dif[1] != 0) for (i in seq_len(dif[1])) piweight <- c(piweight, 0) - c(-1, piweight)
  if (length(dif) > 1) {
    for (i in 2:length(dif)) {
      if (dif[i] > 1) {
        piweight <- c(piweight, rep(0, dif[i])) - c(rep(0, dif[i] - 1), -1, piweight)
      }
    }
  }
  piweight <- piweight[1:m]

  delta <- 0.7
  c_ls <- cumsum(piweight) - 1                 # LS filter
  d_tc <- numeric(m)                           # TC filter
  d_tc[1] <- piweight[1] - delta
  for (i in 2:m) d_tc[i] <- delta * d_tc[i - 1] + piweight[i]

  sum1 <- 1 + sum(piweight^2)
  sum2 <- 1 + sum(c_ls^2)
  sum3 <- 1 + sum(d_tc^2)
  va <- mean(residuals^2)

  num <- type <- wcoeff <- LCrit <- NULL
  atip <- NULL
  maxL <- crit + 1

  while (maxL > crit) {
    ka1 <- sum1; ks1 <- sum2; kt1 <- sum3
    maxL <- 0
    for (i in 1:m) {
      suma1 <- sum(residuals[i:m] * c(1, -piweight[1:(m - i)]))
      suma2 <- sum(residuals[i:m] * c(1, -c_ls[1:(m - i)]))
      suma3 <- sum(residuals[i:m] * c(1, -d_tc[1:(m - i)]))

      ka1 <- ka1 - piweight[m - i + 1]^2
      w_ao <- suma1 / ka1; v_ao <- va / ka1; l_ao <- w_ao / sqrt(v_ao)

      ks1 <- ks1 - c_ls[m - i + 1]^2
      w_ls <- suma2 / ks1; v_ls <- va / ks1; l_ls <- w_ls / sqrt(v_ls)

      kt1 <- kt1 - d_tc[m - i + 1]^2
      w_tc <- suma3 / kt1; v_tc <- va / kt1; l_tc <- w_tc / sqrt(v_tc)

      if (abs(l_ao) > maxL && i > sum(dif)) {
        maxL <- abs(l_ao); t <- i; w <- w_ao; v <- v_ao; ts_type <- "AO"
      }
      if (abs(l_ls) > maxL && LS && i != m && i > sum(dif)) {
        maxL <- abs(l_ls); t <- i; w <- w_ls; v <- v_ls; ts_type <- "LS"
      }
      if (abs(l_tc) > maxL && i != m && i > sum(dif)) {
        maxL <- abs(l_tc); t <- i; w <- w_tc; v <- v_tc; ts_type <- "TC"
      }
    }

    if (maxL > crit) {
      filt <- switch(ts_type, AO = piweight, LS = c_ls, TC = d_tc)
      residuals[t:m] <- residuals[t:m] + w * c(-1, filt[1:(m - t)])
      val <- mean(residuals^2)
      l <- w / sqrt(v * val / va)
      va <- val

      num <- c(num, t); type <- c(type, ts_type)
      wcoeff <- c(wcoeff, w); LCrit <- c(LCrit, abs(l))
      atip <- data.frame(Obs = num, type_detected = type,
                         W_coeff = wcoeff, ABS_L_Ratio = LCrit)
    }
  }
  list(atip = atip, sigma2 = va, resid = residuals)
}

#' Remove detected outlier effects from a series (linearisation).
#' @param serie series on the same scale the model was fitted on (log)
#' @param atip  data.frame returned in outdetec()$atip
lineal <- function(serie, atip) {
  m <- length(serie)
  for (i in seq_len(nrow(atip))) {
    t <- atip$Obs[i]; type <- atip$type_detected[i]; w <- atip$W_coeff[i]
    if (type == "TC") serie[t:m] <- serie[t:m] - w * c(1, 0.7^(1:(m - t)))
    if (type == "LS") serie[t:m] <- serie[t:m] - w
    if (type == "AO") serie[t]   <- serie[t] - w
  }
  serie
}
