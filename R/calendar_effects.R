# =============================================================================
# Calendar-effect regressors for monthly series
#
# Adapted and refactored from helper functions provided in the Time Series
# course at UPC. The original hard-coded Easter table (1901-2050) has been
# replaced by the Meeus/Jones/Butcher algorithm, and the trading-day
# regressor no longer depends on the `chron` package.
# =============================================================================

#' Gregorian Easter Sunday for a vector of years (Meeus/Jones/Butcher).
#' @return Date vector
easter_date <- function(year) {
  a <- year %% 19
  b <- year %/% 100
  c <- year %% 100
  d <- b %/% 4
  e <- b %% 4
  f <- (b + 8) %/% 25
  g <- (b - f + 1) %/% 3
  h <- (19 * a + b - d - g + 15) %% 30
  i <- c %/% 4
  k <- c %% 4
  l <- (32 + 2 * e + 2 * i - h - k) %% 7
  m <- (a + 11 * h + 22 * l) %/% 451
  month <- (h + l - 7 * m + 114) %/% 31
  day <- ((h + l - 7 * m + 114) %% 31) + 1
  as.Date(sprintf("%04d-%02d-%02d", year, month, day))
}

#' Internal: year/month grid for a monthly ts starting at `start` with `n` obs.
.month_grid <- function(start, n) {
  idx <- (start[2] - 1) + 0:(n - 1)
  data.frame(year = start[1] + idx %/% 12, month = idx %% 12 + 1)
}

#' Trading-day regressor.
#'
#' For each month: (#weekdays) - (5/2) * (#weekend days).
#' It is 0 when weekdays and weekend days are in the "ideal" 5:2 ratio,
#' positive when the month has extra working days.
#'
#' @param start c(year, month) of the first observation
#' @param n     number of monthly observations
trading_days <- function(start, n) {
  grid <- .month_grid(start, n)
  values <- mapply(function(y, m) {
    first <- as.Date(sprintf("%04d-%02d-01", y, m))
    last  <- seq(first, by = "month", length.out = 2)[2] - 1
    days  <- seq(first, last, by = "day")
    wd    <- as.POSIXlt(days)$wday          # 0 = Sunday, 6 = Saturday
    weekend <- sum(wd %in% c(0, 6))
    (length(days) - weekend) - 5 / 2 * weekend
  }, grid$year, grid$month)
  ts(values, start = start, frequency = 12)
}

#' Easter regressor.
#'
#' Assumes Easter affects the `duration` days before Easter Sunday.
#' The regressor is +w in March and -w in April, with
#' w = (duration/2 - days_in_April) / duration, so it is 0 when the holiday
#' period is split evenly between March and April.
#'
#' @param start    c(year, month) of the first observation
#' @param n        number of monthly observations
#' @param duration length (days) of the Easter effect window
easter_effect <- function(start, n, duration = 6) {
  grid <- .month_grid(start, n)
  years <- unique(grid$year)
  ed <- easter_date(years)
  easter_month <- as.integer(format(ed, "%m"))
  easter_day   <- as.integer(format(ed, "%d"))
  days_in_april <- ifelse(easter_month == 4, pmin(easter_day, duration), 0)
  weight <- (duration / 2 - days_in_april) / duration
  w <- weight[match(grid$year, years)]
  values <- ifelse(grid$month == 3, w, ifelse(grid$month == 4, -w, 0))
  ts(values, start = start, frequency = 12)
}
