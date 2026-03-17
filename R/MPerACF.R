#' Robust autocorrelation or autocovariance function estimation from the robust M-periodogram
#'
#' This function computer and plots(by default) the robust estimates of the autocovariance or the autocorrelation function for univariate and multivariate time series
#' based on the M-periodogram and the M-cross-periodogram.
#'
#' @param x a numeric vector or matrix.
#' @param lag.max maximum lag at which to calculate the acf. Default is 10*log10(N/m) where
#' N is the number of observations and m the number of series. Will be automatically limited
#' to one less than the number of observations in the series.
#' @param type character string giving the type of acf to be computed. Allowed values are "correlation" (the default) or "covariance".
#' Accepts parcial names.
#' @param plot logical. If TRUE (the default) the acf is plotted.
#' @param na.action function to be called to handle missing values. na.pass can be used.
#' @param demean logical. Should the covariances be about the sample means?
#' @param ... further arguments to be passed to plot.acf.
#' @return An object of class "robacf", which is a list with the following elements:
#' @return \code{lag} A three dimensional array containing the lags at which the acf is estimated.
#' @return \code{acf} An array with the same dimensions as lag containing the estimated acf.
#' @return \code{type} The type of correlation (same as the type argument).
#' @return \code{n.used} The number of observations in the time series.
#' @return \code{series} The name of the series x.
#' @return \code{snames} The series names for a multivariate time series.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Higor Cotta, Valderio Reisen, Pascal Bondon and Céline Lévy-Leduc. Part of the code re-used from the acf() function.
#' @references Fuller, Wayne A. Introduction to statistical time series.  John Wiley & Sons, 2009
#' @export
#' @import stats
#' @examples
#' data.set <- cbind(fdeaths, mdeaths)
#' MPerACF(data.set)
MPerACF <- function(x, lag.max = NULL, type = c("correlation", "covariance"), plot = TRUE, na.action = na.fail, demean = TRUE, ...) {
  # Validate inputs and normalize time-series structure.
  type <- match.arg(type)
  series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  x.freq <- frequency(x)
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (is.na(sampleT) || is.na(nser)) {
    stop("'sampleT' and 'nser' must be integer")
  }
  if (is.null(lag.max)) {
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  }
  lag.max <- as.integer(min(lag.max, sampleT - 1L))
  if (is.na(lag.max) || lag.max < 0) {
    stop("'lag.max' must be at least 0")
  }
  if (demean) {
    # Remove sample means before robust spectral estimation when requested.
    x <- sweep(x, 2, colMeans(x, na.rm = TRUE), check.margin = FALSE)
  }
  lag <- matrix(1, nser, nser)
  lag[lower.tri(lag)] <- -1
  acf.per <- array(1, c(lag.max, nser, nser))
  if (nser == 1L) { # Univariate
    # Diagonal (auto) term: inverse spectral transform of the M-periodogram.
    periodogram <- c(0, MPerioReg(x[, 1L]))
    acf.per[, 1, 1] <- .first_row_spectral_transform(periodogram, lag.max, scale = 2 * pi)
    if (type == "correlation" && lag.max > 0L) {
      acf.per <- acf.per / acf.per[1]
    }
  }
  else { # Multivariate
    # First compute all robust marginal auto-covariance sequences.
    for (i in seq_len(nser)) {
      periodogram <- c(0, MPerioReg(x[, i]))
      acf.per[, i, i] <- .first_row_spectral_transform(periodogram, lag.max, scale = 2 * pi)
    }
    # Then compute each robust cross pair once and fill both orientations.
    for (i in seq_len(nser - 1L)) {
      for (j in seq.int(i + 1L, nser)) {
        crossperi <- MCrossPeriodogram(x[, i], x[, j])
        acf.per[, i, j] <- .first_row_spectral_transform(crossperi$cross.periodxy, lag.max, scale = 1 / 2)
        acf.per[, j, i] <- .first_row_spectral_transform(crossperi$cross.periodyx, lag.max, scale = 1 / 2)
      }
    }
    if (type == "correlation" && lag.max > 0L) {
      # Convert robust covariance sequences to robust correlation sequences.
      var0 <- sqrt(diag(acf.per[1L, , , drop = TRUE]))
      for (i in seq_len(nser)) {
        acf.per[, i, i] <- acf.per[, i, i] / acf.per[1L, i, i]
      }
      for (i in seq_len(nser - 1L)) {
        for (j in seq.int(i + 1L, nser)) {
          scale.ij <- var0[i] * var0[j]
          acf.per[, i, j] <- acf.per[, i, j] / scale.ij
          acf.per[, j, i] <- acf.per[, j, i] / scale.ij
        }
      }
    }
  }
  lag <- outer(0:(lag.max - 1), lag / x.freq)
  acf.out <- structure(list(acf = acf.per, type = type, n.used = sampleT, lag = lag, series = series, snames = colnames(x)),
    class = "robacf"
  )
  if (plot) {
    plot.robacf(acf.out)
    invisible(acf.out)
  } else {
    acf.out
  }
}
