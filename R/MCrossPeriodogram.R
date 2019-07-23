#' Robust M-cross-periodogram
#'
#' This function computes the Robust M-cross-periodogram using M-regression.
#' @param series1 univariate time series
#' @param series2 univariate time series
#' @return a numeric vector containing the estimates of the cross-spectral density
#' @author Higor Cotta, Valdério A. Reisen, Pascal Bondon and Céline Lévy-Leduc
#' @references Fuller, Wayne A. Introduction to statistical time series.  John Wiley & Sons, 2009.
#' @export
MCrossPeriodogram <- function(series1, series2) {
  n <- length(series1)
  period1 <- MPerioRegAux(series1)
  period2 <- MPerioRegAux(series2)
  cross.periodxy <- NULL
  cross.periodyx <- NULL
  for (i in 1:(n - 1)) {
    alpha.period1 <- Re(period1[i])
    beta.period1 <- Im(period1[i])
    alpha.period2 <- Re(period2[i])
    beta.period2 <- Im(period2[i])
    cross.periodxy[i] <- (n / 2) * complex(
      real = (alpha.period1 * alpha.period2 + beta.period1 * beta.period2),
      imaginary = -(alpha.period1 * beta.period2 - alpha.period2 * beta.period1)
    )
    cross.periodyx[i] <- (n / 2) * complex(
      real = (alpha.period1 * alpha.period2 + beta.period1 * beta.period2),
      imaginary = -(alpha.period2 * beta.period1 - alpha.period1 * beta.period2)
    )
  }
  cross.periodxy <- c((2 * n * alpha.period1 * alpha.period2), cross.periodxy)
  cross.periodyx <- c((2 * n * alpha.period1 * alpha.period2), cross.periodyx)
  cross.period <- list(cross.periodxy = cross.periodxy, cross.periodyx = cross.periodyx)
  return(cross.period)
}

MPerioRegAux <- function(series) {
  n <- length(series)
  perior <- FFT <- NULL
  g <- n %/% 2
  for (j in 1:g) {
    X1 <- X2 <- NULL
    w <- 2 * pi * j / n
    for (i in 1:n) {
      X1[i] <- cos(w * i)
      X2[i] <- sin(w * i)
    }
    if (j != (n / 2)) {
      MX <- cbind(X1, X2)
      fitrob <- rlm(series ~ MX - 1, method = "M", psi = MASS::psi.huber)
      FFT[j] <- complex(real = fitrob$coef[1], imaginary = -fitrob$coef[2])
    }
    else {
      MX <- cbind(X1)
      fitrob <- rlm(series ~ MX - 1, method = "M", psi = MASS::psi.huber)
      FFT[j] <- complex(real = fitrob$coef[1], imaginary = -0)
    }
    perior[j] <- (FFT[j])
  }
  if ((n %% 2) != 0) {
    return(c(perior, rev(perior)))
  } else {
    return(c(perior, rev(perior))[-g])
  }
}
