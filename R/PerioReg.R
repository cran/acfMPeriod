#' Periodogram
#'
#' This function computes the univariate periodogram using harmonic regression.
#' @param series univariate time series
#' @return a numeric vector containing the robust estimates of the spectral density
#' @author Higor Cotta, Valdério A. Reisen, Pascal Bondon and Céline Lévy-Leduc.
#' @references Reisen, V. A. and Lévy-Leduc, C. and Taqqu, M. (2017) An M-estimator for the long-memory parameter. \emph{Journal of Statistical Planning and Inference},  187, 44-55.
#' @references Fuller, Wayne A. Introduction to statistical time series.  John Wiley & Sons, 2009.
#' @export
#' @examples
#' PerioReg(ldeaths)
PerioReg <- function(series) {
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
      fit <- lm(series ~ MX - 1)
      FFT[j] <- sqrt(n / (8 * pi)) * complex(real = fit$coef[1], imaginary = -fit$coef[2])
    }
    else {
      MX <- cbind(X1)
      fit <- lm(series ~ MX - 1)
      FFT[j] <- sqrt(n / (2 * pi)) * complex(real = fit$coef[1], imaginary = -0)
    }
    perior[j] <- Mod(FFT[j])^2
  }
  if ((n %% 2) != 0) {
    return(c(perior, rev(perior)))
  } else {
    return(c(perior, rev(perior))[-g])
  }
}
