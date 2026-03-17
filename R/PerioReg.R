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
  # Estimate the harmonic-regression coefficients on positive frequencies.
  n <- length(series)
  fft.half <- .harmonic_fft(series, robust = FALSE)
  g <- length(fft.half)
  if (g == 0L) {
    return(numeric(0L))
  }

  scale.factor <- rep.int(sqrt(n / (8 * pi)), g)
  if ((n %% 2L) == 0L) {
    scale.factor[g] <- sqrt(n / (2 * pi))
  }

  # Convert coefficients into periodogram ordinates and mirror the spectrum.
  perior.half <- Mod(scale.factor * fft.half)^2
  .mirror_half_spectrum(perior.half, n)
}
