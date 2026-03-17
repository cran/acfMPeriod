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
  # Build robust spectra and combine them into robust cross-periodograms.
  n <- length(series1)
  period1 <- MPerioRegAux(series1)
  period2 <- MPerioRegAux(series2)

  .cross_period_from_components(period1, period2, n)
}

MPerioRegAux <- function(series) {
  # Auxiliary robust spectrum used by robust cross-periodogram estimators.
  n <- length(series)
  .mirror_half_spectrum(.harmonic_fft(series, robust = TRUE), n)
}
