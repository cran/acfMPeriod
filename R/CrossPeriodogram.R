#' Cross-periodogram
#'
#' This function computes the cross-periodogram using harmonic regression.
#' @param series1 univariate time series
#' @param series2 univariate time series
#' @return a numeric vector containing the estimates of the cross-spectral density
#' @author Higor Cotta, Valdério A. Reisen, Pascal Bondon and Céline Lévy-Leduc
#' @references Fuller, Wayne A. Introduction to statistical time series.  John Wiley & Sons, 2009.
#' @export
CrossPeriodogram <- function(series1, series2) {
  # Build spectra for each series and combine them into cross-periodograms.
  n <- length(series1)
  period1 <- PerioRegAux(series1)
  period2 <- PerioRegAux(series2)

  .cross_period_from_components(period1, period2, n)
}


PerioRegAux <- function(series) {
  # Auxiliary spectrum used by cross-periodogram estimators.
  n <- length(series)
  .mirror_half_spectrum(.harmonic_fft(series, robust = FALSE), n)
}
