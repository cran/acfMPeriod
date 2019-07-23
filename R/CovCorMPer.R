#' Robust covariance or correlation matrix from the MPer-ACF
#'
#' Wrapper that computes the covariance or correlation matrix of \code{x} at lag 0 obtained from the robust MPer-ACF.
#' @param x a numeric matrix
#' @param type character string giving the type of acf to be computed. Allowed values are "correlation" (the default) or "covariance".
#' @return a numeric matrix
#' @export
#' @examples
#' data.set <- cbind(fdeaths, mdeaths)
#' CovCorMPer(data.set)
CovCorMPer <- function(x, type = c("correlation", "covariance")) {
  type <- match.arg(type)
  mat.comp <- MPerACF(x, lag.max = 1, plot = FALSE, type)$acf
  return(mat.comp[1, , ])
}
