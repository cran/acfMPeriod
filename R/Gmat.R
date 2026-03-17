Gmat <- function(n) {
  # Build the unitary Fourier matrix used in the original formulation.
  idx <- 0:(n - 1)
  phase <- outer(idx, idx, FUN = "*")
  (1 / sqrt(n)) * exp((0 + 1i) * 2 * pi * phase / n)
}
