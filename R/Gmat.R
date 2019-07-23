Gmat <- function(n) {
  Gj <- function(j, n) {
    icomp <- complex(real = 0, imaginary = 1)
    w <- icomp * 2 * pi * (0:(n - 1) * j) / n
  }
  matAux <- matrix(NA, n, n)
  for (i in 0:(n - 1)) {
    matAux[(i + 1), ] <- exp(Gj(i, n))
  }
  return((1 / sqrt(n)) * matAux)
}
