.fit_harmonic_model <- function(series, design, robust = FALSE) {
  # Fit the harmonic regression either classically or with robust M-estimation.
  if (robust) {
    fit <- MASS::rlm(
      x = design,
      y = series,
      method = "M",
      psi = MASS::psi.huber
    )
  } else {
    fit <- stats::lm.fit(x = design, y = series)
  }

  as.numeric(fit$coefficients)
}

.harmonic_fft <- function(series, robust = FALSE) {
  # Estimate Fourier-like harmonic coefficients by frequency-wise regression.
  series <- as.numeric(series)
  n <- length(series)
  g <- n %/% 2L
  if (g == 0L) {
    return(complex(0L))
  }

  idx <- seq_len(n)
  nyquist_j <- n / 2
  fft_half <- complex(length.out = g)

  for (j in seq_len(g)) {
    # Build sine/cosine regressors for frequency j.
    w <- 2 * pi * j / n
    x1 <- cos(w * idx)

    if (j != nyquist_j) {
      design <- cbind(x1, sin(w * idx))
      coef <- .fit_harmonic_model(series, design, robust = robust)
      fft_half[j] <- complex(real = coef[1L], imaginary = -coef[2L])
    } else {
      coef <- .fit_harmonic_model(series, matrix(x1, ncol = 1L), robust = robust)
      fft_half[j] <- complex(real = coef[1L], imaginary = 0)
    }
  }

  fft_half
}

.mirror_half_spectrum <- function(values, n) {
  # Reconstruct the full spectrum from the positive-frequency half.
  g <- length(values)
  if (g == 0L) {
    return(values)
  }

  mirrored <- c(values, rev(values))
  if ((n %% 2L) != 0L) {
    mirrored
  } else {
    mirrored[-g]
  }
}

.cross_period_from_components <- function(period1, period2, n) {
  # Combine two complex spectra into forward/backward cross-periodograms.
  if (length(period1) == 0L || length(period2) == 0L) {
    return(list(cross.periodxy = complex(0L), cross.periodyx = complex(0L)))
  }

  cross.periodxy <- (n / 2) * period1 * Conj(period2)
  cross.periodyx <- (n / 2) * Conj(period1) * period2
  lead.term <- 2 * n * Re(period1[length(period1)]) * Re(period2[length(period2)])

  list(
    cross.periodxy = c(lead.term, cross.periodxy),
    cross.periodyx = c(lead.term, cross.periodyx)
  )
}

.first_row_spectral_transform <- function(spectrum, lag.max, scale = 1) {
  # Inverse FFT recovers the first covariance row used as ACF/ACOVF lags.
  if (lag.max == 0L) {
    return(numeric(0L))
  }

  transformed <- scale * Re(fft(spectrum, inverse = TRUE) / length(spectrum))
  transformed[seq_len(lag.max)]
}
