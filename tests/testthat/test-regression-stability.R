test_that("periodogram-based outputs stay numerically stable", {
  expected_perio <- c(
    54355.5004592646,
    10451.5142545753,
    5911.49039673241,
    7837.75364680223,
    11483.3082372142,
    1690957.8127608
  )

  expected_mperio <- c(
    53096.2442301288,
    10755.4666474105,
    5572.1002296534,
    7764.77153767615,
    10863.7903632865,
    1690957.8127608
  )

  expected_cross <- c(
    -7140.71965174661 + 0i,
    98342.7771705888 + 27062.4297491946i,
    9634.99830810556 - 14936.5146694394i,
    13960.1470963724 + 7397.78764404646i,
    17709.809723643 + 1461.74403390652i,
    28204.5900689397 - 17511.6768138825i
  )

  expected_mcross <- c(
    -12046.2056265516 + 0i,
    87745.378287788 + 33464.1227416662i,
    10484.9192376837 - 9508.87928105146i,
    12928.8542613639 + 3954.9869757276i,
    12791.1633558138 + 2420.96648872139i,
    26127.2685076561 - 12775.5721609452i
  )

  expect_equal(head(PerioReg(as.numeric(ldeaths)), 6), expected_perio, tolerance = 1e-8)
  expect_equal(head(MPerioReg(as.numeric(ldeaths)), 6), expected_mperio, tolerance = 1e-8)

  cross <- CrossPeriodogram(as.numeric(fdeaths), as.numeric(mdeaths))
  mcross <- MCrossPeriodogram(as.numeric(fdeaths), as.numeric(mdeaths))

  expect_equal(head(cross$cross.periodxy, 6), expected_cross, tolerance = 1e-8)
  expect_equal(head(mcross$cross.periodxy, 6), expected_mcross, tolerance = 1e-8)
})


test_that("acf and wrapper outputs remain stable", {
  expected_per_acf <- c(
    1,
    0.749803706229933,
    0.384087582603209,
    -0.0101674551483199,
    -0.402726876390332,
    -0.666980684626753,
    -0.742410586859306,
    -0.663511779086587
  )

  expected_mper_acf <- c(
    1,
    0.755188619733236,
    0.395335523918558,
    0.005546433731047,
    -0.411109270778229,
    -0.679900418902508,
    -0.75801773585186,
    -0.679018117915031
  )

  per_out <- PerACF(as.numeric(ldeaths), lag.max = 8, type = "correlation", plot = FALSE)
  mper_out <- MPerACF(as.numeric(ldeaths), lag.max = 8, type = "correlation", plot = FALSE)

  expect_equal(as.numeric(per_out$acf[, 1, 1]), expected_per_acf, tolerance = 1e-8)
  expect_equal(as.numeric(mper_out$acf[, 1, 1]), expected_mper_acf, tolerance = 1e-8)

  covcor_per_expected <- matrix(c(1, 0.975698292635721, 0.975698292635721, 1), nrow = 2)
  covcor_mper_expected <- matrix(c(1, 0.972436074347327, 0.972436074347327, 1), nrow = 2)

  covcor_per <- CovCorPer(cbind(fdeaths, mdeaths), type = "correlation")
  covcor_mper <- CovCorMPer(cbind(fdeaths, mdeaths), type = "correlation")

  expect_equal(covcor_per, covcor_per_expected, tolerance = 1e-8)
  expect_equal(covcor_mper, covcor_mper_expected, tolerance = 1e-8)
})


test_that("multivariate acf structures are symmetric", {
  x <- cbind(fdeaths, mdeaths)

  per_cov <- PerACF(x, lag.max = 6, type = "covariance", plot = FALSE)
  mper_cov <- MPerACF(x, lag.max = 6, type = "covariance", plot = FALSE)

  expect_equal(dim(per_cov$acf), c(6, 2, 2))
  expect_equal(dim(mper_cov$acf), c(6, 2, 2))

  expect_equal(per_cov$acf[, 1, 2], per_cov$acf[, 2, 1], tolerance = 1e-8)
  expect_equal(mper_cov$acf[, 1, 2], mper_cov$acf[, 2, 1], tolerance = 1e-8)

  expect_s3_class(per_cov, "acf")
  expect_s3_class(mper_cov, "robacf")
})
