context("Test simulation functions work")

library(mvtnorm)
library(purrr)

data(UrineSpectra, package = "MetabolAnalyze")
pilot_data <- UrineSpectra[[1]]
covariates <- UrineSpectra[[2]][,2]

test_that("sim_PPCA_pilot outputs correct dimensions and type", {

  ppca_fit <- MetabolAnalyze::ppca.metabol(pilot_data, printout = FALSE)

  n <- 100
  sims <- sim_PPCA_pilot(n, ppca_fit, colMeans(pilot_data))

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), ncol(pilot_data))
  expect_is(sims, "matrix")
})

test_that("sim_PPCCA_pilot outputs correct dimensions and type", {

  ppcca_fit <- MetabolAnalyze::ppcca.metabol(pilot_data, covariates, printout = FALSE)

  n <- 100
  sims <- sim_PPCCA_pilot(n, ppcca_fit, colMeans(pilot_data))

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), ncol(pilot_data))
  expect_is(sims, "matrix")
})

test_that("sim_DPPCA_pilot outputs correct dimensions and type", {

  ppca_fit <- MetabolAnalyze::ppca.metabol(pilot_data, printout = FALSE)

  n <- 100
  sims <- sim_DPPCA_pilot(n, ppca_fit)

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), ncol(pilot_data))
  expect_is(sims, "matrix")
})

test_that("sim_PPCA outputs correct dimension and type", {

  n <- 100
  n_spectral_bins <- 40
  sims <- sim_PPCA(n, n_spectral_bins)

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), n_spectral_bins)
  expect_is(sims, "matrix")
})

test_that("sim_PPCCA outputs correct dimension and type", {
  n <- 100
  n_spectral_bins <- 40
  n_covars <- 3
  sims <- sim_PPCCA(n, n_spectral_bins, n_covars)

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), n_spectral_bins)
  expect_is(sims, "matrix")
})

test_that("sim_DPPCA outputs correct dimension and type", {
  n <- 100
  n_spectral_bins <- 40
  sims <- sim_DPPCA(n, n_spectral_bins)

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), n_spectral_bins)
  expect_is(sims, "matrix")
})
