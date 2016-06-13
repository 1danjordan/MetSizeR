context("Test simulation functions work")

library(MetabolAnalyze, quietly = TRUE)
data(UrineSpectra)
pilot_data <- UrineSpectra[[1]]
covariates <- UrineSpectra[[2]][,2]
n <- 100

# I want to pull all of these out into a function because
# they are all the same tests but I'm not sure how
# to do it properly



test_that("sim_PPCA_pilot outputs correct dimensions and type", {

  ppca_fit <- ppca.metabol(pilot_data)

  sims <- sim_PPCA_pilot(n, ppca_fit, colMeans(pilot_data))

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), ncol(pilot_data))
  expect_is(sims, "matrix")
})

test_that("sim_PPCCA_pilot outputs correct dimensions and type", {

  ppcca_fit <- ppcca.metabol(pilot_data, covariates)

  sims <- sim_PPCCA_pilot(n, ppcca_fit, colMeans(pilot_data))

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), ncol(pilot_data))
  expect_is(sims, "matrix")
})

test_that("sim_DPPCA_pilot outputs correct dimensions and type", {

  ppca_fit <- ppca.metabol(pilot_data)

  sims <- sim_DPPCA_pilot(n, ppca_fit)

  expect_equal(nrow(sims), n)
  expect_equal(ncol(sims), ncol(pilot_data))
  expect_is(sims, "matrix")
})
