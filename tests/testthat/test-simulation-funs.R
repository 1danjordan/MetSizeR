context("Test simulation functions")

library(MetabolAnalyze)
data(UrineSpectra)

test_that("sim_pilot_PPCA works", {



  pilot_data <- UrineSpectra[[1]]
  ppca_fit <- ppca.metabol(pilot_data, minq = 2, maxq = 2)
  sims <- sim_DPPCA_pilot(100, ppca_fit)

  expect_
})
