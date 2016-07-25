context("Test that we can estimate required sample sizes")

library(MetSizeR)
library(purrr, quietly = TRUE)

test_that("compute_treatment_sizes follows correct logic", {

  inputs <- list(
    equal_and_even      = list(nhat = 100, n1 = 20, n2 = 20),
    unequal_and_even    = list(nhat = 100, n1 = 10, n2 = 20),
    equal_and_uneven    = list(nhat = 101, n1 = 20, n2 = 20),
    unequal_and_uneven  = list(nhat = 101, n1 = 10, n2 = 20)
  )

  outputs <- list(
    equal_and_even     = list(nhat = 100, treatment_1 = 50, treatment_2 = 50),
    unequal_and_even   = list(nhat = 101, treatment_1 = 34, treatment_2 = 67),
    equal_and_uneven   = list(nhat = 102, treatment_1 = 51, treatment_2 = 51),
    unequal_and_uneven = list(nhat = 102, treatment_1 = 34, treatment_2 = 68)
  )

  lifted_fun <- lift_dl(compute_treatment_sizes)
  results <- map(inputs, lifted_fun)

  list(results, outputs) %>% pmap(~ expect_equal(.x, .y, check.names = FALSE))
})


test_that("test est_sample_size runs correctly", {

  df <- data.frame(n_samps = c(10, 20), fdr = c(1, 0))

  result <- est_sample_size(df, 0, 10, 10)
  output <- list(samples = 20, treatment_1 = 10, treatment_2 = 10)
  expect_that(result, output)
})

