context("standardize")

test_that("standardize outputs correctly", {
  x <- rnorm(10)

  expect_equal(range(standardize(x)), c(0,1))
  expect_equal(standardize(1:5), c(0, 0.25, 0.5, 0.75, 1))

})

test_that("standardize is the same as the deprecated function", {
  # Remove Species variable
  iris_data <- iris[, -5]

  # Define the old function
  standardize2 <- function(C) {
    for(i in 1:ncol(C)) {
      rg <- range(C[,i])
      C[,i] <- (C[,i]-min(C[,i]))/(rg[2] - rg[1])
    }
    C
  }

  expect_equal(
    purrr::dmap(iris_data, standardize),
    standardize2(iris_data)
  )
})
