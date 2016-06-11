context("sample_dist")

library(purrr)
library(dplyr)

test_that("sample_dist works", {
  sample_dist(iris[-5], 100, 10, 9, 3) %>%
    expect_is("list")
})
