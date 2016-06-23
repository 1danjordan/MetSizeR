context("Test getting model functions")

library(purrr)

test_that("get_model_fun returns functions", {

  models <- list("ppca", "ppcca", "dppca")

  # Models for without pilot data
  models %>%
    map(get_model_fun %>% expect_type("closure"))

  # Models for with pilot data
  models %>%
    map(get_model_fun(pilot = TRUE) %>% expect_type("closure"))
})
