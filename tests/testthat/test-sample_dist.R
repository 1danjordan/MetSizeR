 context("sample_dist")

 test_that("sample_dist works", {
   sample_dist(iris[-5], 100, 10, 9, 3) %>%
     expect_is("list")
 })
