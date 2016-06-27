context("Test estimate_fdr is working")
library(MetSizeR)

test_that("fdr outputs correctly", {

  n_mets <- 10
  n_resamps <- 10
  # t_stats for 10 metabolites resampled 10 times
  t_stat <- rnorm(n_mets * n_resamps) %>% matrix(n_resamps, n_mets)
  # standard deviations for t_stats
  sd <- rexp(n_mets * n_resamps) %>% matrix(n_resamps, n_mets)

  sig_metabs <- sample(1:n_mets, 1)

  res <- estimate_fdr(t_stat, sd, sig_metabs, p_stat = 0.95, TRUE)

  res %>%
    expect_is("numeric") %>%
    expect_length(1) %>%
    expect_gte(0)
})
