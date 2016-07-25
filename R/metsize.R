#' Select model to generate simulations with
#'
#' @param model  a character vector of model to return
#' @param pilot  an indicator that pilot data exists
#'
#' @return a closure
#'
#' @example
#' # get PPCA simulator when pilot data is not present
#' get_model_fun(model = "ppca")
#'
#' # get DPPCA simulator when pilot data is present
#' get_model_fun(model = "dppca", pilot = TRUE)

get_model_fun <- function(model, pilot = NULL) {

  if(is_null(pilot)) {
    switch(model,
           ppca  = sim_PPCA,
           ppcca = sim_PPCCA,
           dppca = sim_DPPCA)
  } else {
    switch(model,
           ppca  = sim_PPCA_pilot,
           ppcca = sim_PPCCA_pilot,
           dppca = sim_DPPCA_pilot)
  }
}

#' Estimate the False Discovery Rate (FDR)
#'
#' This function takes the output of sample_dist, the resampled
#' test statistics of simulated metabolomic data, and estimates
#' the median false discovery rate according to some level of
#' confidence \code{p_stat}.
#'
#' @param test_stat     a matrix (n x m)
#' @param sd            a matrix (n x m)
#' @param pilot         boolean indicating if pilot data is present
#' @param sig_metabs    a vector indicating significant metabolites
#' @param p_stat        level of confidence (0 < x < 1)
#'
#' @return              the median false discovery rate (numeric)

estimate_fdr <- function(t_stat, sd, sig_metabs, p_stat, pilot){

  delta <- ifelse(pilot, qnorm(0.99), qnorm(0.89))
  ind <- matrix(FALSE, nrow = nrow(t_stat), ncol = ncol(t_stat))
  ind[, sig_metabs] <- TRUE

  t_stat[sig_metabs] <- t_stat[sig_metabs] + (delta/sd)[sig_metabs]

  abs_t_stat <- abs(t_stat)
  crit <- quantile(abs_t_stat, p_stat)

  # number of false positives
  errors <- (colSums(abs_t_stat > crit & !ind)) / (colSums(abs_t_stat > crit))

  # median FDR
  quantile(errors[!is.na(errors)], 0.5)
}

#' Compute the number of samples required for each treatment group
#'
#' This function takes an estimated sample size and the number
#' of samples in each treatment from a pilot experiment and computes
#' the required treatment size for an experiment.
#'
#' @param nhat  the sample size required for an experiment
#' @param n1    the size of the sample in group 1 of the pilot
#' @param n2    the size of the sample in group 2 of the pilot
#'
#' @return      a named list


compute_treatment_sizes <- function(nhat, n1, n2) {

  is_whole <- function(x) { x%%1 == 0 }

  if (is.na(nhat)) {
    print("Sorry: an error has occurred. Please rerun the function.")
    stop()
  }

  # if sample sizes are the same
  if(n1 == n2) {
    if(is_whole(nhat / 2)) {          # if nhat is even
      treatment_1 <- nhat / 2
      treatment_2 <- nhat / 2
    } else {                          # if nhat is uneven
      treatment_1 <- (nhat + 1) / 2
      treatment_2 <- (nhat + 1) / 2
    }
  } else {                            # if sample sizes are not the same
    if(is_whole(n1 * nhat / (n1 + n2))) {
      treatment_1 <- n1 * nhat / (n1 + n2)
      treatment_2 <- nhat - n1
    } else {
      treatment_1 <- ceiling(n1 * nhat / (n1 + n2))
      treatment_2 <- ceiling(n2 * nhat / (n1 + n2))

    }
  }

  return(list(
    samples = treatment_1 + treatment_2,
    treatment_1 = treatment_1,
    treatment_2 = treatment_2))
}

#' Determine the sample size at which the FDR line is equal to 0.05
#'
#' @param df          a data frame of sample sizes and False Discovery Rates
#' @param target_fdr  the target false discovery rate
#' @param n1          the number of samples in treatment 1 of pilot
#' @param n2          the number of samples in treatment 2 of pilot
#'
#' @result            a named list

est_sample_size <- function(df, target_fdr, n1, n2) {

  upp_bound <- df %>%
    filter(fdr < 0.05) %>%
    slice(which.min(n_samps))

  low_bound <- df %>%
    filter(fdr > 0.05) %>%
    slice(which.max(n_samps))

  lin_approx <- lm(fdr ~ n_samps, data = bind_rows(upp_bound, low_bound))
  nhat <- round((target_fdr - lin_approx$coef[1]) / lin_approx$coef[2])

  # Remove "(Intercept)" attr
  nhat <- unname(nhat)
  return(compute_treatment_sizes(nhat, n1, n2))
}



