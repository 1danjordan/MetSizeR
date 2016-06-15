#' Determine the sample size at which the FDR line is equal to 0.05
#'
#'
#' @param sims    simulation of metabolites
#' @param props

# results_sim is a 4D matrix with dimensions:
#  1 - number of increments
#  2 - median FDR of resampled simulations
#  3 - 90th percentile FDR of resampled simulations
#  4 - 10th percentile FDR of resampled simulations

est_sample_size <- function(results_sim, target.fdr) {

  # index 1
  # this is the index of:
  #  1 take all the rows
  #  2 find the rows where the median FDR is less than 0.05
  #  3 get them minimum row number
  # I feel like this should be sorted first...

  ind1 <- min(c(1:nrow(results_sim))[results_sim[, 2] < 0.05])
  # index 2
  ind2 <- max(c(1:ind1)[results_sim[1:ind1,2]>0.05])

  opty <- c(results_sim[ind1,2], results_sim[ind2,2])
  optx <- c(results_sim[ind1,1], results_sim[ind2,1])

  # linear model of opty on optx
  optres <- lm(opty~optx)

  # estimated sample size
  nhat <- round((target.fdr - optres$coef[1])/optres$coef[2])

  # error message
  if(is.na(nhat)){print("Sorry: an error has occurred. Please rerun the function."); stop()}

  # if sample sizes are the same
  if(n1 == n2)
  {
    # if nhat is even
    if(is.wholenumber(nhat/2)) {
      n1 <- n2 <- nhat/2
    } else {                # if nhat is not divisible by 2
        nhat <- nhat + 1    # add one to make it even
        n1 <- n2 <- nhat/2  # n1 and n2 are half that
        }
  } else {                  # if sample sizes are not the same
    if(is.wholenumber(n1 * nhat / (n1 + n2))) { # n1 is whole number
      n1 <- n1 * nhat / (n1 + n2)
      n2 <- nhat - n1
    } else {
      n1.user <- n1
      n2.user <- n2
      n1 <- ceiling(n1.user * nhat / (n1.user + n2.user))
      n2 <- ceiling(n2.user * nhat / (n1.user + n2.user))
      nhat <- n1 + n2
    }
  }

  # return the estimated sample sizes
  est <- c(nhat, n1, n2)
  names(est) <- c("n","n1","n2")
}

#' estimate the False Discovery Rate (FDR)
#'
#' @param test_stat

estimate_fdr <- function(test_stat, ){

}