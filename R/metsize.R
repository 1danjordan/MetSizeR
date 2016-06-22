#' Select model to generate simulations with
#'
#' @param model  a character vector of model to return
#' @param pilot  an indicator that pilot data exists
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

#'


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

# I genuinely feel like this is incorrect. It just isn't making
# sense
#'
#' est_sample_size <- function(results_sim, target.fdr) {
#'
#'   # index 1
#'   # this is the index of:
#'   #  1 take all the rows
#'   #  2 find the rows where the median FDR is less than 0.05
#'   #  3 get them minimum row number
#'   #  why get the minimum row number?...
#'
#'
#'   ind1 <- min(c(1:nrow(results_sim))[results_sim[, 2] < 0.05])
#'   # 1 ... ind1
#'
#'   # index 2
#'   # this is the index where:
#'   #   we take the result where the media FDR is greater than 0.05
#'   #   then we find the index of all of them
#'   #   the we get the highest index
#'
#'   ind2 <- max(c(1:ind1)[results_sim[1:ind1, 2] > 0.05])
#'   # 1 ... ind2
#'
#'   # is this supposed to be some form of linear interpolation?
#'   opty <- c(results_sim[ind1,2], results_sim[ind2,2])
#'   optx <- c(results_sim[ind1,1], results_sim[ind2,1])
#'
#'   # linear model of opty on optx
#'   optres <- lm(opty ~ optx)
#'
#'   # estimated sample size
#'   nhat <- round((target.fdr - optres$coef[1])/optres$coef[2])
#'
#'   # error message
#'   if(is.na(nhat)){print("Sorry: an error has occurred. Please rerun the function."); stop()}
#'
#'   # if sample sizes are the same
#'   if(n1 == n2)
#'   {
#'     # if nhat is even
#'     if(is.wholenumber(nhat/2)) {
#'       n1 <- n2 <- nhat/2
#'     } else {                # if nhat is not divisible by 2
#'         nhat <- nhat + 1    # add one to make it even
#'         n1 <- n2 <- nhat/2  # n1 and n2 are half that
#'         }
#'   } else {                  # if sample sizes are not the same
#'     if(is.wholenumber(n1 * nhat / (n1 + n2))) { # n1 is whole number
#'       n1 <- n1 * nhat / (n1 + n2)
#'       n2 <- nhat - n1
#'     } else {
#'       n1.user <- n1
#'       n2.user <- n2
#'       n1 <- ceiling(n1.user * nhat / (n1.user + n2.user))
#'       n2 <- ceiling(n2.user * nhat / (n1.user + n2.user))
#'       nhat <- n1 + n2
#'     }
#'   }
#'
#'   # return the estimated sample sizes
#'   est <- c(nhat, n1, n2)
#'   names(est) <- c("n","n1","n2")
#' }
#'
#' #' estimate the False Discovery Rate (FDR)
#' #'
#' #' @param test_stat
#'
#' estimate_fdr <- function(test_stat, ){
#'
#'   res.sampdist <- samp.dist(T,S,TS,x,y,n11,n22,in1n2,nn2,cpcf)
#'
#'   TS <- res.sampdist$TS
#'   S  <- res.sampdist$S
#'
#'   in1n2 <- 1/n1star + 1/n2star
#'   Add.sd <- sqrt(in1n2)
#'
#'   # calculating the shift in metabolites in grp 2
#'   vars <- S / Add.sd
#'   Add[,,k] <- delta / (vars * Add.sd)
#'
#'   # estimating the FDR
#'   tsB <- TS
#'
#'   # Add an increment factor to the t.statistic values
#'   # of the truly significant metabolites
#'   tsB[pos,] <- tsB[pos,] + Add[pos,,k]
#'   atsB <- abs(tsB)
#'
#'   # identifying a cut-off point (m-th largest absolute
#'   # value of the p TsB values)
#'   crit <- quantile(atsB, pstat)
#'
#'   # number of false positives
#'   errors <- (colSums(atsB > crit & !ind)) / (colSums(atsB > crit))
#'
#'   # median FDR of the Tstar permutations
#'   fdr_sim[k,s] <- quantile(errors[!is.na(errors)], 0.5)
#' }









