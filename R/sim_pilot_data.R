#' Simulate data using a PPCCA model when pilot data is present
#'
#'@param n                  number of simulated samples
#'@param ppca_obj           the results from \code{ppca.metabol}
#'@param pilot_var_means    the means of the pilot data variables
#'
#' @examples
#' # Simulate 100 samples using the PPCA model with pilot data
#' # using the Urine Spectra dataset
#'
#' library(MetabolAnalyze)
#' data(UrineSpectra)
#'
#' pilot_data <- UrineSpectra[[1]]
#' ppca_fit <- ppca.metabol(pilot_data, minq = 2, maxq = 2)
#'
#' sim_pilot_PPCA(100, ppca_fit, colMeans(pilot_data))

sim_pilot_PPCA <- function(n, ppca_obj, pilot_var_means) {

  n_pilot_vars  <- nrow(ppca_obj$loadings)
  n_latent_dims <- ncol(ppca_obj$loadings)

  u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims))

  x <- rmvnorm(n, rep(0, n_pilot_vars), diag(ppca_obj$sig, n_pilot_vars)) +
         tcrossprod(u, ppca_obj$loadings) +
         matrix(pilot_var_means, n, n_pilot_vars, byrow=TRUE)

  return(x)
}

#' Simulate data using PPCA model when pilot and covariate data present
#'
#' @param n                   number of simulated samples
#' @param ppcca_obj           the results from \code{ppcca.metabol}
#' @param pilot_var_means     the means of the pilot data variables
#'
#' @examples
#' # Simulate 100 samples using the PPCCA model with pilot data
#' # and covariates using the Urine Spectra dataset
#'
#' library(MetabolAnalyze)
#' data(UrineSpectra)
#'
#' pilot_data <- UrineSpectra[[1]]
#' covariates <- UrineSpectra[[2]][,2]
#' ppcca_fit <- ppcca.metabol(pilot_data, covariates, minq = 2, maxq = 2)
#'
#' sim_PPCCA_pilot(100, ppcca_fit, colMeans(pilot_data))

# Sampling 1 obs produces NaNs
# returns matrix with dimnames (as a result of $loadings)

sim_PPCCA_pilot <- function(n, ppcca_obj, pilot_var_means) {

  n_pilot_vars  <- nrow(ppcca_obj$loadings)
  n_latent_dims <- ncol(ppcca_obj$loadings)
  n_covars      <- ncol(ppcca_obj$coefficients) - 1

  c <- rmvnorm(n, rep(0, n_covars), diag(n_covars)) %>%
    apply(2, standardize) %>%
    cbind(1, .)

  u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims)) +
    t(tcrossprod(ppcca_obj$coefficients, c))

  x <- rmvnorm(n, rep(0, n_pilot_vars), diag(ppcca_obj$sig, n_pilot_vars)) +
    tcrossprod(u, ppcca_obj$loadings) +
    matrix(pilot_var_means, n, n_pilot_vars, byrow = TRUE)

  return(x)

 }


#' Simulate data using DPPCA model when pilot and covariate data present
#'
#'@param n                  number of simulated samples
#'@param ppca_obj           the results from \code{ppca.metabol}
#'
#' @examples
#' # Simulate 100 samples using the DPPCA model with pilot data
#' # using the Urine Spectra dataset
#'
#' library(MetabolAnalyze)
#' data(UrineSpectra)
#'
#' pilot_data <- UrineSpectra[[1]]
#' ppca_fit <- ppca.metabol(pilot_data, minq = 2, maxq = 2)
#'
#' sim_DPPCA_pilot(100, ppca_fit)

sim_DPPCA_pilot <- function(n, ppca_obj) {

  v       <- 0.1                       # v2.true
  phi     <- 0.8                       # phi.true
  eta_sd  <- sqrt(v / (1 - phi^2))     # eta.sd | SV model on the errors
  alpha_sigma <- 5                     # ao | alpha.sigma <- ao <- 5
                                       # the scale parameter of the prior
                                       # distribution of the variance.

  n_pilot_vars  <- nrow(ppca_obj$loadings)
  n_latent_dims <- ncol(ppca_obj$loadings)
  mu <- rep(0, n_latent_dims)

  rmvnorm_fun <- partial(rmvnorm, method = "svd")

  u <- rnorm(n_latent_dims, mu, eta_sd) %>%
    diag(n_latent_dims) %>%
    rmvnorm_fun(n, mu, .)

  # No pilot_var_means added like in the other models?
  x <- exp(rnorm(1, 0, eta_sd)) %>%
    diag(n_pilot_vars) %>%
    rmvnorm_fun(n, rep(0, n_pilot_vars), .) +
      tcrossprod(u, ppca_obj$loadings)

  return(x)
}



