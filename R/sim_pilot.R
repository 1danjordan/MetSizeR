
#' Simulate data without pilot data using a PPCA model
#'
#' @param n                 the number of samples to be simulated
#' @param n_spectral_bins   the number of spectral bins considered

sim_PPCA <- function(n, n_spectral_bins) {

  n_latent_vars <- 2

  # the scale parameter of the prior distribution of the variance
  alpha.sigma   <- 5
  # the shape parameter of the prior distribution of the variance
  beta.sigma    <- 2 * (alpha.sigma - 1)

  bo <- c(0.5 * (alpha.sigma - 1), 0.25 * (alpha.sigma - 1))
  sig <- 1 / rgamma(1, alpha.sigma, beta.sigma)

  u <- rmvnorm(n, rep(0, n_latent_vars), diag(n_latent_vars))
  v <- 1 / (rgamma(n_latent_vars, ao, bo))
  w <- rmvnorm(n_spectral_bins, rep(0, n_latent_vars), diag(v))

  x <- rmvnorm(n, rep(0, n_latent_vars), diag(n_latent_vars)) + tcrossprod(u, w)

  return(x)
}


#' Simulate data for DPPCA without pilot data
#'
#' @param n              number of samples to be simulated
#' @param n_pilot_vars   number of variables in the pilot data
#' @param n_latent_dims  number of latent dimension

sim_DPPCA <- function(n, n_pilot_vars, n_latent_dims) {

  v       <- 0.1                       # v2.true
  phi     <- 0.8                       # phi.true
  eta_sd  <- sqrt(v / (1 - phi^2))     # eta.sd | SV model on the errors
  alpha_sigma <- 5                     # ao | alpha.sigma <- ao <- 5
                                       # the scale parameter of the prior
                                       # distribution of the variance.

  rmvnorm_fun <- partial(rmvnorm, method = "svd")

  # TODO: There's definitely something fishy with passing
  # the vector b into rgamma().
  b <- c(0.5 * (alpha_sigma - 1), 0.25 * (alpha_sigma - 1))

  u <- rnorm(n_latent_dims, rep(0, n_latent_dims), eta_sd) %>%
    diag(n_latent_dims) %>%
    rmvnorm_fun(n, rep(0, n_latent_dims), .)

  w <- (1 / rgamma(n_latent_dims, alpha_sigma, b)) %>%
    diag(n_latent_dims) %>%
    rmvnorm_fun(n_pilot_vars, rep(0, n_latent_dims), .)

  x <- exp(rnorm(1, 0, eta_sd)) %>%
    diag(n_pilot_vars) %>%
    rmvnorm_fun(n, rep(0, n_pilot_vars), .)

  return(x + tcrossprod(u, w))

}
