#' Simulate data for DPPCA without pilot data
#'
#' @param n              number of samples to be simulated
#' @param n_pilot_vars   number of variables in the pilot data
#' @param n_latent_dims  number of latent dimension

sim_pilot_DPPCA <- function(n, n_pilot_vars, n_latent_dims) {

  v       <- 0.1                       # v2.true
  phi     <- 0.8                       # phi.true
  mu      <- rep(0, n_latent_dims)     # zero vector for rmvnorm mean arg
  eta_sd  <- sqrt(v / (1 - phi^2))     # eta.sd | SV model on the errors
  alpha_sigma <- 5                     # ao | alpha.sigma <- ao <- 5
                                       # the scale parameter of the prior
                                       # distribution of the variance.

  rmvnorm_fun <- partial(rmvnorm, method = "svd")

  # TODO: There's definitely something fishy with passing
  # the vector b into rgamma().
  b <- c(0.5 * (alpha_sigma - 1), 0.25 * (alpha_sigma - 1))

  u <- rnorm(n_latent_dims, mu, eta_sd) %>%
    diag(n_latent_dims) %>%
    rmvnorm_fun(n, mu, .)

  w <- (1 / rgamma(n_latent_dims, alpha_sigma, b)) %>%
    diag(n_latent_dims) %>%
    rmvnorm_fun(n_pilot_vars, mu, .)

  x <- exp(rnorm(1, 0, eta_sd)) %>%
    diag(n_pilot_vars) %>%
    rmvnorm_fun(n, rep(0, n_pilot_vars), .) + tcrossprod(u, w)

  return(x)

}
