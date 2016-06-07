 #' Simulate data for PPCA when pilot data is present
 #'
 #' @param n_samps number of samples
 #' @param n_pilot_vars the number of variables in the pilot data
 #' @param n_covars the number of covariates
 #' @param n_latent_dims the number of latent dimensions (set as 2 always)

 # This name of the function along with the argument names are
 # obviously changing  I think there might be a mistake in here,
 # given this is supposed to be PPCA but the PPCCA coefficients
 # end up in there...?

 sim_PPCA_data_with_pilot <- function(n_samps,
                                      n_pilot_vars,
                                      n_covars,
                                      n_latent_dims,
                                      ppcca_coefs) {

   u <- rmvnorm(n, rep(0, n_covars), diag(n_covars)) %>%
     map(standardize) %>%
     rbind(rep(1, n), t(.))

   x <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims)) +
     t(ppcca_coefs %*% u)

   return(x)

 }

 #' Simulate data for PPCCA when pilot data is present
 #'
 #'@param n                  number of samples
 #'@param n_latent_dims      number of latent dimensions
 #'@param ppca_sigma         the posterior mode estimate of the variance
 #'                          of the error terms (MetabolAnalyze docs)
 #'@param ppca_loadings      a loadings matrix
 #'@param pilot_var_means    the means of the variables from the pilot data


 # I think I might pass the whole MetabolAnalyze list object
 # into this function

 sim_PPCCA_data_with_pilot <- function(n, n_pilot_vars,   # p
                                       n_latent_dims,     # q
                                       ppca_sigma,        # sig
                                       ppca_loadings,     # W
                                       pilot_var_means) { # mu

   u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims))

   x <- rmvnorm(n, rep(0, n_pilot_vars), ppca_sigma * diag(n_pilot_vars)) +
     tcrossprod(u, ppca_loadings) +
     matrix(pilot_var_means, n, n_pilot_vars, byrow=TRUE)

   return(x)
 }

# The variable names are a total mess and I don't know what
# any of them refer to or their meaning really so I'm just trying
# to back track them and simplify from there

sim_DPPCA_data_with_pilot <- function(){

  v              <- 0.1                 # v2.true
  phi            <- 0.8                 # phi.true
  alpha_sigma    <- 5                   # ao | alpha.sigma <- ao <- 5
                                        # the scale parameter of the prior
                                        # distribution of the variance.
  mu <- rep(0, n_latent_dims)           # zero vector for rmvnorm mean arg

  # SV model on the errors
  eta_sd   <- sqrt(v / (1 - phi^2))                 # eta.sd
  # a square matrix with values only on the diagonal
  eta_sc_sd <- eta_sd %>% diag(n_latent_dims)


  u_sigma <- rmvnorm(1, mu, eta_sc_sd) %>% diag(n_latent_dims)
  u <- rmvnorm(n, mu, diag(eta_sc_true, n_latent_dims))

  # TODO: There's definitely something fishy with passing
  # the vector b into rgamma().

  b <- c(0.5 * (alpha_sigma - 1), 0.25 * (alpha_sigma - 1)) # bo
  w_sigma <- 1 / rgamma(n_latent_dims, alpha_sigma, b) %>% diag(n_latent_dims)
  w <- rmvnorm(n_pilot_vars, mu, w_sigma)


  x_sigma <- rnorm(1, 0, eta_sd) %>% exp() %>% diag(n_latent_dims)
  x <- rmvnorm(n, mu, x_sigma) + tcrossprod(u, w)

  return(x)
}









