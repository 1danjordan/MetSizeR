#' Simulating data using pilot data
#'
#' \code{sim_pilot_data} simulates data using pilot data provided by the user
#'

sim_pilot_data <- function


#' Simulating data using pilot data with the PPCA model
#'
#'
#'
#' 1 = PPCA; 2 = PPCCA; 3 = DPPCA

# PPCA
   C <- rmvnorm(n, ZeroL, IL)
   # Standardize covariates for stability
   C <- standardize(C)
   C <- rbind(rep(1, n), t(C))
   u <- rmvnorm(n, Zeroq, Iq) + t(Alpha %*% C)
   }#ifppca

 x <- rmvnorm(n, Zerop, sig * Ip) +
   tcrossprod(u, W) +
   matrix(mu, n, p, byrow=TRUE)

# If PPCCA
 u <- rmvnorm(n,Zeroq,Iq)
 x <- rmvnorm(n, Zerop, sig * Ip) + tcrossprod(u, W) + matrix(mu, n, p, byrow=TRUE)

 #' Simulate data for PPCA when pilot data is present
 #'
 #' @param n_samps number of samples
 #' @param n_pilot_vars the number of variables in the pilot data
 #' @param n_covars the number of covariates
 #' @param n_latent_dims the number of latent dimensions (set as 2 always)

 # This name of the function along with the argument names are obviously changing
 # I think there might be a mistake in here, given this is supposed to be PPCA but
 # the PPCCA coefficients end up in there...?

 sim_PPCA_data_with_pilot <- function(n_samps,
                                      n_pilot_vars,
                                      n_covars,
                                      n_latent_dims,
                                      ppcca_coefs) {

   C <- rmvnorm(n, rep(0, n_covars), diag(n_covars)) %>%
     map(standardize) %>%
     rbind(rep(1, n), t(.))

   u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims)) +
     t(ppcca_coefs %*% C)

 }

 #' Simulate data for PPCCA when pilot data is present
 #'
 #'@param n                  number of samples
 #'@param n_latent_dims      number of latent dimensions
 #'@param ppca_sigma         the posterior mode estimate of the variance
 #'                          of the error terms (from the MetabolAnalyze docs)
 #'@param ppca_loadings      a loadings matrix
 #'@param pilot_var_means    the means of the variables from the pilot data

 sim_PPCCA_data_with_pilot <- function(n, n_pilot_vars,   # p
                                       n_latent_dims,     # q
                                       ppca_sigma,        # sig
                                       ppca_loadings,     # W
                                       pilot_var_means) {  # mu

   u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims))

   x <- rmvnorm(n, rep(0, n_pilot_vars), ppca_sigma * diag(n_pilot_vars)) +
     tcrossprod(u, ppca_loadings) +
     matrix(pilot_var_means, n, n_pilot_vars, byrow=TRUE)

   return(x)
 }












