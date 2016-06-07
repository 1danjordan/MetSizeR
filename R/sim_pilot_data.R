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

 sim_PPCA_data_with_pilot <- function(n_samps, n_pilot_vars, n_covars, n_latent_dims, ppcca_coefs){

   C <- rmvnorm(n, rep(0, n_covars), diag(n_covars)) %>%
     map(standardize) %>%
     rbind(rep(1, n), t(.))

   u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims)) + t(ppcca_coefs %*% C)

 }
