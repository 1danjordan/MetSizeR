#' Simulate data for PPCCA when pilot data is present
#'
#'@param n                  number of samples
#'@param n_latent_dims      number of latent dimensions
#'@param ppca_sigma         the posterior mode estimate of the variance
#'                          of the error terms (MetabolAnalyze docs)
#'@param ppca_loadings      a loadings matrix
#'@param pilot_var_means    the means of the variables from the pilot data

sim_PPCA_data_with_pilot <- function(n, ppca_obj, pilot_var_means) {

  n_pilot_vars  <- nrow(ppca_obj$loadings)
  n_latent_dims <- ncol(ppca_obj$loadings)

  u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims))

  x <- rmvnorm(n, rep(0, n_pilot_vars), ppca_obj$sig * diag(n_pilot_vars)) +
    tcrossprod(u, ppca_obj$loadings) +
    matrix(pilot_var_means, n, n_pilot_vars, byrow=TRUE)

  return(x)
}

#' Simulate data for PPCA when pilot data is present
 #'
 #' @param n_samps number of samples
 #' @param n_pilot_vars the number of variables in the pilot data
 #' @param n_covars the number of covariates
 #' @param n_latent_dims the number of latent dimensions (set as 2 always)

 # This name of the function along with the argument names are  changing

 sim_PPCCA_data_with_pilot <- function(n_samps,
                                      n_pilot_vars,
                                      n_covars,
                                      n_latent_dims,
                                      ppca_obj) {

   u <- rmvnorm(n, rep(0, n_covars), diag(n_covars)) %>%
     map(standardize) %>%
     rbind(rep(1, n), t(.))

   u <- rmvnorm(n, rep(0, n_latent_dims), diag(n_latent_dims)) +
     t(ppca_obj$coefs %*% u)

   x <- rmvnorm(n, rep(0, n_pilot_vars), ppca_sigma * diag(n_pilot_vars)) +
     tcrossprod(u, ppca_loadings) +
     matrix(pilot_var_means, n, n_pilot_vars, byrow=TRUE)

   return(x)

 }







