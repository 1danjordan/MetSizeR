#' Scale atomic vector to lie between 0 and 1
#'
#' This function takes an atomic vector of data and scales it to
#' lie between the values 0 and 1. It is intended as a helper function
#' and can be used with \code{map}
#'
#' @param x An atomic vector
#' @examples
#' x <- rnorm(100, 3, 5)
#' standardize(x)
#'
#' iris %>%
#'  select(-Species) %>%
#'  dmap(standardize)

standardize <- function(x) {
  .min <- min(x)
  .max <- max(x)

  (x - .min) / (.max - .min)
}
