#' Scale atomic vector to lie between 0 and 1
#'
#' This function takes a vector of data and scales it to
#' lie between the values 0 and 1. It is intended as a helper function
#' and can be used with \code{map}
#'
#' @param x A numeric vector
#'
#' @examples
#' # Standardize a vector
#' x <- rnorm(100, 3, 5)
#' standardize(x)
#'
#' # Standardize the columns of a data frame
#' iris %>%
#'  select(-Species) %>%
#'  purrr::dmap(standardize)

standardize <- function(x) {
  .min <- min(x)
  .max <- max(x)

  (x - .min) / (.max - .min)
}

#' Split data into random groups
#'
#'
#'@param n Number of observations to be grouped
#'@param probs A named vector of probabilities for groups
#'

random_group <- function(n, probs) {
  probs <- probs / sum(probs)
  g <- findInterval(seq(0, 1, length = n), c(0, cumsum(probs)),
                    rightmost.closed = TRUE)
  names(probs)[sample(g)]
}

#' Generate n random splits of a data frame
#'
#' Generate n random splits of a data frame
#'
#' @param df a data frame
#' @inheritParams random_group

partition <- function(df, n, probs) {
  replicate(n, split(df, random_group(nrow(df), probs)), FALSE) %>%
    transpose() %>%
    as_data_frame()
}

#' Sample from the null distribution
#'
#'

sample_dist <- function(x, n_sims, n1, n2){
  probs <- c(A = n1, B = n2)

  partition(x, n_sims, probs) %>% transmutate(
    sd = map2(A, B, ~  sqrt((1/n1 + 1/n2) * ((n1 - 1)*diag(var(.x)) + (n2 - 1)*diag(var(.y))) / (n1 + n2 - 2)))
  )
}
