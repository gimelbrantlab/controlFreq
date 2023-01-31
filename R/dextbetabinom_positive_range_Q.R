#' Positivity Range of Extended Beta-Binomial Distribution in Terms of Q
#'
#' Returns the interval of Q where the probability distribution function of the extended beta-binomial distribution is positive. Useful later for the gradient descend procedures.
#'
#' @param x The number of successful draws.
#' @param size The total number of draws.
#' @param AI The proportion between alpha and beta.
#'
#' @return A tuple representing the infimum and the supremum of all values of Q where the pdf is positive.
#'
dextbetabinom_positive_range_Q <- function(x, size, AI){
  d = dextbetabinom_main_ranges_d(x, size, AI, 1-AI)[2,1]
  if(d == Inf){
    return(c(Inf, -Inf))
  } else if(d == -Inf){
    return(c(-Inf, Inf))
  } else {
    return(c(1+d*(size-1)/(1+d), size))
  }
}
