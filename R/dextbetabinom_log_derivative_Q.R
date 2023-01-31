#' The Log Derivative of the Extended Beta-Binomial Distribution
#'
#' Returns the log derivative of the pdf of the extended beta-binomial distribution with regards to Q.
#'
#' @param x The number of successful draws.
#' @param size The total number of draws.
#' @param AI The proportion between alpha and beta.
#' @param Q The overdispersion parameter.
#' @param force_calculation Optional (default=FALSE) Whether to skip the sanity checks for whether the parameters are in the definition range and calculate the derivative directly. Do the checks by default.
#'
#' @return The value of the logarithm of the pdf of the extended beta-binomial distribution.
#'
dextbetabinom_log_derivative_Q <- function(x, size, AI, Q,
                                           force_calculation=FALSE){
  if(!force_calculation){
    positive_range = dextbetabinom_positive_range_Q(x, size, AI)
    if(positive_range[1] == -Inf){
      return(0)
    } else
      if(is.na(Q) | (Q <= positive_range[1]) | (Q >= positive_range[2])){
        return(NA)
      }
  }

  a_part = (size-1)*(1-AI)*sum(vapply(seq_len(x)-1, function (i) {
    i/((i-AI)*Q + (AI*size-i))/((i-1)*Q + (size-i))
  }, numeric(1)))

  b_part = (size-1)*sum(vapply(seq_len(size-x)-1, function (j) {
    (j*AI - x*(1-AI))/((j+AI-1)*Q + (size*(1-AI)-j))/((x+j-1)*Q + (size-x-j))
  }, numeric(1)))

  return(a_part+b_part)
}
