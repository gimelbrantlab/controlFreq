#' The Log Likelihood Derivative using the Extended Beta-Binomial Distribution
#'
#' Returns the log likelihood derivative using the extended beta-binomial distribution for an array of parameters.
#'
#' @param mat_counts A vector of the maternal counts.
#' @param pat_counts A vector of the paternal counts.
#' @param AIs A vector of true allelic imbalance proportions.
#' @param Q The overdispersion parameter.
#' @param force_calculation Optional (default=FALSE) Whether to skip the sanity checks for whether the parameters are in the definition range and calculate the derivative directly. Do the checks by default.
#' @param na.rm Optional (default=FALSE) Whether to remove the parameter rows with at least one NA. Do not remove by default.
#'
#' @return The value of the pdf of the extended beta-binomial distribution.
#'
dextbetabinom_log_likelihood_derivative_Q <- function(mat_counts, pat_counts, AIs, Q, force_calculation=FALSE, na.rm = F){
  sum(sapply(seq_along(along.with = AIs), function(i){
    dextbetabinom_log_derivative_Q(x=mat_counts[i], size=mat_counts[i]+pat_counts[i], AI=AIs[i], Q=Q, force_calculation=force_calculation)
  }), na.rm=na.rm)
}
