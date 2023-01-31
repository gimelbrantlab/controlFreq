#' The Likelihood using the Extended Beta-Binomial Distribution
#'
#' Returns the product of the pdfs of the extended beta-binomial distribution for an array of parameters.
#'
#' @param mat_counts A vector of the maternal counts.
#' @param pat_counts A vector of the paternal counts.
#' @param AIs A vector of true allelic imbalance proportions.
#' @param Q Optional (default=NA) The overdispersion parameter.
#' @param log.p Optional (default=FALSE) Whether to return the likelihood or the logarithm of it. Get the unchanged value of the likelihood by default.
#' @param na.rm Optional (default=FALSE) Whether to remove the parameter rows with at least one NA. Do not remove by default.
#'
#' @return The value of the pdf of the extended beta-binomial distribution.
#'
dextbetabinom_likelihood_Q <- function(mat_counts, pat_counts, AIs, Q, log.p = FALSE, na.rm = F){
  logres = sum(sapply(seq_along(along.with = AIs), function(i){
    dextbetabinom(x=mat_counts[i], size=mat_counts[i]+pat_counts[i], AI=AIs[i], Q=Q, log.p=TRUE)
  }), na.rm=na.rm)
  if(log.p){
    return(logres)
  } else {
    return(exp(logres))
  }
}
