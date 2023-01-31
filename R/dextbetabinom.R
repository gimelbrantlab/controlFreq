#' The Probability Density Function of the Extended Beta-Binomial Distribution
#'
#' Returns the value of the pdf of the extended beta-binomial distribution. The function will accept the parameters either using the alpha, beta and d parameters, or using the AI and Q parameters.
#'
#' @param x The number of successful draws.
#' @param size The total number of draws.
#' @param alpha Optional (default=NA) The prior number of successful draws in Polya urn model.
#' @param beta Optional (default=NA) The prior number of unsuccessful draws in Polya urn model.
#' @param d Optional (default=NA) The prior number of unsuccessful draws in Polya urn model.
#' @param AI Optional (default=NA) The proportion between alpha and beta.
#' @param Q Optional (default=NA) The overdispersion parameter.
#' @param log.p Optional (default=FALSE) Whether to return the pdf or the logarithm of it. Get the unchanged value of the pdf by default.
#'
#' @return The value of the pdf of the extended beta-binomial distribution.
#'
#' @export
#'
dextbetabinom <- function(x, size,
                          alpha=NA, beta=NA, d=NA,
                          AI=NA, Q=NA,
                          log.p=FALSE){
  if(is.na(alpha) | is.na(beta)){
    alpha = AI*(size-Q)
    beta = (1-AI)*(size-Q)
    d = Q-1
  } else if(!is.na(AI) & (AI*beta != alpha*(1-AI))){
    return(NA)
  }

  ranges = dextbetabinom_main_ranges_d(x, size, alpha, beta, AI)
  defined_range = ranges[1,]
  positive_range = ranges[2,]
  if((defined_range[1] == Inf) & (defined_range[2] == -Inf)){
    return(NA)
  }
  if((defined_range[1] == -Inf) & (defined_range[2] == Inf)){
    if(log.p){
      return(log(defined_range[3]))
    } else {
      return(defined_range[3])
    }
  }
  if(is.na(d) | (d <= defined_range[1]) | (d >= defined_range[2])){
    return(NA)
  }
  if((d <= positive_range[1]) | (d >= positive_range[2])){
    if(log.p){
      return(-Inf)
    } else {
      return(0)
    }
  }
  if(abs(alpha) > abs(beta)){
    a = alpha
    b = beta
    m = x
    p = size-x
  } else {
    a = beta
    b = alpha
    m = size-x
    p = x
  }

  a_part = -sum(vapply(seq_len(m)-1, function (i) {
    log1p(b/(a + i*d))
  }, numeric(1)))

  b_part = -sum(vapply(seq_len(p)-1, function (j) {
    log1p(((j+1)*a - m*b + m*d)/(m + j + 1)/(b + j*d))
  }, numeric(1)))

  logres = a_part+b_part
  if(log.p){
    return(logres)
  } else {
    return(exp(logres))
  }
}
