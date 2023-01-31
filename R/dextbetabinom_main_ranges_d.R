#' Definition and Positivity Range of Extended Deta-Binomial Distribution In Terms of d
#'
#' Returns the intervals where the probability density function of the extended beta-binomial distribution is defined, and where it is positive, in terms of the intermediate variable d.
#' Here we strive to consider as many edge cases carefully as we can, e.g. if x = 0 and size = 0, we aknowledge that the probability value is 1 even if alpha and beta are NA.
#'
#' @param x The number of successful draws.
#' @param size The total number of draws.
#' @param alpha The prior number of successful draws in Polya urn model.
#' @param beta The prior number of unsuccessful draws in Polya urn model.
#' @param AI Optional (default=NA) A parameter to indicate the proportion between alpha in beta in the edge cases like alpha = beta = 0.
#'
#' @return A matrix 2x3 with the basic information about the values of the probability density function as a function of d. The elements result(1,1) and result(1,2) denote the interval of the values of d where the pdf is defined. The elements result(2,1) and result(2,2) denote the interval of the values of d where the pdf is positive. The entry result(1,3) contains the value of the pdf if it does not depend on d.
#'
#' @importFrom stats dbinom
#'
dextbetabinom_main_ranges_d <- function(x, size, alpha, beta, AI=NA){
  total_set = c(-Inf, Inf, NA)
  empty_set = c(Inf,-Inf, NA)

  res = matrix(nrow=2, ncol=3)
  mat_count = x
  pat_count = size-x

  # counts should be defined and non-negative, otherwise this is not a point the distribution
  if (!is.finite(mat_count) | !is.finite(pat_count) | (mat_count < 0) | (pat_count < 0)){
    res[1,] = empty_set
    res[2,] = empty_set
  } else
    # the only case which doesn't care about alpha and beta parameters is (0,0)
    # otherwise the denominator contains alpha+beta
    if((mat_count == 0) & (pat_count == 0)){
      res[1,] = total_set
      res[1,3] = 1
      res[2,] = total_set
    } else
      # so now alpha and beta should be defined, and shouldn't be of different signs
      if(is.na(alpha) | is.na(beta) | (sign(alpha)*sign(beta) == -1)){
        res[1,] = empty_set
        res[2,] = empty_set
      } else
        # alpha and beta are both zero
        if(alpha+beta == 0){
          if(is.na(AI) | (AI < 0) | (AI > 1)){
            # if we don't know AI from different means, NA
            res[1,] = empty_set
            res[2,] = empty_set
          } else {
            res[1,] = total_set
            # the distribution here is AI on (n,0) and (1-AI) on (0,n)
            if(pat_count == 0){
              res[1,3] = AI
              if(AI > 0){
                res[2,] = total_set
              } else {
                res[2,] = empty_set
              }
            } else
              if(mat_count == 0){
                res[1,3] = 1-AI
                if(AI < 1){
                  res[2,] = total_set
                } else {
                  res[2,] = empty_set
                }
              } else {
                res[1,3] = 0
                res[2,] = empty_set
              }
          }
        } else
          if(is.infinite(alpha) & is.infinite(beta)){
            if(is.na(AI) | (AI < 0) | (AI > 1)){
              # if we don't know AI from different means, NA
              res[1,] = empty_set
              res[2,] = empty_set
            } else {
              # the distribution here is binomial with the probability AI
              res[1,] = total_set
              res[1,3] = dbinom(mat_count, mat_count+pat_count, AI)
              if(res[1,3] > 0){
                res[2,] = total_set
              } else {
                res[2,] = empty_set
              }
            }
          } else
            # the case where only one of alpha and beta is infinite doesn't even need AI
            if(is.infinite(alpha)){
              res[1,] = total_set
              if(pat_count == 0){
                res[1,3] = 1
                res[2,] = total_set
              } else {
                res[1,3] = 0
                res[2,] = empty_set
              }
            } else
              if(is.infinite(beta)){
                res[1,] = total_set
                if(mat_count == 0){
                  res[1,3] = 1
                  res[2,] = total_set
                } else {
                  res[1,3] = 0
                  res[2,] = empty_set
                }
              } else
                # at this point, all four parameters are defined, are finite, and alpha+beta != 0
                # therefore the distribution is always well-defined at least for d=0

                # the last case where d doesn't matter
                if(mat_count+pat_count == 1){
                  res[1,] = total_set
                  if(is.na(AI)){
                    AI = alpha/(alpha+beta)
                  }
                  if(mat_count == 1){
                    res[1,3] = AI
                  } else {
                    res[1,3] = 1-AI
                  }
                  if(res[1,3] > 0){
                    res[2,] = total_set
                  } else {
                    res[2,] = empty_set
                  }
                } else {
                  s = sign(alpha+beta)
                  definition_bound = binom_progression_bound(s*(alpha+beta), mat_count+pat_count)
                  if(s == 1){
                    res[1,] = c(definition_bound, Inf, NA)
                  } else {
                    res[1,] = c(-Inf, -definition_bound, NA)
                  }

                  positivity_bound = max(
                    binom_progression_bound(s*(alpha+beta), mat_count+pat_count),
                    binom_progression_bound(s*alpha, mat_count),
                    binom_progression_bound(s*beta, pat_count)
                  )
                  if(positivity_bound == Inf){
                    end = -Inf
                  } else {
                    end = Inf
                  }
                  if(s == 1){
                    res[2,] = c(positivity_bound, end, NA)
                  } else {
                    res[2,] = c(-end, -positivity_bound, NA)
                  }
                }

  return(res)
}
