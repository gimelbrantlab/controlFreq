#' Bound for Binomial Progression
#'
#' Returns the bound for x allowed in the expression 1/(start * (start+x) * ... * (start+(count-1)*x)).
#' A small subfunction to make the calculations in the next function a bit easier.
#' The edge cases are as follows. If the count is equal to zero, all x are allowed, so the infimum is -Inf.
#' If the count is greater than zero, and the starting term is equal to zero, then all x are disallowed.
#' Then the infimum is Inf.
#'
#' @param start Non-negative real number, the starting term of the arithmetic progression in the denominator.
#' @param count The number of terms in the denominator.
#'
#' @return The infimum of all values of x for which the expression is well-defined, i.e. the denominator is non-zero.
#'
binom_progression_bound <- function(start, count) {
  # we expect non-negative start and integer non-negative count
  if(count == 0){
    return(-Inf)
  } else if(start == 0){
    return(Inf)
  } else {
    return(-start/(count-1))
  }
}
