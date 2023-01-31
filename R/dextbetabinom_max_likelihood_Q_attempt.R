#' Attempt To Find Q Maximizing The Likelihood
#'
#' Finds a maximizing Q using one attempt at gradient descend.
#'
#' @param df A dataframe with the maternal and paternal counts and allelic imbalance proportions.
#' @param starting_Q Starting Q for the graident descent.
#' @param inf_Q The lowest Q allowed for the gradient descent.
#' @param sup_Q The highest Q allowed for the gradient descent.
#' @param derivative_to_diff_const The constant for transitioning from the derivative to the step of the gradient descend.
#' @param streak_base The exponent base speeding up gradient descent going consecutively in one direction.
#' @param oob_tax The tax on any step going out of the predefined bounds.
#' @param max_step The bound on the amount of steps in the gradient descent.
#' @param final_accuracy The target accuracy to finish the gradient step.
#'
#' @importFrom stats median
#'
#' @return The list with the final Q of the gradient descent and the success status.
#'
dextbetabinom_max_likelihood_Q_attempt <- function(df,
                                                   starting_Q, inf_Q, sup_Q,
                                                   derivative_to_diff_const, streak_base, oob_tax,
                                                   max_step, final_accuracy){
  # initializing main loop variables
  Q = starting_Q
  step = 0
  streak = 0
  bound_bumps = 0

  # print(df)
  # print(sup_Q)

  # calibrating derivative to diff constant
  probes = 5
  derivatives = sapply(1:probes, function(i){
    dextbetabinom_log_likelihood_derivative_Q(df$mat_count, df$pat_count, df$AI, (i*inf_Q + (probes+1-i)*sup_Q)/(probes+1), force_calculation=TRUE)
  })
  # print(derivatives)
  median_derivative_value = median(abs(derivatives))
  if(median_derivative_value == 0){
    return(list(Q=starting_Q, status="const"))
  } else {
    diff = sign(derivatives[abs(derivatives) == median_derivative_value])
  }

  d2d_const = derivative_to_diff_const/median_derivative_value

  status = NA
  while(abs(diff) >= final_accuracy) {
    step = step+1
    if(step >= max_step){
      return(list(Q=Q, status="max_step"))
    }

    derivative = dextbetabinom_log_likelihood_derivative_Q(df$mat_count, df$pat_count, df$AI, Q, force_calculation=TRUE)

    derivative_moment = sign(diff)*sign(derivative)
    if(derivative_moment == -1){
      status = "accurate"
    }
    if(sign(derivative_moment) == sign(streak)){
      streak = streak+derivative_moment
    } else {
      streak = derivative_moment
    }
    diff = d2d_const*(Q-inf_Q)*(streak_base**streak)*derivative

    if(Q+diff < inf_Q){
      if(is.na(status)){
        status = "inf_bump"
      }
      diff = oob_tax*(inf_Q-Q)
    } else if(Q+diff > sup_Q){
      if(is.na(status)){
        status = "sup_bump"
      }
      diff = oob_tax*(sup_Q-Q)
    }

    #print(paste("Q", Q, "derivative", derivative, "diff", diff))
    Q = Q + diff
  }
  if(is.na(status)){
    status = "accurate"
  }
  return(list(Q=Q, status=status))
}
