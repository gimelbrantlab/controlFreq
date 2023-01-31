#' Best Fitting Q According to the Extended Beta-Binomial Model
#'
#' Finds a Q maximizing the likelihood using several attempts of gradient descend. One attempt may not succeed because of data rows with low coverage. Then we automatically filter a few rows with the lowest coverage and repeat the attempt. We do so until an attempt of the gradient descend succeeds.
#'
#' @param allele_counts A dataframe with the maternal and paternal counts.
#' @param AI_vect A vector of true allelic imbalance proportions.
#' @param inf_Q Optional (default=1) Only use the data rows for which inf_Q yields a well-defined distribution. The default inf_Q=1 simplifies to the binomial distribution, therefore it does not restrict the input.
#' @param sup_Q Optional (default=1) Only use the data rows for which sup_Q yields a well-defined distribution. The default sup_Q=1 does not restrict the input. The parameter can be used to effectively filter the rows with a low coverage, e.g. setting sup_Q=30 filters out the data rows with the coverage less than 30.
#'
#' @return The fitted value of Q.
#'
#' @export
#'
dextbetabinom_max_likelihood_Q <- function(allele_counts, AI_vect,
                                           inf_Q=1, sup_Q=1){
  # Note: it might be slowed down due to low sup_Q,
  # you may choose sup_Q=30, for example, and it'll start with coverage 30.

  df = data.frame(mat_count = allele_counts[, 1],
                  pat_count = allele_counts[, 2],
                  AI = AI_vect)

  df[,c("inf", "sup")] = t(sapply(seq_len(nrow(df)), function(i){
    drow = df[i,]
    return(
      dextbetabinom_positive_range_Q(x=drow$mat_count, size=drow$mat_count+drow$pat_count, drow$AI)
    )
  }))

  # here we remove all genes with the positive range (-Inf, Inf)
  # since their derivative is 0
  df = df[(df$inf > -Inf) | (df$sup < Inf), ]

  # now we apply the restrictions
  df = df[(df$inf < inf_Q) & (df$sup > sup_Q), ]
  if(nrow(df) == 0){
    print("initial inf_Q and sup_Q are too strict, no genes left")
    return(NA)
  }
  inf_Q = max(df$inf, -1)
  sup_Q = min(df$sup, 100)

  starting_Q = 1
  derivative_to_diff_const = 1/10
  status = "start"
  while(TRUE){
    if(status == "accurate"){
      return(Q)
    }
    if(status != "start"){
      #print(paste("got status", status))
      #flush.console()
    }
    if(status == "inf_bump"){
      # this probably comes from the situation where
      # there is no local maximum near the lowest bound
      return(-Inf)
    }
    if(status == "sup_bump"){
      if(is.finite(min(df$sup))){
        sup_Q = 1.1*min(df$sup)+1
        df = df[df$sup > sup_Q,]
        if(nrow(df) == 0){
          return(Inf)
        }
      } else {
        return(Inf)
      }
    }
    if(status == "const"){
      return(Q)
    }
    if(status == "max_step"){
      starting_Q = Q
      derivative_to_diff_const = derivative_to_diff_const/2
    }
    #print(df)
    res = dextbetabinom_max_likelihood_Q_attempt(df,
                                                 starting_Q=starting_Q, inf_Q=inf_Q, sup_Q=sup_Q,
                                                 derivative_to_diff_const=derivative_to_diff_const,
                                                 streak_base=1.1, oob_tax=0.95,
                                                 max_step=100, final_accuracy=1/40000)
    Q = res$Q
    status = res$status
  }

  return(Q)
}
