#' Overdispersion measures (iQCCs)
#'
#' Calculates iQCC values for selected samples in allelic table.
#'
#' @param df Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param reps Optional (default=NA, all replicates), a vector of replicate numbers for which the analysis should be applied
#' @param inf_Q Optional (default=1) Only use the data rows for which inf_Q yields a well-defined distribution. The default inf_Q=1 simplifies to the binomial distribution, therefore it does not restrict the input.
#' @param sup_Q Optional (default=1) Only use the data rows for which sup_Q yields a well-defined distribution. The default sup_Q=1 does not restrict the input. The parameter can be used to effectively filter the rows with a low coverage, e.g. setting sup_Q=30 filters out the data rows with the coverage less than 30.
#'
#' @return A list of iQCC triplets for each sample within selected: lower, upper and exact estimate.
#'
#' @export
#'
compute_iQCC_for_selected_samples <- function(df, reps, inf_Q=1, sup_Q=1){
  res = lapply(reps, function(r){
    pack_in = merge(thresholdingCounts(df, reps = r), countsToAI(df = df, reps = reps), by = "ID")
    pack_ex = merge(thresholdingCounts(df, reps = r), countsToAI(df = df, reps = reps[reps!=r]), by = "ID")
    iQCC_in = sqrt(dextbetabinom_max_likelihood_Q(allele_counts = pack_in[,2:3], AI_vect = pack_in[,4], sup_Q = sup_Q))
    iQCC_ex = sqrt(dextbetabinom_max_likelihood_Q(allele_counts = pack_ex[,2:3], AI_vect = pack_ex[,4], sup_Q = sup_Q))
    print(sqrt(iQCC_in*iQCC_ex))
    return(list(
      iQCC_inclused_inAI = iQCC_in,
      iQCC_excluded_inAI = iQCC_ex,
      iQCC_geom_mean = sqrt(iQCC_in*iQCC_ex)
    ))
  })
  names(res) = reps
  return(res)
}
