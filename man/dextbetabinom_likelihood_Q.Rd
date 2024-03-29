% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dextbetabinom_likelihood_Q.R
\name{dextbetabinom_likelihood_Q}
\alias{dextbetabinom_likelihood_Q}
\title{The Likelihood using the Extended Beta-Binomial Distribution}
\usage{
dextbetabinom_likelihood_Q(
  mat_counts,
  pat_counts,
  AIs,
  Q,
  log.p = FALSE,
  na.rm = F
)
}
\arguments{
\item{mat_counts}{A vector of the maternal counts.}

\item{pat_counts}{A vector of the paternal counts.}

\item{AIs}{A vector of true allelic imbalance proportions.}

\item{Q}{Optional (default=NA) The overdispersion parameter.}

\item{log.p}{Optional (default=FALSE) Whether to return the likelihood or the logarithm of it. Get the unchanged value of the likelihood by default.}

\item{na.rm}{Optional (default=FALSE) Whether to remove the parameter rows with at least one NA. Do not remove by default.}
}
\value{
The value of the pdf of the extended beta-binomial distribution.
}
\description{
Returns the product of the pdfs of the extended beta-binomial distribution for an array of parameters.
}
