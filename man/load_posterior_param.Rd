% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MegaLMM_functions.R
\name{load_posterior_param}
\alias{load_posterior_param}
\title{load the posterior samples of a single parameter from all saved chunks}
\usage{
load_posterior_param(MegaLMM_state, param, chunks = NULL)
}
\arguments{
\item{param}{Name of parameter to load}

\item{folder}{folder to find Posterior chunk files}

\item{samples}{vector of sample indices to load. If NULL, all samples loaded}
}
\value{
array with all chuncks of Posterior samples appended together
}
\description{
load the posterior samples of a single parameter from all saved chunks
}
