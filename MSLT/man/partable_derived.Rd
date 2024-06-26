% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table.R
\name{partable_derived}
\alias{partable_derived}
\title{Return table with derived model parameters}
\usage{
partable_derived(run, sdrep_bc = NULL, CI_level = 0.95)
}
\arguments{
\item{run}{run returned by \code{\link{fitMSLT}} to extract model parameters from}

\item{sdrep_bc}{the biascorrected sdreport object for the run to extract model parameters from}

\item{CI_level}{Confidence level to accompany the point estimates with}
}
\value{
Table with derived model parameters
}
\description{
Return table with derived model parameters
}
