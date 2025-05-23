% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{fitMSLT}
\alias{fitMSLT}
\title{fitMSLT}
\usage{
fitMSLT(data, par, conf, rel.tol = 1e-10, map = setMap(conf, par), ...)
}
\arguments{
\item{data}{data list returned by \code{\link{setData}}}

\item{par}{parameter list  returned by \code{\link{setPar}}}

\item{conf}{model configurations list  returned by \code{\link{defConf}}}

\item{rel.tol}{rel.tol in nlminb}

\item{map}{map object used by TMB to couple and remove parameters}

\item{...}{additional parameters passed to sdreport}
}
\value{
fitted model
}
\description{
fitMSLT
}
