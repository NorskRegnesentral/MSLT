% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{setData}
\alias{setData}
\title{setData}
\usage{
setData(d, predAreaUTM, conf)
}
\arguments{
\item{d}{Data provided by user}

\item{predAreaUTM}{prediction area provided by user}

\item{conf}{model configurations list  returned by \code{\link{defConf}}}
}
\value{
data list used by TMB
}
\description{
Define data used internally based on data porvided and configuration file
}
