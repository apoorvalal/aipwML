% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fns.R
\name{ate_aipw}
\alias{ate_aipw}
\title{compute AIPW}
\usage{
ate_aipw(d, psrange = c(0, 1))
}
\arguments{
\item{d}{output from fit_me function}

\item{psrange}{2-vector of pscore to trim within}
}
\value{
AIPW estimate of ATE
}
\description{
compute AIPW
}
