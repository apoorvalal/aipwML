% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fns.R
\name{ate_ipw}
\alias{ate_ipw}
\title{estimate ATE by IPW}
\usage{
ate_ipw(est, w, y, df, fml, ret_est = T, psrange = c(0, 1), ...)
}
\arguments{
\item{est}{name of estimator}

\item{w}{name of treatment in dataframe}

\item{y}{name of outcome}

\item{df}{dataframe}

\item{fml}{formula for propensity score}

\item{ret_est}{boolean for whether to return ATE estimate, returns pscore vectors when false}

\item{psrange}{range of propensity scores to use for ATE estimation}
}
\value{
vector of treatment effects or ATE estimate
}
\description{
estimate ATE by IPW
}
