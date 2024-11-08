% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_p_val}
\alias{get_p_val}
\title{Get p value from diagnostic}
\usage{
get_p_val(diagnostic.calib, diagnostic.test)
}
\arguments{
\item{diagnostic.calib}{a matrix of calibration scores, size Ncalib x Nscores.}

\item{diagnostic.test}{a matrix of test scores, size Ntest x Nscores.}
}
\value{
A matrix of p values, size Ntest x Nscores.
}
\description{
Computes the empirical p values for a set of scores.
}
\keyword{internal}
