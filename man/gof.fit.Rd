% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfit.R
\name{gof.fit}
\alias{gof.fit}
\title{"Prior" goodness of fit function}
\usage{
gof.fit(
  data.test,
  data.ref,
  data.calib,
  k = c(2, 3, 5, 8, 14, 20),
  k_range = range(k),
  score = c("lof", "kNN"),
  ncores = 1,
  kdist = NULL,
  kid = NULL
)
}
\arguments{
\item{data.test}{a normalized matrix of observations, with dimension Nobs x Nstats.}

\item{data.ref}{a normalized matrix of reference summary statistics from the training set, with dimension Ntrain x Nstats.}

\item{data.calib}{a normalized matrix of reference summary statistics from the training set, with dimension Ncalib x Nstats.}

\item{k}{The kth-distance to be used to calculate LOFs of kNNs.
k can be a vector which contains multiple k values based on which the score needs to be calculated.}

\item{k_range}{A vector of min and max values for k when using the "LOF max" score (see Details).}

\item{score}{the score to use for calling outliers. Can be "lof", "kNN", or both.}

\item{ncores}{number of cores for parallel computations.}

\item{kdist}{(optional) a pre-computed matrix of k-distances:
\code{kdist[i,k]} should contain the distance of particle i in the \code{data.ref} matrix
to its k-th neighbor. See \code{\link[dbscan]{kNN}}.}

\item{kid}{(optional) a pre-computed matrix of k-ids:
\code{kid[i,k]} should contain the id of the k-th -neighbor of particle in the \code{data.ref} matrix.
See \code{\link[dbscan]{kNN}}.}
}
\value{
An object of class \code{gfit}. For each score, it contains a list with:
\itemize{
\item \code{score.calib} a matrix of LOF scores for the calibration data (rows), for each possible value of k (if \code{nboot=1}).
\item \code{score.test} a matrix of LOF scores for the test data (rows), for each possible value of k (if \code{nboot=1}).
\item \code{pval} a matrix of p-values for the test data (rows), for each possible value of k (if \code{nboot=1}),
or an array of such matrices (if \code{nboot>1}).
\item \code{n.ref} the number of particles in the reference dataset (number of rows of \code{sumstat}).
\item \code{n.calib} the number of particles used for calibration (equal to \code{nb.replicate}).
\item \code{n.boot} the number of bootstrap replicates.
\item \code{k_range} the min and max values for k when using the "LOF max" score.
}
}
\description{
Internal function to compute the GoF using the lof and mean kNN scores.
The data should already be formatted and normalized.
Unless you are an advanced user, you should probably use
\code{\link{gfit}}.
}
\details{
Matrices \code{kdist} and \code{kid} can be obtained through a call to
\code{\link[dbscan]{kNN}} on \code{data.ref}.
These are only needed for the lof score.
}
