#' @title Prior Goodness of Fit test
#'
#' @description
#' This function performs a prior Goodness of Fit (GoF) test,
#' that checks whether we can reject the hypothesis that the `target`
#' is from the same model as the particles from the `sumstat` matrix.
#'
#' It works by sub-sampling `nb.replicate` points from the
#' reference table `sumstat`, and then by comparing the score of the
#' `target` with the distribution of scores of the calibration points.
#'
#' @param target a matrix of observations, with dimension Nobs x Nstats.
#' @param sumstat a matrix of summary statistics from the training set,
#' with dimension Ntrain x Nstats.
#' @param norm the normalization function. Default to \code{\link{sd}}.
#' @param k The kth-distance to be used to calculate LOFs of kNNs.
#' k can be a vector which contains multiple k values based on which the score needs to be calculated.
#' @param k_range A vector of min and max values for k when using the "LOF max" score (see Details).
#' @param nb.replicate number of replicates to draw from the training set for the calibration of the p value.
#' @param score the score to use for calling outliers. Can be "lof", "kNN", or both.
#' @param ncores number of cores for parallel computations.
#' @param nboot number of bootstrap replicates (see Details). Default to 1: no bootstrap replicate.
#'
#' @details
#' The LOF and kNN scores are computations rely on the
#' function \code{\link[dbscan]{kNN}} from the `dbscan` package.
#'
#' `target` can be a matrix, in which case the GoF test is performed independently
#' for each observation vector (rows of the matrix).
#'
#' When using the "LOF max" statistics, the maximum LOF score over `k_range`
#' is taken for all calibration and observation points.
#'
#' When `nboot>1`, the test is performed several times, by drawing `nboot`
#' different sets of `nb.replicates` calibration points.
#' Highest Density Intervals are then computed using the
#' \code{\link[HDInterval]{hdi}} function from package \code{HDInterval}.
#'
#' @return An object of class `gfit`. For each score, it contains a list with:
#'  * `score.calib` a matrix of LOF scores for the calibration data (rows), for each possible value of k (if `nboot=1`).
#'  * `score.test` a matrix of LOF scores for the test data (rows), for each possible value of k (if `nboot=1`).
#'  * `pval` a matrix of p-values for the test data (rows), for each possible value of k (if `nboot=1`),
#'  or an array of such matrices (if `nboot>1`).
#'  * `n.ref` the number of particles in the reference dataset (number of rows of `sumstat`).
#'  * `n.calib` the number of particles used for calibration (equal to `nb.replicate`).
#'  * `n.boot` the number of bootstrap replicates.
#'  * `k_range` the min and max values for k when using the "LOF max" score.
#'
#' @examples
#' data(gaussian_laplace)
#' ## 1000 simulations from the Laplace dataset
#' sumstat <- gaussian_laplace$dataset_laplace$sim
#' ## 10 observations from the Gaussian dataset
#' target <- gaussian_laplace$dataset_gaussian$sim[1:10, ]
#'
#' ## Apply gfit
#' resgfit <- gfit(target, sumstat)
#' resgfit
#'
#' ## Apply gfit with bootstrap
#' resgfit <- gfit(target, sumstat, nboot = 10) # nboot should be larger
#' resgfit
#'
#' @references
#' Le Mailloux G., Bastide, P., Marin, J-M., Estoup, A. (2024+),
#' Goodness of Fit for Bayesian Generative Models with Applications in Population Genetics.
#'
#' Markus M. Breunig, Hans-Peter Kriegel, Raymond T. Ng, and Jörg Sander. 2000.
#' LOF: identifying density-based local outliers. SIGMOD Rec. 29, 2 (June 2000), 93–104.
#' https://doi.org/10.1145/335191.335388
#'
#' @seealso
#' [summary.gfit()] for summarizing the results;
#' [gof.fit()] for prior GoF computation with fixed calibration points;
#' [hpgfit()] for post-inference holdout GoF test.
#'
#' @export
#'
gfit <- function(target, sumstat,
                 nb.replicate = 100,
                 score = c("lof", "kNN"),
                 k = c(2, 3, 5, 8, 14, 20),
                 k_range = range(k),
                 norm = sd,
                 ncores = 1,
                 nboot = 1) {
  ## check dimensions
  target <- check_target(target)
  sumstat <- check_sumstat(target, sumstat)
  if (nb.replicate >= nrow(sumstat)) stop("sumstat must have more rows than nb.replicate.")
  ## Normalization
  if (is.function(norm)) {
    all_stats <- norm_stats(target, sumstat, norm)
    sumstat <- all_stats$sumstat
    target <- all_stats$target
    rm(all_stats)
  }
  ## scores
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  ## Boot
  if (nboot == 1) { # no bootstrap replicate
    ## Split ref
    calibid <- sample(seq_len(nrow(sumstat)), size = nb.replicate)
    data.calib <- sumstat[calibid, ]
    data.ref <- sumstat[-calibid, ]
    res <- gof.fit(target, data.ref, data.calib, k, k_range, score, ncores)
    res <- lapply(res, function(rr) c(rr, n.boot = 1))
  } else { # with bootstrap
    # parallel
    if (ncores > 1) {
      cl <- parallel::makeCluster(ncores) # , outfile = ""
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
    } else {
      foreach::registerDoSEQ()
    }
    iboot <- 0
    res_boot <- foreach::foreach(iboot = 1:nboot) %dopar% {
      ## Split ref
      calibid <- sample(seq_len(nrow(sumstat)), size = nb.replicate)
      data.calib <- sumstat[calibid, ]
      data.ref <- sumstat[-calibid, ]
      return(gof.fit(target, data.ref, data.calib, k, k_range, score, 1))
    }
    # merge and summarise results
    res <- lapply(score, function(sc) list(pval = sapply(res_boot, function(rr) rr[[sc]]$pval, simplify = "array"),
                                            n.ref = nrow(sumstat), n.calib = nb.replicate, k_range = k_range, n.boot = nboot))
    names(res) <- score
  }
  class(res) <- "gfit"
  return(res)
}

#' @export
#' @inheritParams base::print
#' @method print gfit
#' @rdname gfit
print.gfit <- function(x, ...) {
  print(summary(x))
}

#' @export
#' @inheritParams base::summary
#' @param level the level for the confidence interval
#' @method summary gfit
#' @rdname gfit
summary.gfit <- function(object, score = "lof", k = "max", level = 0.95, ...) {
  res <- summary_internal(object, score, k, level, ...)
  class(res) <- "gfitsummary"
  return(res)
}

summary_internal <-  function(object, score = "lof", k = "max", level = 0.95, ...) {
  if (!(score %in% names(object))) stop("The requested score is not in the fitted object.")
  res <- object[[score]]
  res$k <- colnames(res$pval)
  k_ind <- which(k == res$k)
  if (length(k_ind) == 0) stop("The requested k is not in the fitted object.")
  if (res$n.boot > 1) {
    pval <- t(apply(res$pval[, k_ind, , drop = FALSE], 1, summarise_function, level = level))
  } else {
    resbis <- compute_pval_sd(list(score = res))
    pval <- summarise_function_sd(resbis$score$pval[, k_ind], resbis$score$pval_sd[, k_ind], level)
  }
  res <- c(res[!(names(res) == "pval")], list(pval = pval))
  res$score <- score
  res$level <- level
  res$k_sel <- k
  return(res)
}

summarise_function <- function(x, level) {
  hpi <- unname(HDInterval::hdi(as.vector(x), credMass = level))
  return(c(median = median(x),
           lwr = hpi[1],
           upr = hpi[2]))
}

summarise_function_sd <- function(x, sdx, level) {
  qq <- qnorm((1 + level) / 2)
  return(cbind(estim = x,
           lwr = pmax(0, x - qq * sdx),
           upr = pmin(1, x + qq * sdx)))
}

#' @export
#' @inheritParams base::print
#' @method print gfitsummary
#' @rdname gfit
print.gfitsummary <- function(x, ...) {
  # Call
  cat(paste0("Prior GoF analysis\n"))
  cat(paste0("Using the ", x$score, " with k = ", x$k_sel, "\n"))
  cat(paste0("Number of lines in the reference table: ", x$n.ref, " ;\n"))
  cat(paste0("Number of calibration point: ", x$n.calib, " ;\n"))
  cat(paste0("Range of k values in the original object: k in [", paste(x$k_range, collapse = ", "), "] ;\n"))
  cat(paste0("Number of target observations: ", nrow(x$pval), " ;\n"))
  if (x$n.boot > 1) {
    cat(paste0("Table of median, lower and upper ", x$level*100, "% HPD pvalues on ", x$n.boot, " bootstrap replicates:\n"))
    print(x$pval)
  } else {
    cat(paste0("Table of estimate, lower and upper ", x$level*100, "% confidance interval pvalues:\n"))
    print(x$pval)
    cat(paste0("Confidence intervals are based on asymptotic standard error estimation."))
  }
}

#' @title "Prior" goodness of fit function
#'
#' @description
#' Internal function to compute the GoF using the lof and mean kNN scores.
#' The data should already be formatted and normalized.
#' Unless you are an advanced user, you should probably use
#' \code{\link{gfit}}.
#'
#' @inheritParams gfit
#' @param data.test a normalized matrix of observations, with dimension Nobs x Nstats.
#' @param data.ref a normalized matrix of reference summary statistics from the training set, with dimension Ntrain x Nstats.
#' @param data.calib a normalized matrix of reference summary statistics from the training set, with dimension Ncalib x Nstats.
#' @param kdist (optional) a pre-computed matrix of k-distances:
#' \code{kdist[i,k]} should contain the distance of particle i in the `data.ref` matrix
#' to its k-th neighbor. See \code{\link[dbscan]{kNN}}.
#' @param kid (optional) a pre-computed matrix of k-ids:
#' \code{kid[i,k]} should contain the id of the k-th -neighbor of particle in the `data.ref` matrix.
#' See \code{\link[dbscan]{kNN}}.
#'
#' @details
#' Matrices `kdist` and `kid` can be obtained through a call to
#' \code{\link[dbscan]{kNN}} on `data.ref`.
#' These are only needed for the lof score.
#'
#' @inherit gfit return
#'
#' @export
#'
gof.fit <- function(data.test, data.ref, data.calib,
                    k = c(2, 3, 5, 8, 14, 20),
                    k_range = range(k),
                    score = c("lof", "kNN"),
                    ncores = 1,
                    kdist = NULL, kid = NULL) {
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  # precomputations
  if (is.null(kdist) && "lof" %in% score) {
    Tab.dbscan <- dbscan::kNN(data.ref, max(k))
    kdist <- Tab.dbscan$dist
    kid <- Tab.dbscan$id
    rm(Tab.dbscan)
  }
  # score functions
  score_fun <- NULL
  for (sc in score) {
    if ("lof" == sc) score_fun[[sc]] <- lof.internal
    if ("kNN" == sc) score_fun[[sc]] <- meankNN
  }
  # gof function
  n.calib = nrow(data.calib)
  get_lof_prior <- function(pods) {
    res <- sapply(score, function(sc) score_fun[[sc]](pods, data.ref, k, kdist, kid))
    return(res)
  }
  # parallel
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores) # , outfile = ""
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }
  irep <- 0
  # calib
  diagnostic.calib <- foreach::foreach(irep = 1:nrow(data.calib), .combine = rbind) %dopar% {
    return(get_lof_prior(data.calib[irep, ]))
  }
  diagnostic.calib <- format_diagnostic(diagnostic.calib, k)
  diagnostic.calib <- add_score_max(diagnostic.calib, k_range)
  # test
  diagnostic.test <- foreach::foreach(irep = 1:nrow(data.test), .combine = rbind) %dopar% {
    return(get_lof_prior(data.test[irep, ]))
  }
  diagnostic.test <- format_diagnostic(diagnostic.test, k)
  diagnostic.test <- add_score_max(diagnostic.test, k_range)
  # pvals
  pval.test <- mapply(get_p_val, diagnostic.calib, diagnostic.test, SIMPLIFY = FALSE)
  ## result
  tmpfun <- function(x,y,z) {
    res <- list(score.calib = x, score.test = y, pval = z, n.ref = nrow(data.ref), n.calib = nrow(x), k_range = k_range)
    colnames(res$score.calib) <- colnames(res$score.test) <- colnames(res$pval) <- c(k, "max")
    return(res)
  }
  res <- mapply(tmpfun,
                diagnostic.calib, diagnostic.test, pval.test,
                SIMPLIFY = FALSE)
  return(res)
}
