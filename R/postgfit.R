#' @title Posterior goodness of fit function
#'
#' @description
#' This function is likely to be DEPRECATED in stable versions of the package.
#'
#' This function performs a prior Goodness of Fit (GoF) test,
#' that checks whether we can reject the hypothesis that the `target`
#' is from the same model as the particles from the `sumstat` matrix.
#'
#' It works by sub-sampling `nb.replicate` points from a localized subset of the
#' reference table `sumstat`, and then by comparing the score of the
#' `target` with the distribution of scores of the calibration points.
#'
#' @param sumstat.replica a matrix of replicate summary statistics from the training set,
#' with dimension Ntrain x Nstats.
#' Each line i of the replicate summary statistic should have been simulated
#' with the same parameters as the one used for the same line i of the reference table.
#' @param eps the proportion of data points used to approximate the posterior distribution
#' of the reference table around each observation point. Default to 0.01.
#' @inheritParams gfit
#'
#' @return For each score, it a list with:
#'  * `score.calib` a matrix of LOF scores for the calibration data (rows), for each possible value of k.
#'  * `score.test` a matrix of LOF scores for the test data (rows), for each possible value of k.
#'  * `pval` a matrix of p-values for the test data (rows), for each possible value of k.
#'  * `n.ref` the number of particles in the reference dataset (number of rows of `sumstat`).
#'  * `n.calib` the number of particles used for calibration (equal to `nb.replicate`).
#'  * `k_range` the min and max values for k when using the "LOF max" score.
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
#' [gfit()] for pre-inference prior GoF test.
#'
#' @export
#'
postgfit <- function(target, sumstat, sumstat.replica,
                     nb.replicate = 100,
                     score = c("lof", "kNN"),
                     k = c(2, 3, 5, 8, 14, 20),
                     k_range = range(k),
                     eps = 0.01,
                     norm = sd,
                     ncores = 1) {
  ## check dimensions
  target <- check_target(target)
  sumstat <- check_sumstat(target, sumstat)
  sumstat.replica <- check_sumstat(target, sumstat.replica)
  check_replica(sumstat, sumstat.replica)
  if (nb.replicate >= nrow(sumstat)) stop("sumstat must have more rows as nb.replicate.")
  ## Normalization
  n.target <- nrow(target)
  if (is.function(norm)) {
    all_stats <- norm_stats(rbind(target, sumstat.replica), sumstat, norm)
    sumstat <- all_stats$sumstat
    target <- all_stats$target[1:n.target, , drop = FALSE]
    sumstat.replica <- all_stats$target[-(1:n.target), , drop = FALSE]
    rm(all_stats)
  }
  ## Split ref
  calibid <- sample(seq_len(nrow(sumstat)), size = nb.replicate)
  data.calib <- sumstat[calibid, ]
  data.ref <- sumstat[-calibid, ]
  data.ref.replicat <- sumstat.replica[-calibid, ]
  ## gfit
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  return(gof.fit.post(target, data.calib, data.ref, data.ref.replicat, k, k_range, eps, score, ncores))
}

#' @title "Posterior" goodness of fit function
#'
#' @description
#' Internal function to compute the posterior GoF.
#' The data should already be formatted and normalized.
#' Unless you are an advanced user, you should probably use
#' \code{\link{postgfit}}.
#'
#' @inheritParams gof.fit
#' @inheritParams postgfit
#' @param data.ref.replica a normalized matrix of reference replicate summary statistics
#' from the training set, with dimension Ntrain x Nstats.
#' Each line i of the replicate summary statistic should have been simulated
#' with the same parameters as the one used for the same line i of the reference table.
#'
#' @inherit postgfit return
#'
#' @export
#'
gof.fit.post <- function(data.test, data.calib, data.ref, data.ref.replica,
                         k = c(2, 3, 5, 8, 14, 20),
                         k_range = range(k),
                         eps = 1 / 100,
                         score = c("lof", "kNN"),
                         ncores = 1) {
  # scores
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  score_fun <- NULL
  for (sc in score) {
    if ("lof" == sc) score_fun[[sc]] <- lof.internal
    if ("kNN" == sc) score_fun[[sc]] <- meankNN
  }
  # gof post function
  n.calib = nrow(data.calib)
  n.post <- round(eps * nrow(data.ref)) # number of particles in the "posterior"
  get_lof_post <- function(pods) {
    dists.pods.ref <- calDIST(pods, data.ref) # distances from pods
    index.best <- order(dists.pods.ref)[1:n.post] # index of posterior particles
    data.ref.posterior <- data.ref.replica[index.best, ]
    res <- sapply(score, function(sc) score_fun[[sc]](pods, data.ref.posterior, k))
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
    return(get_lof_post(data.calib[irep, ]))
  }
  diagnostic.calib <- format_diagnostic(diagnostic.calib, k)
  diagnostic.calib <- add_score_max(diagnostic.calib, k_range)
  # test
  diagnostic.test <- foreach::foreach(irep = 1:nrow(data.test), .combine = rbind) %dopar% {
    return(get_lof_post(data.test[irep, ]))
  }
  diagnostic.test <- format_diagnostic(diagnostic.test, k)
  diagnostic.test <- add_score_max(diagnostic.test, k_range)
  # pvals
  pval.test <- mapply(get_p_val, diagnostic.calib, diagnostic.test, SIMPLIFY = FALSE)
  ## result
  tmpfun <- function(x,y,z) {
    res <- list(score.calib = x, score.test = y, pval = z, n.ref = nrow(data.ref), n.calib = nrow(x), eps = eps, k_range = k_range)
    colnames(res$score.calib) <- colnames(res$score.test) <- colnames(res$pval) <- c(k, "max")
    return(res)
  }
  res <- mapply(tmpfun,
                diagnostic.calib, diagnostic.test, pval.test,
                SIMPLIFY = FALSE)
  return(res)
}
