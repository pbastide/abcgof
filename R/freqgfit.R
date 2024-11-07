#' @title Pre-Computed Holdout Frequentist Goodness of Fit Test
#'
#' @description
#' This function is similar to [hpgfit()], but it only works for the
#' rejection method, and for datasets with pre-computed replicates.
#'
#' It is likely to be DEPRECATED in stable versions of the package.
#'
#' This function performs a post-inference holdout Goodness of Fit (GoF) test,
#' that checks whether we can reject the hypothesis that the `target`
#' is from the same model as the particles from the `sumstat` matrix.
#'
#' It works by first approximating the posterior distribution using a
#' localization parameter `eps` the rejection
#' method (see documentation from `abc::abc`).
#' It then draws a fraction `split` of the reference table `sumstat`
#' for calibration, and then compares the score of the
#' `target` with the distribution of scores of the calibration points.
#'
#' @param split proportion of the posterior to be used for calibration. Default to 0.5.
#' @param target.replica a matrix of observations, with dimension Nobs x Nstats.
#' Each line i of the replicate target should come from the same model
#' as the one used for the same line i of the reference table.
#' @inheritParams postgfit
#'
#' @return For each score, a list with:
#'  * `pval` a matrix of p-values for the test data (rows), for each possible value of k (if `nboot=1`),
#'  or an array of such matrices (if `nboot>1`).
#'  * `n.ref` the number of particles in the reference dataset (number of rows of `sumstat`).
#'  * `n.calib` the number of particles used for calibration.
#'  * `eps` the localization fraction.
#'  * `split` the fraction of data used for calibration.
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
#' [hpgfit()] for post-inference holdout GoF test.
#'
#' @export
#'
freqgfit <- function(target, target.replica, sumstat, sumstat.replica,
                     score = c("lof", "kNN"),
                     k = c(2, 3, 5, 8, 14, 20),
                     k_range = range(k),
                     eps = 0.01, split = 0.5,
                     norm = sd, ncores = 1) {
  ## check dimensions
  target <- check_target(target)
  target.replica <- check_target(target.replica)
  target.replica <- check_sumstat(target, target.replica)
  check_replica(target, target.replica)
  sumstat <- check_sumstat(target, sumstat)
  sumstat.replica <- check_sumstat(target, sumstat.replica)
  check_replica(sumstat, sumstat.replica)
  ## Normalization
  n.target <- nrow(target)
  if (is.function(norm)) {
    all_stats <- norm_stats(rbind(target, target.replica, sumstat.replica), sumstat, norm)
    sumstat <- all_stats$sumstat
    target <- all_stats$target[1:n.target, , drop = FALSE]
    target.replica <- all_stats$target[(n.target+1):(2*n.target), , drop = FALSE]
    sumstat.replica <- all_stats$target[-(1:(2*n.target)), ]
    rm(all_stats)
  }
  ## gfit
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  return(gof.fit.freq(target, target.replica, sumstat, sumstat.replica, k, k_range, eps, split, score, ncores))
}

#' @title "Frequentist" goodness of fit function
#'
#' @description
#' Internal function to compute the GoF using the LOF scores.
#' The data should already be formatted and normalized.
#' Unless you are an advanced user, you should probably use
#' \code{\link{freqgfit}}.
#'
#' @inheritParams gof.fit
#' @inheritParams gof.fit.post
#' @inheritParams freqgfit
#' @param data.test.obs a normalized matrix of observations, with dimension Nobs x Nstats.
#' @param data.test.new a normalized matrix of new (replicated) observations, with dimension Nobs x Nstats.
#' @param ncores number of cores for parallel computations.
#'
#' @return For each score, a list, with
#'  * `pval` a matrix of p-values for the test data (rows), for each possible value of k.
#'
#' @export
#'
gof.fit.freq <- function(data.test.obs, data.test.new,
                         data.ref, data.ref.replica,
                         k = c(2, 3, 5, 8, 14, 20),
                         k_range = range(k),
                         eps = 1 / 100, split = 0.5,
                         score = c("lof", "kNN"),
                         ncores = 1) {

  # scores
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  score_fun <- NULL
  for (sc in score) {
    if ("lof" == sc) score_fun[[sc]] <- lof.internal
    if ("kNN" == sc) score_fun[[sc]] <- meankNN
  }

  n.post <- round(eps * nrow(data.ref)) # number of particles in the "posterior"
  n.calib <- round(n.post * split)

  get_lof_freq <- function(i) {
    pods <- data.test.obs[i, ]
    pnew <- data.test.new[i, ]
    # post
    dists.pods.ref <- calDIST(pods, data.ref) # distances from pods
    index.best <- order(dists.pods.ref)[1:n.post] # index of posterior particles
    data.ref.posterior <- data.ref.replica[index.best, ]
    # split
    index.split <- sample(n.post, n.calib)
    data.ref.posterior.split.calib <- data.ref.posterior[index.split,]
    data.ref.posterior.split.cloud <- data.ref.posterior[-index.split,]
    Tab.dbscan.post <- dbscan::kNN(data.ref.posterior.split.cloud, max(k)) #Pre-computation
    # lof
    diagnostic.calib <- sapply(score, function(sc) apply(data.ref.posterior.split.calib, 1,
                                                           function(pcalib) score_fun[[sc]](pcalib, data.ref.posterior.split.cloud, k, Tab.dbscan.post$dist, Tab.dbscan.post$id)))
    diagnostic.calib <- format_diagnostic(diagnostic.calib, k)
    diagnostic.calib <- add_score_max(diagnostic.calib, k_range)
    diagnostic.test <- sapply(score, function(sc) score_fun[[sc]](pnew, data.ref.posterior.split.cloud, k,
                                                                  Tab.dbscan.post$dist, Tab.dbscan.post$id))
    diagnostic.test <- format_diagnostic(diagnostic.test, k)
    diagnostic.test <- add_score_max(diagnostic.test, k_range)
    # pval
    pval.test <- as.data.frame(mapply(get_p_val_vec, diagnostic.calib, diagnostic.test, SIMPLIFY = FALSE))
    return(pval.test)
  }

  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores) # , outfile = ""
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }

  irep <- 0
  pval.test <- foreach::foreach(irep = 1:nrow(data.test.obs), .combine = rbind) %dopar% {
    get_lof_freq(irep)
  }
  pval.test <- format_diagnostic(pval.test, c(k, "max"))

  ## result
  pval.test <- lapply(pval.test, function(x) {colnames(x) <- c(k, "max"); rownames(x) <- NULL; return(list(pval = x, n.ref = nrow(data.ref), n.calib = nrow(x), eps = eps, split = split, k_range = k_range))})
  return(pval.test)
}
