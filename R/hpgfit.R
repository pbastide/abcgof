#' @title Holdout Frequentist Goodness of Fit Test
#'
#' @description
#' This function performs a post-inference holdout Goodness of Fit (GoF) test,
#' that checks whether we can reject the hypothesis that the `target`
#' is from the same model as the particles from the `sumstat` matrix.
#'
#' It works by first approximating the posterior distribution using a
#' localization parameter `eps` and one of the rejection, local linear or ridge regression
#' methods (see documentation from `abc::abc`).
#' It then draws a fraction `split` of the reference table `sumstat`
#' for calibration, and then compares the score of the
#' `target` with the distribution of scores of the calibration points.
#'
#' @param method one of "rejection" (default), "loclinear" or "ridge".
#' @param kernel for "loclinear" and "ridge", the kernel to be used.
#' @param lambda for "ridge", the lambda penalties to be used in the ridge regression. The median of all predictions with various lambdas will be used. Default to `lambda = c(0.0001, 0.001, 0.01)`.
#' @param param the matrix of parameters matching the sumstat.
#' @param sim.fun a function that, for a matrix of parameters (size nsim x nparams), returns a simulated matrix of sumstats (size nsim x nsumstats).
#' @param ... further arguments to be passed to `sim.fun`.
#' @param param_transform a named vector of function, that to each parameter attributes a transformation function.
#' Each function must be one of "none" (no transformation), "log" or "logit".
#' If "logit", upper and lower values are taken in `param_lower_bound` `param_upper_bound` arguments.
#' If the vector is of length 1, it will be recycled to all parameters.
#' Default to "none": no transformation on any parameter.
#' @param param_lower_bound a named vector of lower bounds for each parameter in `param`.
#' If the vector is of length 1, it will be recycled to all parameters.
#' Default to `-Inf`: no lower bound.
#' @param param_upper_bound a named vector of lower bounds for each parameter in `param`.
#' If the vector is of length 1, it will be recycled to all parameters.
#' Default to `Inf`: no upper bound.
#' @inheritParams freqgfit
#' @inheritParams gfit
#' @param trans.fun A custom transformation function, that takes as input the whole vector of parameters.
#' If specified, it will override the default logit transform. Both `trans.fun` and `back.trans.fun` must be specified.
#' Default to NULL.
#' @param back.trans.fun A custom back transformation function, that takes as input the whole vector of parameters.
#' If specified, it will override the default logit transform. Both `trans.fun` and `back.trans.fun` must be specified.
#' Default to NULL.
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
#' different sets of calibration points from the (fixed) posterior.
#' Highest Density Intervals are then computed using the
#' \code{\link[HDInterval]{hdi}} function from package \code{HDInterval}.
#'
#' @return An object of class `hpgfit`. For each score, it contains a list with:
#'  * `pval` a matrix of p-values for the test data (rows), for each possible value of k (if `nboot=1`),
#'  or an array of such matrices (if `nboot>1`).
#'  * `n.ref` the number of particles in the reference dataset (number of rows of `sumstat`).
#'  * `n.calib` the number of particles used for calibration.
#'  * `eps` the localization fraction.
#'  * `split` the fraction of data used for calibration.
#'  * `n.boot` the number of bootstrap replicates.
#'  * `k_range` the min and max values for k when using the "LOF max" score.
#'
#' @examples
#' data(gaussian_laplace)
#'
#' ## 1000 simulations from the Laplace dataset
#' sumstat <- gaussian_laplace$dataset_laplace$sim
#' ## Matching parameters and simulation function
#' param <- gaussian_laplace$dataset_laplace$param
#' sim_laplace <- gaussian_laplace$dataset_laplace$sim.fun
#'
#' ## 10 observations from the Gaussian dataset
#' target <- gaussian_laplace$dataset_gaussian$sim[1:10, ]
#' ## 10 replicate obsertvations from the Gaussian dataset
#' param_target <- gaussian_laplace$dataset_gaussian$param[1:10, ]
#' target.rep <- gaussian_laplace$dataset_gaussian$sim.fun(param_target)
#'
#' ## Apply hpgfit
#' res <- hpgfit(target, target.rep, param, sumstat,
#'               sim.fun = sim_laplace,
#'               method = "rejection", eps = 0.1) # rejection method
#' res
#'
#' ## Apply hpgfit with bootstrap
#' res <- hpgfit(target, target.rep, param, sumstat,
#'               sim.fun = sim_laplace,
#'               method = "rejection", eps = 0.1,
#'               nboot = 10) # nboot should be larger
#' res
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
#' [summary.hpgfit()] for summarizing the results;
#' [gof.fit.holdout()] for prior GoF computation with fixed calibration points;
#' [gfit()] for pre-inference prior GoF test.
#'
#' @export
#'
hpgfit <- function(target, target.replica,
                   param, sumstat,
                   sim.fun, method = c("rejection", "loclinear", "ridge"),
                   kernel = c("epanechnikov", "rectangular", "gaussian", "triangular", "biweight", "cosine"),
                   lambda = c(0.0001, 0.001, 0.01),
                   param_transform = "none",
                   param_lower_bound = -Inf,
                   param_upper_bound = Inf,
                   trans.fun = NULL,
                   back.trans.fun = NULL,
                   score = c("lof", "kNN"),
                   k = c(2, 3, 5, 8, 14, 20),
                   k_range = range(k),
                   eps = 0.01, split = 0.5,
                   norm = sd, ncores = 1, nboot = 1, ...) {
  ## check dimensions
  target <- check_target(target)
  target.replica <- check_target(target.replica)
  target.replica <- check_sumstat(target, target.replica)
  check_replica(target, target.replica)
  sumstat <- check_sumstat(target, sumstat)
  check_replica(sumstat, param)
  ## Normalization
  n.target <- nrow(target)
  normvec <- rep(1, ncol(sumstat))
  if (is.function(norm)) {
    all_stats <- norm_stats(rbind(target, target.replica), sumstat, norm)
    sumstat <- all_stats$sumstat
    target <- all_stats$target[1:n.target, , drop = FALSE]
    target.replica <- all_stats$target[(n.target+1):(2*n.target), , drop = FALSE]
    normvec <- all_stats$normvec
    rm(all_stats)
  }
  ## check transforms
  if ((is.null(trans.fun) && !is.null(back.trans.fun)) || (!is.null(trans.fun) && is.null(back.trans.fun))) {
    stop("One of trans.fun or back.trans.fun is NULL. Please specify both if you need to use a custom transformation.")
  }
  if (is.null(trans.fun) || is.null(back.trans.fun)) {
    trans_and_back <- check_param_transform(param, param_transform, param_lower_bound, param_upper_bound)
    trans.fun <- trans_and_back$transform
    back.trans.fun <- trans_and_back$back_transform
  }
  ## gfit
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  res <- gof.fit.holdout(target, target.replica, param, sumstat,
                         sim.fun, method, kernel, lambda,
                         trans.fun, back.trans.fun,
                         k, k_range, eps, split, normvec, score, ncores, nboot, ...)
  class(res) <- "hpgfit"
  return(res)
}

#' @export
#' @inheritParams base::print
#' @method print hpgfit
#' @rdname hpgfit
print.hpgfit <- function(x, ...) {
  print(summary(x))
}

#' @export
#' @inheritParams base::summary
#' @param level the level for the confidence interval
#' @method summary hpgfit
#' @rdname hpgfit
summary.hpgfit <- function(object, score = "lof", k = "max", level = 0.95, ...) {
  res <- summary_internal(object, score, k, level, ...)
  class(res) <- "hpgfitsummary"
  return(res)
}

#' @export
#' @inheritParams base::print
#' @method print hpgfitsummary
#' @rdname hpgfit
print.hpgfitsummary <- function(x, ...) {
  # Call
  cat(paste0("Holdout Posterior GoF analysis\n"))
  cat(paste0("Using the ", x$score, " with k = ", x$k_sel, "\n"))
  cat(paste0("Number of lines in the reference table: ", x$n.ref, " ;\n"))
  cat(paste0("Number of lines in the posterior table: ", x$eps * x$n.ref, " (eps=", x$eps * 100, "%);\n"))
  cat(paste0("Number of calibration point: ", x$n.calib, " (split=", x$split * 100, "%);\n"))
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

#' @title "Frequentist" goodness of fit function
#'
#' @description
#' Internal function to compute the GoF using the LOF scores.
#' The data should already be formatted and normalized.
#' Unless you are an advanced user, you should probably use
#' \code{\link{hpgfit}}.
#'
#' @inheritParams gof.fit
#' @inheritParams gof.fit.post
#' @inheritParams gof.fit.freq
#' @inheritParams hpgfit
#' @param data.param the matrix of parameters matching the sumstat.
#' @param normvec the column normalization values.
#' @param trans.fun a transformation function, that takes a matrix of parameters, and return a transformed matrix of parameters.
#' @param back.trans.fun a transformation function, that takes a transformed matrix of parameters, and return a matrix of parameters in the original space.
#' @return For each score, a list, with
#'  * `pval` a matrix of p-values for the test data (rows), for each possible value of k.
#'
#' @export
#'
gof.fit.holdout <- function(data.test.obs, data.test.new,
                            data.param, data.ref,
                            sim.fun, method = c("loclinear", "ridge", "rejection"),
                            kernel = c("epanechnikov", "rectangular", "gaussian", "triangular", "biweight", "cosine"),
                            lambda = c(0.0001, 0.001, 0.01),
                            trans.fun = function(x) return(x),
                            back.trans.fun = function(x) return(x),
                            k = c(2, 3, 5, 8, 14, 20),
                            k_range = range(k),
                            eps = 1 / 100, split = 0.5,
                            normvec,
                            score = c("lof", "kNN"),
                            ncores = 1,
                            nboot = 1, ...) {

  # scores
  score <- match.arg(score, c("lof", "kNN"), several.ok = TRUE)
  score_fun <- NULL
  for (sc in score) {
    if ("lof" == sc) score_fun[[sc]] <- lof.internal
    if ("kNN" == sc) score_fun[[sc]] <- meankNN
  }
  # method
  method <- match.arg(method)
  kernel <- match.arg(kernel)

  n.post <- round(eps * nrow(data.ref)) # number of particles in the "posterior"
  n.calib <- round(n.post * split)
  n.param <- ncol(data.param)

  if (ncores > 1 && nboot == 1) { # use core for pods only if no bootstrap replicates
    cl <- parallel::makeCluster(ncores) # , outfile = ""
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
    ncoresboot <- ncores
  }

  irep <- 0
  pval.test <- foreach::foreach(irep = 1:nrow(data.test.obs), .combine = rbind, .export = ...names()) %dopar% {
    get_lof_freq(data.test.obs[irep, ], data.test.new[irep, ],
                 data.param, data.ref,
                 sim.fun, method, kernel, lambda,
                 trans.fun, back.trans.fun,
                 k, k_range, n.post, n.calib, n.param, normvec,
                 score, score_fun,
                 nboot, ncoresboot,
                 ...)
  }
  if (nboot == 1) {
    pval.test <- format_diagnostic(pval.test, c(k, "max"))
  } else {
    if (nrow(data.test.obs) == 1) pval.test <- t(pval.test)
    pval.test <- apply(pval.test, 2, function(x) format_diagnostic(do.call(rbind, x), c(k, "max")))
    pval.test <- lapply(score, function(sc) sapply(pval.test, function(rr) rr[[sc]], simplify = "array"))
    names(pval.test) <- score
  }
  ## result
  pval.test <- lapply(pval.test, function(x) {colnames(x) <- c(k, "max"); rownames(x) <- NULL; return(list(pval = x, n.ref = nrow(data.ref), n.calib = floor(nrow(data.ref) * eps * split), eps = eps, split = split, k_range = k_range, n.boot = nboot))})
  return(pval.test)
}

get_lof_freq <- function(pods, pnew,
                         data.param, data.ref,
                         sim.fun, method, kernel, lambda,
                         trans.fun, back.trans.fun,
                         k, k_range, n.post, n.calib, n.param, normvec,
                         score, score_fun,
                         nboot, ncoresboot,
                         ...) {
  # post
  dists.pods.ref <- calDIST(pods, data.ref) # distances from pods
  index.best <- order(dists.pods.ref)[1:n.post] # index of posterior particles
  if (method == "rejection") {
    data.param.adj <- data.param[index.best, ]
  } else {
    # weights
    dist <- dists.pods.ref[index.best]
    ds <- max(dist)
    if(kernel == "epanechnikov") weights <- 1 - (dist/ds)^2
    if(kernel == "rectangular") weights <- dist/ds
    if(kernel == "gaussian") weights <- 1/sqrt(2*pi)*exp(-0.5*(dist/(ds/2))^2)
    if(kernel == "triangular") weights <- 1 - abs(dist/ds)
    if(kernel == "biweight") weights <- (1 - (dist/ds)^2)^2
    if(kernel == "cosine") weights <- cos(pi/2*dist/ds)
    # transform
    data.param.trans <- data.param[index.best, ]
    data.param.trans <- trans.fun(data.param.trans)
    if (method == "loclinear") {
      # local linear adjustment
      data.param.adj <- linear_adjust(pods, data.ref[index.best, ], data.param.trans, weights)
    } else if (method == "ridge") {
      data.param.adj <- linear_adjust_ridge(pods, data.ref[index.best, ], data.param.trans, weights, lambda)
    }
    # re-sampling parameters according to weights
    data.param.adj <- data.param.adj[sample(seq_len(n.post), n.post, replace = TRUE, prob = weights), ]
    # back transform
    data.param.adj <- back.trans.fun(data.param.adj)
    colnames(data.param.adj) <- colnames(data.param.trans)
  }
  colnames(data.param.adj) <- colnames(data.param)
  data.ref.posterior <- sim.fun(data.param.adj, ...)
  data.ref.posterior <- norm_vec(data.ref.posterior, normvec)
  ## pvals
  if (nboot == 1) { # no bootstrap replicate
    pval.test <- get_pvals_post(pnew, data.ref.posterior, score, score_fun, k, k_range, n.post, n.calib)
    return(pval.test)
  } else { # with bootstrap
    # parallel
    if (ncoresboot > 1) {
      cl <- parallel::makeCluster(ncoresboot) # , outfile = ""
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
    } else {
      foreach::registerDoSEQ()
    }
    iboot <- 0
    res_boot <- foreach::foreach(iboot = 1:nboot, .export = c("score", "score_fun", "k", "k_range", "n.post", "n.calib")) %dopar% {
      get_pvals_post(pnew, data.ref.posterior, score, score_fun, k, k_range, n.post, n.calib)
    }
    # merge and summarise results
    pval.test <- sapply(res_boot, function(x) as.matrix(x), simplify = "array")
    return(res_boot)
  }
}

get_pvals_post <- function(pnew, data.ref.posterior, score, score_fun, k, k_range, n.post, n.calib) {
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

linear_adjust <- function(pods, data.ref.best, data.param.trans, weights) {
  n.param <- ncol(data.param.trans)
  n.post <- nrow(data.ref.best)
  # local linear adjustment
  fit1 <- lsfit(data.ref.best, data.param.trans, wt = weights)
  # pred <- t(fit1$coefficients) %*% c(1, pods)
  # pred <- rep(1, n.post) %*% t(pred)
  # residuals <- fit1$residuals
  # data.param.adj <- pred + residuals
  coeff_fit <- t(structure(cbind(fit1$coefficients)[fit1$qr$pivot,], names = names(fit1$coefficients)))
  pred <- coeff_fit %*% c(1, pods)
  pred <- matrix(pred, ncol = n.param, nrow = n.post, byrow = TRUE)
  residuals <- data.param.trans - cbind(1, data.ref.best) %*% t(coeff_fit)
  residuals <- cbind(residuals)
  # the_m <- apply(residuals, FUN = mean, 2)
  # residuals <- sapply(1:n.param, FUN = function(x){residuals[, x] - the_m[x]})
  # pred <- sapply(1:n.param, FUN = function(x){pred[, x] + the_m[x]})
  data.param.adj <- pred + residuals
  return(data.param.adj)
}

linear_adjust_ridge <- function(pods, data.ref.best, data.param.trans, weights, lambda, norm = sd) {
  n.param <- ncol(data.param.trans)
  n.post <- nrow(data.ref.best)
  # local ridge linear adjustment
  fitridgenet <- glmnet(data.ref.best, data.param.trans, weights = weights,
                        family = "mgaussian", alpha = 0, lambda = sort(lambda, decreasing = TRUE),
                        standardize = TRUE)
  fitted.values.all <- predict(fitridgenet, newx = data.ref.best)
  pred.all <- predict(fitridgenet, newx = matrix(pods, ncol = length(pods)))
  pred.med <- apply(pred.all, c(1,2), median)
  pred.med <- matrix(pred.med, nrow = n.post, ncol = n.param, byrow = T)
  fitted.values.med <- apply(fitted.values.all, c(1,2), median)
  residuals <- data.param.trans - fitted.values.med # median of fitted values for each accepted point and parameter
  data.param.adj <- pred.med + residuals
  return(data.param.adj)
}

