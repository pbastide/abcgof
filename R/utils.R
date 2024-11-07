#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
NULL

#' @importFrom stats mad predict quantile sd lsfit coef median qnorm
#' @importFrom glmnet glmnet predict.glmnet
#'
NULL

#' @title Normalize the reference table
#'
#' @description
#' This function applies the `norm` function to each of the Ntrain columns
#' of the training `sumstat` matrix. It then normalize each column of both the
#' `target` and `sumstat` matrices by dividing them by the norm:
#'
#' \code{target[, j] = target[, j] / norm(sumstat[, j])} and
#' \code{sumstat[, j] = sumstat[, j] / norm(sumstat[, j])}
#'
#' @param target a matrix of observations, with dimension Nobs x Nstats.
#' @param sumstat a matrix of summary statistics from the training set,
#' with dimension Ntrain x Nstats.
#' @param norm the normalization function. Default to \code{\link{sd}}.
#'
#' @return A list, with the normalized `target` and `sumstat` matrices.
#'
#' @export
#'
norm_stats <- function(target, sumstat, norm = sd) {
  normvec <- apply(sumstat, 2, norm)
  N.stat <- ncol(sumstat)
  for (j in 1:N.stat) {
    if (normvec[j] > 0) {
      sumstat[,j] <- sumstat[,j] / normvec[j]
      target[,j] <- target[,j] / normvec[j]
    }
  }
  return(list(target = target,
              sumstat = sumstat,
              normvec = normvec))
}

norm_vec <- function(sumstat, normvec) {
  for (j in seq_along(normvec)) {
    if (normvec[j] > 0) {
      sumstat[,j] <- sumstat[,j] / normvec[j]
    }
  }
  return(sumstat)
}

#' @title Get p value from diagnostic
#'
#' @description
#' Computes the empirical p values for a set of scores.
#'
#' @param diagnostic.calib a matrix of calibration scores, size Ncalib x Nscores.
#' @param diagnostic.test a matrix of test scores, size Ntest x Nscores.
#'
#' @return A matrix of p values, size Ntest x Nscores.
#'
#' @keywords internal
#'
get_p_val <- function(diagnostic.calib, diagnostic.test) {
  pval.test <- matrix(NA, nrow(diagnostic.test), ncol(diagnostic.test))
  for(i in seq_len(nrow(diagnostic.test))) {
    for(j in seq_len(ncol(diagnostic.test))){
      pval.test[i, j] <- mean(diagnostic.calib[, j] > diagnostic.test[i, j])
    }
  }
  return(pval.test)
}

#' @title Get p value from diagnostic
#'
#' @description
#' Computes the empirical p values for a set of scores.
#'
#' @param diagnostic.calib a matrix of calibration scores, size Ncalib x Nscores.
#' @param diagnostic.test a vector of test scores, length Nscores.
#'
#' @return A vector of p values of length Nscores.
#'
#' @keywords internal
#'
get_p_val_vec <- function(diagnostic.calib, diagnostic.test) {
  pval.test <- rep(NA, length(diagnostic.test))
  for(j in seq_along(diagnostic.test)) {
    pval.test[j] <- mean(diagnostic.calib[, j] > diagnostic.test[j])
  }
  return(pval.test)
}

#' @title Format a diagnostic table
#'
#' @description
#' Format a diagnostic table, and add the max score.
#'
#' @param diagnostic.table a diagnostic table, each column with a different diagnostics.
#' @param scores the different scores per diagnostic.
#'
#' @return The formatted diagnostic table
#'
#' @keywords internal
#'
format_diagnostic <- function(diagnostic.table, scores) {
  # format
  diagnostic.table <- apply(diagnostic.table, 2, function(x) t(matrix(x, nrow = length(scores))), simplify = FALSE)
  diagnostic.table <- lapply(diagnostic.table, function(x) {colnames(x) <- scores; return(x)})
  return(diagnostic.table)
}

#' @title Format a diagnostic table
#'
#' @description
#' Format a diagnostic table, and add the max score.
#'
#' @param diagnostic.table a diagnostic table, each column with a different diagnostics.
#' @param k_range the range of k values on which to restrict the max
#'
#' @return The formatted diagnostic table
#'
#' @keywords internal
#'
add_score_max <- function(diagnostic.table, k_range) {
  # add max
  diagnostic.table <- mapply(cbind,
                             diagnostic.table,
                             lapply(diagnostic.table, function(x) apply(x[, get_range(x, k_range), drop = FALSE], 1, max)),
                             SIMPLIFY = FALSE)
  diagnostic.table <- lapply(diagnostic.table, function(x) {colnames(x)[length(colnames(x))] <- "max"; return(x)})
  return(diagnostic.table)
}

get_range <- function(x, k_range) {
  return(as.numeric(colnames(x)) <= k_range[2] & as.numeric(colnames(x)) >= k_range[1])
}

#' @title Compute power
#'
#' @description
#' Computes the power of a test.
#'
#' @param gfit_obj an object from the \code{\link{gfit}} or \code{\link{postgfit}} functions.
#' @param level the level of the test. Default to 0.05.
#'
#' @return The same object, with an added power slot.
#'
#' @export
#'
compute_power <- function(gfit_obj, level = 0.05) {
  # compute power
  pow <- lapply(gfit_obj, function(x) apply(x$pval, 2, function(y) mean(y < level)))
  powname <- paste0("pow_", level*100, "%")
  # bind with object
  for (sc in names(gfit_obj)) {
    gfit_obj[[sc]][[powname]] <- pow[[sc]]
  }
  return(gfit_obj)
}

#' @title Compute power
#'
#' @description
#' Computes the approximate standard deviation of pvalues
#'
#' @param gfit_obj an object from the \code{\link{gfit}} or \code{\link{postgfit}} functions.
#'
#' @return The same object, with an added pval_sd slot.
#'
#' @export
#'
compute_pval_sd <- function(gfit_obj) {
  # compute sd
  pvalsd <- lapply(gfit_obj, function(x) sqrt(x$pval * (1 - x$pval) / x$n.calib))
  sdname <- paste0("pval_sd")
  # bind with object
  for (sc in names(gfit_obj)) {
    gfit_obj[[sc]][[sdname]] <- pvalsd[[sc]]
  }
  return(gfit_obj)
}

#' @title Compute power
#'
#' @description
#' Computes the approximate standard deviation of the power.
#' This estimate ignores the variability of the pvalues themselves, and is probably too small.
#'
#' @param gfit_obj an object from the \code{\link{gfit}} or \code{\link{postgfit}} functions.
#'
#' @return The same object, with an added pow_sd slot.
#'
#' @export
#'
compute_power_sd <- function(gfit_obj) {
  message("Conditional standard deviation of the power estimate, ignoring the variability of pvalues. This estimate is too narrow.")
  # compute sd
  powsd <- lapply(gfit_obj, function(x) lapply(x[grep("pow_", names(x))], function(y) sqrt(y * (1 - y) / nrow(x$pval))))
  # bind with object
  for (sc in names(gfit_obj)) {
    for (pp in names(gfit_obj[[sc]])[grepl("pow_", names(gfit_obj[[sc]]))]) {
      gfit_obj[[sc]][[paste0("sd_", pp)]] <- powsd[[sc]][[pp]]
    }
  }
  return(gfit_obj)
}


#' @title Check target
#'
#' @description
#' Format the target matrix.
#'
#' @param target the target matrix
#'
#' @return The formatted target matrix.
#'
#' @keywords internal
#'
check_target <- function(target) {
  if (is.vector(target)) return(t(target))
  return(target)
}

#' @title Check sumstat
#'
#' @description
#' Format the sumstat matrix.
#'
#' @param target the target matrix
#' @param sumstat the sumstat matrix
#'
#' @return The formatted sumstat matrix.
#'
#' @keywords internal
#'
check_sumstat <- function(target, sumstat) {
  if (ncol(target) != ncol(sumstat)) stop("target and sumstat must have the same number of columns.")
  if (is.null(colnames(target)) || is.null(colnames(sumstat))) {
    warning("Columns of target and or sumstat are unnamed. I'm assuming that they were given in the correct order.")
  } else {
    targstat <- match(colnames(target), colnames(sumstat))
    if (anyNA(targstat)) stop("target and sumstat must have the same column names.")
    if (!all(targstat == seq_len(ncol(target)))) {
      warning("target and sumstat have the same column names, but not in the same order. I am re-ordering them.")
      sumstat <- sumstat[, targstat, drop = FALSE]
    }
  }
  return(sumstat)
}

#' @title Check replicate
#'
#' @description
#' Check a matrix and its replicate.
#'
#' @param target the target matrix
#' @param target.replica the sumstat matrix
#'
#' @return nothing
#'
#' @keywords internal
#'
check_replica <- function(target, target.replica) {
  if (nrow(target) != nrow(target.replica)) stop("A table and its replicate must have the same number of lines.")
}

#' @title logit transform
#'
#' @param x vector of values
#' @param a lower bound
#' @param b upper bound
#'
#' @return transformed values
#'
#' @keywords internal
#'
logit <- function(x, a = 0, b = 1) {
  if (a >= b) stop("Lower bound a must be lower than upper bound b.")
  if (any(x <= a) || any(x >= b)) stop("All values of x must be in the interval [a,b].")
  x_norm <- (x - a) / (b - a)
  y <- log(x_norm / (1 - x_norm))
  return(y)
}

#' @title logistic transform
#'
#' @param x vector of values
#' @param a lower bound
#' @param b upper bound
#'
#' @return transformed values
#'
#' @keywords internal
#'
logistic <- function(x, a = 0, b = 1) {
  if (a >= b) stop("Lower bound a must be lower than upper bound b.")
  y_norm <- 1 / (1 + exp(-x))
  y <- (b - a) * y_norm + a
  return(y)
}

#' @title logit transform
#'
#' @param x vector of values
#' @param a lower bound
#' @param b vector of upper bound
#'
#' @return transformed values
#'
#' @keywords internal
#'
logitVec <- function(x, a = 0, b = rep(1, length(x))) {
  # if (any(a >= b)) stop("Lower bound a must be lower than upper bound b.")
  # if (any(x <= a) || any(x >= b)) stop("All values of x must be in the interval [a,b].")
  x_norm <- (x - a) / (b - a)
  y <- log(x_norm / (1 - x_norm))
  y[is.infinite(y)] <- sign(y[is.infinite(y)]) * .Machine$double.xmax^(0.5)
  return(y)
}

#' @title logistic transform
#'
#' @param x vector of values
#' @param a lower bound
#' @param b vector of upper bound
#'
#' @return transformed values
#'
#' @keywords internal
#'
logisticVec <- function(x, a = 0, b = rep(1, length(x))) {
  # if (any(a >= b)) stop("Lower bound a must be lower than upper bound b.")
  y_norm <- 1 / (1 + exp(-x))
  y <- (b - a) * y_norm + a
  return(y)
}


#' @title find transform functions
#'
#' @inheritParams hpgfit
#'
#' @return a list, with:
#'  * `transform` a vector of transform for each parameter
#'  * `back_transform` a vector of back transform for each parameter.
#'
#' @keywords internal
#'
check_param_transform <- function(param, param_transform, param_lower_bound, param_upper_bound) {
  # check format
  names_param_transform <- names(param_transform)
  param_transform <- match.arg(param_transform, c("none", "log", "logit"), several.ok = T)
  names(param_transform) <- names_param_transform
  param_transform <- check_one_param_arg(param, param_transform, "param_transform")
  param_lower_bound <- check_one_param_arg(param, param_lower_bound, "param_lower_bound")
  param_upper_bound <- check_one_param_arg(param, param_upper_bound, "param_upper_bound")
  # check bounds
  if (any(param_lower_bound > param_upper_bound, na.rm = TRUE)) stop("param lower bounds must be smaller than upper bounds.")
  if (any(param_transform == "logit" & (is.infinite(param_lower_bound) | is.infinite(param_upper_bound)))) stop("For the logit transform, bound must be finite values.")
  if (any(param_transform == "log" & (param_lower_bound != 0 | param_upper_bound != Inf))) stop("For the log transform, the lower bound must be 0, and upper bound Inf.")
  # build transform
  transf_par <- function(x) {
    return(sapply(1:ncol(x), function(i) return(switch(param_transform[i],
                                                       "none" = x[, i],
                                                       "log" = log(x[, i]),
                                                       "logit" = logit(x[, i], a = param_lower_bound[i], b = param_upper_bound[i])))))
  }
  back_transf_par <- function(x) {
    return(sapply(1:ncol(x), function(i) return(switch(param_transform[i],
                                                       "none" = x[, i],
                                                       "log" = exp(x[, i]),
                                                       "logit" = logistic(x[, i], a = param_lower_bound[i], b = param_upper_bound[i])))))
  }
  return(list(transform = transf_par,
              back_transform = back_transf_par))
}

check_one_param_arg <- function(param, param_arg, param_arg_name) {
  if (length(param_arg) == 1) {
    param_arg <- rep(param_arg, ncol(param))
    if (!is.null(colnames(param))) {
      names(param_arg) <- colnames(param)
    } else {
      names(param_arg) <- colnames(param) <- 1:ncol(param)
    }
  }
  if (ncol(param) != length(param_arg)) stop(paste0("Argument ", param_arg_name, " must be a single value, or a vector with as many entries as the number of parameters in `param`."))
  if (is.null(colnames(param)) || is.null(names(param_arg))) {
    warning(paste0("Columns of param and or ", param_arg_name, " are unnamed. I'm assuming that they were given in the correct order."))
  } else {
    targstat <- match(colnames(param), names(param_arg))
    if (anyNA(targstat)) stop(paste0("param and ", param_arg_name, " must have the same names."))
    if (!all(targstat == seq_len(ncol(param)))) {
      warning(paste0("param and ", param_arg_name, " have the same column names, but not in the same order. I am re-ordering them."))
      param_arg <- param_arg[targstat]
    }
  }
  return(param_arg)
}
