test_that("Test holdout posterior GoF LOF computations", {

  ## Data
  nsimu <- 1001
  n <- 150
  n.moments <- 2
  n.noise <- 10

  set.seed(20200430)
  mu <- runif(nsimu,-10,10)
  sigma2 <- runif(nsimu,1,4)
  comp_momments <- function(x){c(mean(x),var(x))}
  sumgaussian <- t(sapply(1:nsimu,function(x){comp_momments(rnorm(n,mu[x],sd=sqrt(sigma2[x])))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumgaussian <- cbind(sumgaussian,noise)
  sumgaussian.rep <- t(sapply(1:nsimu,function(x){comp_momments(rnorm(n,mu[x],sd=sqrt(sigma2[x])))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumgaussian.rep <- cbind(sumgaussian.rep, noise)
  mul <- runif(nsimu,-10,10)
  b2 <- runif(nsimu,qnorm(3/4),4*qnorm(3/4))
  sumcauchy <- t(sapply(1:nsimu,function(x){c(comp_momments(rcauchy(n,mul[x],b2[x])), runif(n.noise))}))
  sumcauchy.rep <- t(sapply(1:nsimu,function(x){c(comp_momments(rcauchy(n,mul[x],b2[x])), runif(n.noise))}))

  target <- rbind(sumgaussian[1, ], sumcauchy[1, ])
  target.rep <- rbind(sumgaussian.rep[1, ], sumcauchy.rep[1, ])
  sumstat <- sumgaussian[-1, ]
  sumstat.rep <- sumgaussian.rep[-1, ]
  param <- cbind(mu[-1], sigma2[-1])

  colnames(target) <- colnames(target.rep) <- colnames(sumstat) <- colnames(sumstat.rep) <- 1:ncol(target)

  ## Test that rejection and freq give the same results
  set.seed(20200430)
  find_oneparam_index <- function(oneparam) which(apply(param, 1, function(x) all(sapply(seq_along(x), function(y) x[y] == oneparam[y]))))
  find_params_index <- function(params) apply(params, 1, find_oneparam_index)
  resgfithp <- hpgfit(target, target.rep, param, sumstat,
                      sim.fun = function(params) sumstat.rep[find_params_index(params),],
                      method = "rejection",
                      eps = 0.1)
  # freq
  set.seed(20200430)
  resgfitfreq <- freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1)
  expect_equal(resgfitfreq$lof$pval, resgfithp$lof$pval)

  # sd
  resgfithp <- compute_pval_sd(resgfithp)
  expect_equal(resgfithp$lof$pval_sd[1, "max"], 0.0459565, ignore_attr = TRUE, tolerance = 1e-4)
  # sd pow
  resgfithp <- compute_power(resgfithp, level = 0.01)
  resgfithp <- compute_power(resgfithp, level = 0.05)
  expect_message(resgfithp <- compute_power_sd(resgfithp))
  expect_equal(resgfithp$lof$`sd_pow_5%`["max"], 0.3535534, ignore_attr = TRUE, tolerance = 1e-5)
  expect_equal(resgfithp$lof$`sd_pow_1%`["max"], 0.3535534, ignore_attr = TRUE, tolerance = 1e-5)

  # resimulate
  set.seed(20200430)
  sim.rnorm <- function(oneparam) c(comp_momments(rnorm(n, oneparam[1], sqrt(oneparam[2]))), runif(n.noise))
  sim.fun.all <- function(params) t(apply(params, 1, sim.rnorm))
  resgfithp <- hpgfit(target, target.rep, param, sumstat,
                      sim.fun = sim.fun.all,
                      method = "rejection",
                      eps = 0.1)
  expect_equal(resgfithp$lof$pval[, "max"] <= 0.05, resgfitfreq$lof$pval[, "max"] <= 0.05)

  # resimulate - locallinear
  set.seed(20200430)
  resgfithploclin <- hpgfit(target, target.rep, param, sumstat,
                            sim.fun = sim.fun.all,
                            method = "loclinear",
                            eps = 0.1)
  expect_equal(resgfithp$lof$pval[, "max"] <= 0.05, resgfithploclin$lof$pval[, "max"] <= 0.05)
  expect_equal(resgfithploclin$lof$pval[1, ], c(0.32, 0.08, 0.46, 0.14, 0.2, 0.32, 0.16), ignore_attr = TRUE)

  # one obs
  set.seed(20200430)
  resgfithploclinone <- hpgfit(target[1, ], target.rep[1, ], param, sumstat,
                               sim.fun = sim.fun.all,
                               method = "loclinear",
                               eps = 0.1)
  expect_equal(resgfithploclin$lof$pval[1, "max"], resgfithploclinone$lof$pval[1, "max"])

  # resimulate - locallinear ridge
  set.seed(20200430)
  resgfithpridge <- hpgfit(target, target.rep, param, sumstat,
                           sim.fun = sim.fun.all,
                           method = "ridge",
                           eps = 0.1,
                           lambda = 10)
  expect_equal(resgfithpridge$lof$pval[1, ],
               c(0.32, 0.92, 0.48, 0.58, 0.46, 0.36, 0.40),
               ignore_attr = TRUE)

  # one obs - lambda = 0 is the same as loclin
  set.seed(20200430)
  resgfithpridgeone <- hpgfit(target[1, ], target.rep[1, ], param, sumstat,
                              sim.fun = sim.fun.all,
                              method = "ridge",
                              eps = 0.1,
                              lambda = 0)
  expect_equal(resgfithpridgeone$lof$pval[1, ], resgfithploclinone$lof$pval[1, ])

  ## bootstrap
  set.seed(20200430)
  resgfithp <- hpgfit(target, target.rep, param, sumstat,
                      sim.fun = sim.fun.all,
                      method = "rejection",
                      eps = 0.1, nboot = 10)
  resgfithp
  sumgfit <- summary(resgfithp)
  expect_equal(sumgfit$pval[1, ], c(0.28, 0.14, 0.54), ignore_attr = TRUE, tolerance = 1e-5)
  # one obs
  set.seed(20200430)
  resgfithp <- hpgfit(target[1, ], target.rep[1, ], param, sumstat,
                      sim.fun = sim.fun.all,
                      method = "rejection",
                      eps = 0.1, nboot = 10)
  resgfithp
  sumgfit <- summary(resgfithp)
  expect_equal(sumgfit$pval[1, ], c(0.28, 0.14, 0.54), ignore_attr = TRUE, tolerance = 1e-5)
  # parallel
  set.seed(20200430)
  resgfithp <- hpgfit(target[1, ], target.rep[1, ], param, sumstat,
                      sim.fun = sim.fun.all,
                      method = "rejection",
                      eps = 0.1, nboot = 10, ncores = 2)
  resgfithp
  sumgfit <- summary(resgfithp)
  expect_true(sumgfit$pval[1, 3] >= sumgfit$pval[1, 2])

  ## summary no boot
  set.seed(20200430)
  resgfithp <-  hpgfit(target, target.rep, param, sumstat,
                       sim.fun = sim.fun.all,
                       method = "rejection",
                       eps = 0.1)
  resgfithp
  sumgfit <- summary(resgfithp)
  expect_equal(sumgfit$pval[1, ], c(0.36, 0.2269532, 0.4930468), ignore_attr = TRUE, tolerance = 1e-5)
  # one obs
  set.seed(20200430)
  resgfithp <-  hpgfit(target[1, ], target.rep[1, ], param, sumstat,
                       sim.fun = sim.fun.all,
                       method = "rejection",
                       eps = 0.1)
  resgfithp
  sumgfit <- summary(resgfithp)
  expect_equal(sumgfit$pval[1, ], c(0.36, 0.2269532, 0.4930468), ignore_attr = TRUE, tolerance = 1e-5)

  # # timing k
  # microbenchmark::microbenchmark(hpgfit(target, target.rep, param, sumstat,
  #                                       sim.fun = sim.fun.all,
  #                                       method = "rejection",
  #                                       eps = 0.1,
  #                                       k = c(2, 5, 10, 15, 20)),
  #                                hpgfit(target, target.rep, param, sumstat,
  #                                       sim.fun = sim.fun.all,
  #                                       method = "rejection",
  #                                       eps = 0.1,
  #                                       k = 2:20),
  #                                times = 10)

})

test_that("Test holdout posterior transforms", {

  ## Data
  nsimu <- 1001
  n <- 150
  n.moments <- 2
  n.noise <- 10

  set.seed(20200430)
  mu <- runif(nsimu,-10,10)
  sigma2 <- runif(nsimu,1,4)
  comp_momments <- function(x){c(mean(x),var(x))}
  sumgaussian <- t(sapply(1:nsimu,function(x){comp_momments(rnorm(n,mu[x],sd=sqrt(sigma2[x])))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumgaussian <- cbind(sumgaussian,noise)
  sumgaussian.rep <- t(sapply(1:nsimu,function(x){comp_momments(rnorm(n,mu[x],sd=sqrt(sigma2[x])))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumgaussian.rep <- cbind(sumgaussian.rep, noise)
  mul <- runif(nsimu,-10,10)
  b2 <- runif(nsimu,qnorm(3/4),4*qnorm(3/4))
  sumcauchy <- t(sapply(1:nsimu,function(x){c(comp_momments(rcauchy(n,mul[x],b2[x])), runif(n.noise))}))
  sumcauchy.rep <- t(sapply(1:nsimu,function(x){c(comp_momments(rcauchy(n,mul[x],b2[x])), runif(n.noise))}))

  target <- rbind(sumgaussian[1, ], sumcauchy[1, ])
  target.rep <- rbind(sumgaussian.rep[1, ], sumcauchy.rep[1, ])
  sumstat <- sumgaussian[-1, ]
  sumstat.rep <- sumgaussian.rep[-1, ]
  param <- cbind(mu[-1], sigma2[-1])

  colnames(target) <- colnames(target.rep) <- colnames(sumstat) <- colnames(sumstat.rep) <- 1:ncol(target)
  colnames(param) <- c("mu", "sigma2")

  # transforms
  expect_error(transmu <- logit(param[, 1], a = -5, b = 10), "interval")
  expect_error(transmu <- logit(param[, 1], a = 15, b = 10), "lower")
  expect_error(transmu <- logistic(param[, 1], a = 15, b = 10), "lower")

  transmu <- logit(param[, 1], a = -10, b = 10)
  expect_equal(logistic(transmu, -10, 10), param[, 1])

  transsig <- logit(param[, 2], a = 1, b = 4)
  expect_equal(logistic(transsig, 1, 4), param[, 2])

  ## check parameters
  # default
  tt <- check_param_transform(param, "none", -Inf, Inf)
  dum <- matrix(rnorm(10), ncol = 2)
  expect_true(all(tt$transform(dum) == dum))
  expect_true(all(tt$back_transform(dum) == dum))
  # errors
  expect_error(check_param_transform(param, c("none", "log", "none"), -Inf, Inf), "a single value")
  expect_error(expect_warning(
    check_param_transform(param, c("none", "log"), -Inf, Inf), "unnamed"),
    "lower bound must be 0")
  expect_error(expect_warning(
    check_param_transform(param, c("none", "logit"), -Inf, Inf), "unnamed"),
    "finite values")
  expect_error(expect_warning(
    check_param_transform(param, c("none", "logit"), 10, 8), "unnamed"),
    "smaller than upper")

  # hpgfit
  set.seed(20200430)
  sim.rnorm <- function(oneparam) c(comp_momments(rnorm(n, oneparam[1], sqrt(oneparam[2]))), runif(n.noise))
  sim.fun.all <- function(params) t(apply(params, 1, sim.rnorm))

  expect_warning(expect_warning(expect_warning(
    expect_error(
      resgfithploclin <- hpgfit(target, target.rep, param, sumstat,
                                sim.fun = sim.fun.all,
                                method = "loclinear",
                                param_transform = c("logit", "logit"),
                                param_upper_bound = c(10, 3),
                                param_lower_bound = c(-10, 1),
                                eps = 0.1),
      "interval"),
    "unnamed"), "unnamed"), "unnamed")

  expect_warning(expect_warning(expect_warning(
    resgfithploclin <- hpgfit(target, target.rep, param, sumstat,
                              sim.fun = sim.fun.all,
                              method = "loclinear",
                              param_transform = c("logit", "logit"),
                              param_upper_bound = c(10, 4),
                              param_lower_bound = c(-10, 1),
                              eps = 0.1),
    "unnamed"), "unnamed"), "unnamed")
  expect_equal(resgfithploclin$lof$pval[1, ], c(0.38, 0.10, 0.72, 0.16, 0.18, 0.34, 0.18), ignore_attr = TRUE)
})
