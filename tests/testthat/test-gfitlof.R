test_that("Test Prior GoF LOF computations", {

  ## Data
  nsimu <- 1000
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
  mul <- runif(nsimu,-10,10)
  b2 <- runif(nsimu,qnorm(3/4),4*qnorm(3/4))
  sumcauchy <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumcauchy <- cbind(sumcauchy,noise)

  target <- sumgaussian[1:100, ]
  sumstat <- sumcauchy
  colnames(target) <- colnames(sumstat) <- 1:ncol(target)

  # gfitlof
  set.seed(20200430)
  resgfit <- gfit(target, sumstat)

  # normalize
  set.seed(20200430)
  calibid <- sample(seq_len(nrow(sumstat)), size = 100)
  normdata <- norm_stats(target, sumstat, sd)
  data.test <- normdata$target
  data.ref <- normdata$sumstat[-calibid, ]
  data.calib <- normdata$sumstat[calibid, ]

  # lof gfit
  resfit <- gof.fit(data.test, data.ref, data.calib, score = "lof")
  # equal
  expect_equal(resgfit$lof$pval, resfit$lof$pval)

  # one obs
  set.seed(20200430)
  resgfitone <- gfit(target[1, ], sumstat)
  expect_equal(resgfitone$lof$pval[1], resgfit$lof$pval[1])

  # manual version
  Tab.dbscan = dbscan::kNN(data.ref,20) #Pre-computation
  dist = Tab.dbscan$dist
  id = Tab.dbscan$id
  k_all <- c(2, 3, 5, 8, 14, 20)
  n.calib = dim(data.calib)[1]
  Tab.diagnostic.calib = matrix(0,n.calib,length(k_all))
  for(i in 1:n.calib){
    pods = data.calib[i,]
    Tab.diagnostic.calib[i,] = lof.internal(pods,data.ref,k_all,dist,id)
  }
  n.test = dim(data.test)[1]
  pval.test = matrix(0,n.test,length(k_all))
  for(i in 1:n.test){
    pods = data.test[i,]
    for(j in 1:length(k_all)){
      diagnostic.test = lof.internal(pods,data.ref,k_all[j],dist,id)
      pval.test[i,j] = mean(Tab.diagnostic.calib[,j]>diagnostic.test)
    }
  }
  # Power
  power.test = matrix(0,length(k_all))
  for(j in 1:length(k_all)){
    power.test[j] = mean(pval.test[,j]<0.05)
  }

  # lof gfit pre-computed
  resfitpre <- gof.fit(data.test, data.ref, data.calib, k_all, range(k_all), "lof", 1, dist, id)

  # equal ?
  expect_equal(pval.test, resfit$lof$pval[, -ncol(resfit$lof$pval)], ignore_attr = TRUE)
  expect_equal(resfitpre, resfit)

  # power
  resgfit <- compute_power(resgfit, 0.05)
  expect_equal(resgfit$lof$`pow_5%`[-length(resgfit$lof$`pow_5%`)], power.test, ignore_attr = TRUE)

  # sd
  resgfit <- compute_pval_sd(resgfit)
  expect_equal(resgfit$lof$pval_sd[1, "max"], 0.047696960, ignore_attr = TRUE)

  # sd pow
  resgfit <- compute_power(resgfit, level = 0.01)
  expect_message(resgfit <- compute_power_sd(resgfit))
  expect_equal(resgfit$lof$`sd_pow_5%`["max"], 0.02179449, ignore_attr = TRUE, tolerance = 1e-5)
  expect_equal(resgfit$lof$`sd_pow_1%`["max"], 0.0, ignore_attr = TRUE, tolerance = 1e-5)

  ## kNN manual
  n.calib = dim(data.calib)[1]
  Tab.diagnostic.calib = matrix(0,n.calib,length(k_all))
  for(i in 1:n.calib){
    pods = data.calib[i,]
    sort.dist.pods.ref = sort(calDIST(pods,data.ref))
    for(j in 1:length(k_all)){
      Tab.diagnostic.calib[i,j] = mean(sort.dist.pods.ref[1:k_all[j]])
    }
  }
  n.test = dim(data.test)[1]
  pval.blum.prior = matrix(0,n.test,length(k_all))
  for(i in 1:n.test){
    pods = data.test[i,]
    sort.dist.pods.ref = sort(calDIST(pods,data.ref))
    for(j in 1:length(k_all)){
      diagnostic.test = mean(sort.dist.pods.ref[1:k_all[j]])
      pval.blum.prior[i,j] = mean(Tab.diagnostic.calib[,j]>diagnostic.test)
    }
  }

  expect_equal(pval.blum.prior, resgfit$kNN$pval[, -ncol(resgfit$kNN$pval)], ignore_attr = TRUE)

  # gfitlof boot
  set.seed(20200430)
  resgfit <- gfit(target[1:2, ], sumstat, k = 1:20, k_range = c(5, 20), nboot = 10)
  resgfit
  sumgfit <- summary(resgfit)
  expect_equal(sumgfit$pval[1, ], c(0.545, 0.33, 0.75), ignore_attr = TRUE)
  sumgfit <- summary(resgfit, score = "kNN", k = 1)
  expect_equal(sumgfit$pval[1, ], c(0.34, 0.17, 0.46), ignore_attr = TRUE, tolerance = 1e-5)
  # one obs
  set.seed(20200430)
  resgfit <- gfit(target[1, ], sumstat, k = 1:20, k_range = c(5, 20), nboot = 10)
  resgfit
  sumgfit <- summary(resgfit)
  expect_equal(sumgfit$pval[1, ], c(0.545, 0.33, 0.75), ignore_attr = TRUE)
  sumgfit <- summary(resgfit, score = "kNN", k = 1)
  expect_equal(sumgfit$pval[1, ], c(0.34, 0.17, 0.46), ignore_attr = TRUE, tolerance = 1e-5)
  # parallel
  set.seed(20200430)
  resgfit <- gfit(target[1, ], sumstat, k = 1:20, k_range = c(5, 20), nboot = 10, ncores = 2)
  sumgfit <- summary(resgfit)
  expect_true(sumgfit$pval[1, 3] >= sumgfit$pval[1, 2])

  # gfitlof summary no boot
  set.seed(20200430)
  resgfit <- gfit(target[1:2, ], sumstat, k = 1:20, k_range = c(5, 20))
  resgfit
  sumgfit <- summary(resgfit)
  expect_equal(sumgfit$pval[1, ], c(0.33000, 0.23784, 0.42216), ignore_attr = TRUE, tolerance = 1e-5)
  sumgfit <- summary(resgfit, score = "kNN", k = 1)
  expect_equal(sumgfit$pval[1, ], c(0.2800000, 0.1919978, 0.3680022), ignore_attr = TRUE, tolerance = 1e-5)

  set.seed(20200430)
  resgfit <- gfit(target[1, ], sumstat, k = 1:20, k_range = c(5, 20))
  resgfit
  sumgfit <- summary(resgfit)
  expect_equal(sumgfit$pval[1, ], c(0.33000, 0.23784, 0.42216), ignore_attr = TRUE, tolerance = 1e-5)
  sumgfit <- summary(resgfit, score = "kNN", k = 1)
  expect_equal(sumgfit$pval[1, ], c(0.2800000, 0.1919978, 0.3680022), ignore_attr = TRUE, tolerance = 1e-5)

  # benchmarck
  # microbenchmark::microbenchmark(gof.fit(data.test, data.ref, data.calib, k_all, "lof", dist, id),
  #                                {  n.calib = dim(data.calib)[1]
  #                                Tab.diagnostic.calib = matrix(0,n.calib,length(k_all))
  #                                for(i in 1:n.calib){
  #                                  pods = data.calib[i,]
  #                                  Tab.diagnostic.calib[i,] = lof.internal(pods,data.ref,k_all,dist,id)
  #                                }
  #                                n.test = dim(data.test)[1]
  #                                pval.test = matrix(0,n.test,length(k_all))
  #                                for(i in 1:n.test){
  #                                  pods = data.test[i,]
  #                                  for(j in 1:length(k_all)){
  #                                    diagnostic.test = lof.internal(pods,data.ref,k_all[j],dist,id)
  #                                    pval.test[i,j] = mean(Tab.diagnostic.calib[,j]>diagnostic.test)
  #                                  }
  #                                }}, times = 10)


})

test_that("Test Posterior GoF LOF computations", {

  ## Data
  nsimu <- 1000
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
  mul <- runif(nsimu,-10,10)
  b2 <- runif(nsimu,qnorm(3/4),4*qnorm(3/4))
  sumcauchy <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumcauchy <- cbind(sumcauchy,noise)
  sumcauchy.rep <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise), nsimu, n.noise)
  sumcauchy.rep <- cbind(sumcauchy.rep, noise)

  target <- sumgaussian[1:100, ]
  sumstat <- sumcauchy
  sumstat.rep <- sumcauchy.rep
  colnames(target) <- colnames(sumstat) <- colnames(sumstat.rep) <- 1:ncol(target)

  # gfitlof
  set.seed(20200430)
  resgfit <- postgfit(target, sumstat, sumstat.rep, eps = 0.1)

  # sd
  resgfit <- compute_pval_sd(resgfit)
  expect_equal(resgfit$lof$pval_sd[1, "max"], 0.009949874, ignore_attr = TRUE, tolerance = 1e-4)

  # normalize
  set.seed(20200430)
  calibid <- sample(seq_len(nrow(sumstat)), size = 100)
  normdata <- norm_stats(target, sumstat, sd)
  data.test <- normdata$target
  data.ref <- normdata$sumstat[-calibid, ]
  data.calib <- normdata$sumstat[calibid, ]
  normdata.rep <- norm_stats(sumstat.rep, sumstat, sd)
  data.ref.replica <- normdata.rep$target[-calibid, ]

  # lof gfit
  resfit <- gof.fit.post(data.test, data.calib, data.ref, data.ref.replica, eps = 0.1, score = "lof")
  # equal
  expect_equal(resgfit$lof$pval, resfit$lof$pval)

  # one obs
  set.seed(20200430)
  resgfitone <- postgfit(target[1, ], sumstat, sumstat.rep, eps = 0.1)
  expect_equal(resgfitone$lof$pval[1], resgfit$lof$pval[1])

  # manual version
  eps <- 0.1
  k_all <- c(2, 3, 5, 8, 14, 20)
  # Calibration
  n.calib = dim(data.calib)[1]
  Tab.diagnostic.calib = matrix(0,n.calib,length(k_all))
  for(i in 1:n.calib){
    pods = data.calib[i,]
    #Posterior
    N = round(eps * dim(data.ref)[1]) # number of particles in the "posterior"
    dists.pods.ref = calDIST(pods,data.ref) # distances from pods
    index.best = order(dists.pods.ref)[1:N] # index of posterior particles
    Posterior = data.ref.replica[index.best,]

    Tab.dbscan.post = dbscan::kNN(Posterior,20) #Pre-computation

    Tab.diagnostic.calib[i,] = lof.internal(pods,Posterior,k_all,Tab.dbscan.post$dist,Tab.dbscan.post$id)
  }
  # Test
  n.test = dim(data.test)[1]
  pval.test = matrix(0,n.test,length(k_all))
  for(i in 1:n.test){
    pods = data.test[i,]
    #Posterior
    N = round(eps * dim(data.ref)[1]) # number of particles in the "posterior"
    dists.pods.ref = abcgof:::calDIST(pods,data.ref) # distances from pods
    index.best = order(dists.pods.ref)[1:N] # index of posterior particles
    Posterior = data.ref.replica[index.best,]

    Tab.dbscan.post = dbscan::kNN(Posterior,20) #Pre-computation

    diagnostic.test = lof.internal(pods,Posterior,k_all,Tab.dbscan.post$dist,Tab.dbscan.post$id)
    for(j in 1:length(k_all)){
      pval.test[i,j] = mean(Tab.diagnostic.calib[,j]>diagnostic.test[j])
    }
  }
  # Power
  power.test = matrix(0,length(k_all))
  for(j in 1:length(k_all)){
    power.test[j] = mean(pval.test[,j]<0.05)
  }

  # equal ?
  expect_equal(pval.test, resfit$lof$pval[, -ncol(resfit$lof$pval)], ignore_attr = TRUE)

  # power
  resgfit <- compute_power(resgfit, 0.05)
  expect_equal(resgfit$lof$`pow_5%`[-length(resgfit$lof$`pow_5%`)], power.test, ignore_attr = TRUE)

  # # benchmarck
  # microbenchmark::microbenchmark(lof.gfit.post(data.test, data.ref, data.ref.replica, eps = 0.1),
  #                                {
  #                                  # manual version
  #                                  eps <- 0.1
  #                                  # Calibration
  #                                  n.calib = dim(data.calib)[1]
  #                                  Tab.diagnostic.calib = matrix(0,n.calib,length(k_all))
  #                                  for(i in 1:n.calib){
  #                                    pods = data.calib[i,]
  #                                    #Posterior
  #                                    N = round(eps * dim(data.ref)[1]) # number of particles in the "posterior"
  #                                    dists.pods.ref = abcgof:::calDIST(pods,data.ref) # distances from pods
  #                                    index.best = order(dists.pods.ref)[1:N] # index of posterior particles
  #                                    Posterior = data.ref.replica[index.best,]
  #
  #                                    Tab.dbscan.post = dbscan::kNN(Posterior,20) #Pre-computation
  #
  #                                    Tab.diagnostic.calib[i,] = lof.internal(pods,Posterior,k_all,Tab.dbscan.post$dist,Tab.dbscan.post$id)
  #                                  }
  #                                  # Test
  #                                  n.test = dim(data.test)[1]
  #                                  pval.test = matrix(0,n.test,length(k_all))
  #                                  for(i in 1:n.test){
  #                                    pods = data.test[i,]
  #                                    #Posterior
  #                                    N = round(eps * dim(data.ref)[1]) # number of particles in the "posterior"
  #                                    dists.pods.ref = abcgof:::calDIST(pods,data.ref) # distances from pods
  #                                    index.best = order(dists.pods.ref)[1:N] # index of posterior particles
  #                                    Posterior = data.ref.replica[index.best,]
  #
  #                                    Tab.dbscan.post = dbscan::kNN(Posterior,20) #Pre-computation
  #
  #                                    diagnostic.test = lof.internal(pods,Posterior,k_all,Tab.dbscan.post$dist,Tab.dbscan.post$id)
  #                                    for(j in 1:length(k_all)){
  #                                      pval.test[i,j] = mean(Tab.diagnostic.calib[,j]>diagnostic.test[j])
  #                                    }
  #                                  }
  #                                }, times = 10)


})

test_that("Test Frequentist GoF LOF computations", {

  ## Data
  nsimu <- 1000
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
  sumcauchy <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumcauchy <- cbind(sumcauchy,noise)
  sumcauchy.rep <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise), nsimu, n.noise)
  sumcauchy.rep <- cbind(sumcauchy.rep, noise)

  target <- sumgaussian[1:10, ]
  target.rep <- sumgaussian.rep[1:10, ]
  sumstat <- sumcauchy
  sumstat.rep <- sumcauchy.rep

  colnames(target) <- colnames(target.rep) <- colnames(sumstat) <- colnames(sumstat.rep) <- 1:ncol(target)

  # gfitlof
  set.seed(20200430)
  resgfit <- freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1)

  # normalize
  normdata <- norm_stats(rbind(target, target.rep), sumstat, sd)
  data.test.obs <- normdata$target[1:10, ]
  data.test.new <- normdata$target[-(1:10), ]
  data.ref <- normdata$sumstat
  normdata.rep <- norm_stats(sumstat.rep, sumstat, sd)
  data.ref.replica <- normdata.rep$target

  # lof gfit
  set.seed(20200430)
  resfit <- gof.fit.freq(data.test.obs, data.test.new, data.ref, data.ref.replica, eps = 0.1, split = 0.5, score = "lof")
  # equal
  expect_equal(resgfit$lof, resfit$lof)

  # sd
  resgfit <- compute_pval_sd(resgfit)
  expect_equal(resgfit$lof$pval_sd[1, "max"], 0.1097269, ignore_attr = TRUE, tolerance = 1e-4)

  # one obs
  set.seed(20200430)
  resgfitone <- freqgfit(target[1, ], target.rep[1, ], sumstat, sumstat.rep, eps = 0.1)
  expect_equal(resgfitone$lof$pval[1], resgfit$lof$pval[1])

  # manual version
  eps <- 0.1
  k_all <- c(2, 3, 5, 8, 14, 20)
  N = round(eps * dim(data.ref)[1]) # number of particles in the "posterior"
  split <- 0.5
  n.calib = round(N*split)
  n.test <- nrow(data.test.obs)
  set.seed(20200430)
  pvals <- foreach::foreach (i = 1:n.test) %do% {
    # init tables
    Tab.lof.calib = matrix(0,n.calib,length(k_all))
    Tab.kNN.calib = matrix(0,n.calib,length(k_all))
    Tab.lofmax.calib = matrix(0,n.calib)
    # init pvals
    pval.test.lofi = rep(0, length(k_all))
    pval.test.lofmaxi = 0
    pval.test.kNNi = rep(0, length(k_all))
    pods = data.test.obs[i,]
    pnew = data.test.new[i,]
    dists.pods.ref = calDIST(pods,data.ref) # distances from pods
    index.best = order(dists.pods.ref)[1:N] # index of posterior particles
    Posterior = data.ref.replica[index.best,]
    index.split = sample(N,n.calib)
    Posterior.split.calib = Posterior[index.split,]
    Posterior.split.cloud = Posterior[-index.split,]
    Tab.dbscan.post = dbscan::kNN(Posterior.split.cloud,20) #Pre-computation
    for(c in 1:n.calib){
      Tab.lof.calib[c,] = lof.internal(Posterior.split.calib[c,],Posterior.split.cloud,k_all,Tab.dbscan.post$dist,Tab.dbscan.post$id)
      Tab.lofmax.calib[c] = max(Tab.lof.calib[c,])

      sort.dist.pcalib.ref = sort(calDIST(Posterior.split.calib[c,],Posterior.split.cloud))
      for(j in 1:length(k_all)){
        Tab.kNN.calib[c,j] = mean(sort.dist.pcalib.ref[1:k_all[j]])
      }
    }
    Tab.lof.test = lof.internal(pnew,Posterior.split.cloud,k_all,Tab.dbscan.post$dist,Tab.dbscan.post$id)
    Tab.kNN.test = matrix(0,length(k_all))
    sort.dist.ptestnew.ref = sort(abcgof:::calDIST(pnew,Posterior.split.cloud))
    for(j in 1:length(k_all)){
      Tab.kNN.test[j] = mean(sort.dist.ptestnew.ref[1:k_all[j]])
    }
    for(j in 1:length(k_all)){
      pval.test.lofi[j] = mean(Tab.lof.calib[,j]>Tab.lof.test[j])
    }
    for(j in 1:length(k_all)){
      pval.test.kNNi[j] = mean(Tab.kNN.calib[,j]>Tab.kNN.test[j])
    }
    pval.test.lofmaxi = mean(Tab.lofmax.calib>max(Tab.lof.test))
    return(list(pval.test.lof = pval.test.lofi,
                pval.test.lofmax = pval.test.lofmaxi,
                pval.test.kNN = pval.test.kNNi))
  }
  ## format results
  pval.test.lof <- do.call(rbind, lapply(pvals, function(pp) pp$pval.test.lof))
  pval.test.lofmax <- sapply(pvals, function(pp) pp$pval.test.lofmax)
  pval.test.kNN <- do.call(rbind, lapply(pvals, function(pp) pp$pval.test.kNN))
  # power
  power.lof.freq = matrix(0,length(k_all))
  for(j in 1:length(k_all)){
    power.lof.freq[j] = mean(pval.test.lof[,j]<0.05)
  }

  ## Equal ?
  expect_equal(resgfit$lof$pval, cbind(pval.test.lof, pval.test.lofmax), ignore_attr = TRUE)
  expect_equal(resgfit$kNN$pval[, -6], pval.test.kNN, ignore_attr = TRUE)

  ## power
  resgfit <- compute_power(resgfit, 0.05)
  expect_equal(resgfit$lof$`pow_5%`[-length(resgfit$lof$`pow_5%`)], power.lof.freq, ignore_attr = TRUE)

})

test_that("Test errors and warnings", {

  ## Data
  nsimu <- 1000
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
  sumcauchy <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumcauchy <- cbind(sumcauchy,noise)
  sumcauchy.rep <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise), nsimu, n.noise)
  sumcauchy.rep <- cbind(sumcauchy.rep, noise)

  target <- sumgaussian[1:2, ]
  target.rep <- sumgaussian.rep[1:2, ]
  sumstat <- sumcauchy
  sumstat.rep <- sumcauchy.rep

  resgfit <- expect_warning(gfit(target, sumstat), "unnamed")
  resgfit <- expect_warning(expect_warning(postgfit(target, sumstat, sumstat.rep, eps = 0.1), "unnamed"), "unnamed")
  resgfit <- expect_warning(expect_warning(expect_warning(freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1), "unnamed"), "unnamed"), "unnamed")

  colnames(target) <- colnames(target.rep) <- colnames(sumstat) <- colnames(sumstat.rep) <- LETTERS[1:ncol(target)]
  colnames(target) <- colnames(target)[c(2, 1, 3:ncol(target))]
  resgfit <- expect_warning(gfit(target, sumstat), "order")
  resgfit <- expect_warning(expect_warning(postgfit(target, sumstat, sumstat.rep, eps = 0.1), "order"), "order")
  resgfit <- expect_warning(expect_warning(expect_warning(freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1), "order"), "order"), "order")

  colnames(target) <- c(colnames(target)[-1], "toto")
  expect_error(gfit(target, sumstat), "names")
  expect_error(postgfit(target, sumstat, sumstat.rep, eps = 0.1), "names")
  expect_error(freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1), "names")

  target <- target[, -1]
  expect_error(gfit(target, sumstat), "same number of columns")
  expect_error(postgfit(target, sumstat, sumstat.rep, eps = 0.1), "same number of columns")
  expect_error(freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1), "same number of columns")

  target <- sumgaussian[1:3, ]
  colnames(target) <- LETTERS[1:ncol(target)]
  expect_error(freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1), "same number of lines")

  target <- target[1:2, ]
  sumstat.rep <- sumstat.rep[-1, ]
  expect_error(postgfit(target, sumstat, sumstat.rep, eps = 0.1), "same number of lines")
  expect_error(freqgfit(target, target.rep, sumstat, sumstat.rep, eps = 0.1), "same number of lines")

})
