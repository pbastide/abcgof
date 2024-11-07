test_that("Test LOF computations", {
  set.seed(20200430)
  skip_if_not_installed("parallelDist")

  ## Data
  nsimu <- 1000
  n <- 150
  n.moments <- 2
  n.noise <- 10

  mu <- runif(nsimu,-10,10)
  sigma2 <- runif(nsimu,1,4)
  # comp_momments <- function(x){lmom::samlmu(x, nmom = n.moments)}
  comp_momments <- function(x){c(mean(x),var(x))}
  sumgaussian <- t(sapply(1:nsimu,function(x){comp_momments(rnorm(n,mu[x],sd=sqrt(sigma2[x])))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumgaussian <- cbind(sumgaussian,noise)
  mul <- runif(nsimu,-10,10)
  b2 <- runif(nsimu,qnorm(3/4),4*qnorm(3/4))
  sumcauchy <- t(sapply(1:nsimu,function(x){comp_momments(rcauchy(n,mul[x],b2[x]))}))
  noise <- matrix(runif(nsimu*n.noise),nsimu,n.noise)
  sumcauchy <- cbind(sumcauchy,noise)

  sumstats <- rbind(sumgaussian, sumcauchy)

  ## LOF
  ktest <- c(1, 2, 4, 6, 8, 15, 20, 25)
  lof1 <- t(sapply(ktest, function(k) dbscan::lof(sumstats, minPts = k+1)))

  for (obs_ind in c(sample(1:nsimu, 5), nsimu+sample(1:nsimu, 5))) {
    all_dists <- parallelDist::parDist(sumstats[-obs_ind, ], diag = TRUE, upper = TRUE)
    all_dists <- as.matrix(all_dists)
    diag(all_dists) <- Inf
    ordered_dist <- apply(all_dists, 1, sort, decreasing = FALSE)

    Tab.dbscan <- dbscan::kNN(sumstats[-obs_ind, ], max(ktest))
    lof3 <- lof.internal(sumstats[obs_ind, ], sumstats[-obs_ind, ], k_all = c(ktest),
                     Tab.dbscan$dist, Tab.dbscan$id)
    for (k in 1:length(ktest)) expect_equal(lof1[k, obs_ind], lof3[k])
  }

  ## replicated points
  ktest <- 20
  for (i in 5:110) sumstats[i, ] <- sumstats[1, ]
  lof1 <- dbscan::lof(sumstats, minPts = ktest+1)
  obs_ind <- 1
  Tab.dbscan <- dbscan::kNN(sumstats[-obs_ind, ], ktest)
  lof3 <- lof.internal(sumstats[obs_ind, ], sumstats[-obs_ind, ], k_all = c(ktest), Tab.dbscan$dist, Tab.dbscan$id)
  expect_equal(lof3, lof1[obs_ind])

})
