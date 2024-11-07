## Code to simulate a toy "Laplace - Gaussian" dataset

library(VGAM)
library(lmom)
library(doParallel)
library(foreach)

n.ref <- 1000 # number of particles in each dataset
set.seed(1289)

## Simulate Gaussian dataset
simgaussian <- function(n) {
  mu <- runif(n, -10, 10)
  sigma <- runif(n, 1, 4)
  param <- cbind(mu, sigma)
  colnames(param) <- c("mu", "sigma")
  sim.fun <- function(params) {
    t(apply(params, 1, function(oneparam) lmom::samlmu(stats::rnorm(350, mean = oneparam[1], sd = oneparam[2]), nmom = 20)))
  }
  sim <- sim.fun(param)
  simrep <- sim.fun(param)
  return(list(param = param, sim = sim, sim.fun = sim.fun))
}
dataset_gaussian <- simgaussian(n.ref)

## Simulate Laplace dataset
simlaplace <- function(n) {
  mu <- runif(n, -10, 10)
  sigma <- runif(n, 1, 4)
  param <- cbind(mu, sigma/sqrt(2))
  colnames(param) <- c("mu", "scale")
  sim.fun <- function(params, ...) {
    t(apply(params, 1, function(oneparam) lmom::samlmu(VGAM::rlaplace(350, location = oneparam[1], scale = oneparam[2]), nmom = 20)))
  }
  sim <- sim.fun(param)
  simrep <- sim.fun(param)
  return(list(param = param, sim = sim, sim.fun = sim.fun))
}
dataset_laplace <- simlaplace(n.ref)

## Dataset
gaussian_laplace <- list(dataset_gaussian = dataset_gaussian,
                         dataset_laplace = dataset_laplace)

usethis::use_data(gaussian_laplace, overwrite = TRUE)
