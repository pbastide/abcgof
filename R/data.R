#' Gaussian Laplace Toy dataset
#'
#' A dataset containing two reference tables, simulated using either the
#' a Gaussian or a Laplace model.
#'
#' @details
#' Each dataset contains 1000 lines.
#' Parameters `mu` and `sigma` were drawn uniformly between -10 and 10, and 1 and 4, respectively.
#' For each set of parameter:
#'    * a Gaussian vector of size 350 was simulated with mean `mu` and standard deviation `sigma` using base function `stats::rnorm`.
#'    * a Laplace vector of size 350 was simulated with location `mu` and scale `sigma/sqrt(2)` using function `VGAM::rlaplace`.
#' Then, for each vector, a set of 20 summary statistics was computed, taking the 20 first normalized L-moments using function `lmom::samlmu`.
#'
#'
#' @format Two lists `dataset_gaussian` and `dataset_laplace`, with entries:
#' \describe{
#'   \item{param}{A matrix with 100,000 rows and 2 columns with the parameters used for the simulations}
#'   \item{sim}{A data frame with the simulated summary statistics (350 columns) for each of the 100,000 rows of the param matrix.}
#'   \item{sim.fun}{The simulation function used.}
#' }
#'
"gaussian_laplace"

# ## Plot PCA
# trainall <- rbind(data.ref, data.test.obs)
# res.pca <- FactoMineR::PCA(trainall, graph = FALSE)
# colo <-c(rep(1,n.ref),rep(2,n.test))
# mix <- c(sample(1:n.ref, 2000), (n.ref+1):(n.ref+n.test))
# pdf(file = here("data", "datasim", paste0(results_name, "_pca.pdf")))
# plot(res.pca$ind$coord[mix,c(1,2)],col=colo[mix],pch="*")
# dev.off()
