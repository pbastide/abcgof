# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Distances of a vector to rows of a matrix
#'
#' @param y vector of numeric values
#' @param x matrix of row vectors
#' @return vector of distances between y and each row of x
#' @keywords internal
calDIST <- function(y, x) {
    .Call(`_abcgof_calDIST`, y, x)
}

