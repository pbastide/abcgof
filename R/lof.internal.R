#' @title LOF one new point
#'
#' @description
#' Compute the LOF for the new point with respect to the data set of reference for various values of k.
#' You need to have pre-computed k_dist for points in the reference data set.
#'
#' @param new The vector that represents the new point
#' @param set.ref The ensemble of points in the reference data set  (row by row : \code{set.ref[j,]} represents the j_th point)
#' @param Tab.dist \code{Tab.dist[i,k]} k-distance of particles i
#' @param Tab.id \code{Tab.id[i,1:k]} k-neighborhood of particles i
#' @param k_all specify the values of k
#'
#' @return vector of LOF values for the new point with the respect to the reference data set, for each k in k_all.
#' @keywords internal
#'
lof.internal <- function(new, set.ref, k_all, Tab.dist = NULL, Tab.id = NULL) {

  if (is.null(Tab.dist)) {
    Tab.dbscan <- dbscan::kNN(set.ref, max(k_all))
    Tab.dist <- Tab.dbscan$dist
    Tab.id <- Tab.dbscan$id
    rm(Tab.dbscan)
  }

  result = numeric(length(k_all)) # will contain LOF values
  dist.new.ref = calDIST(new,set.ref) # distances between the new point and particles in ref
  dist.new.ref.sorted = sort(dist.new.ref) # increasing order.

  i = 0
  for (k in k_all){
    i <- i+1

    ### lrd for the new point :
    Bs <- which(dist.new.ref <= dist.new.ref.sorted[k]) # k-Neighborhood of the new points

    #kdist.new.local = function(j) kdist.new(k.dist.ref[k,j], k.dist.ref[k-1,j], dist.new.ref[j])
    kdist.new.local = function(j) kdist.new(Tab.dist[j,k], Tab.dist[j,k-1], dist.new.ref[j])
    new.k.dist = sapply(Bs, FUN = kdist.new.local) # New k-dist for ref points in the k-Neighborhood of the new points

    inv.lrd.new = mean(pmax.int(new.k.dist, dist.new.ref[Bs])) # inverse of local reachable density for the new point.
    ###

    ### sum lrd for the k-Neighborhood of the new point :

    sum.lrd.neighborhood <- 0
    for (b in Bs){
      if (Tab.dist[b,k] < dist.new.ref[b]) # The new point is not in the k_Neighborhoord of point b.
      {
        kNeighborhood.b = Tab.id[b,1:k] # !! which(dist.ref[b,] <= k.dist.ref[k,b])
        new.k.dist = sapply(kNeighborhood.b, FUN = kdist.new.local) # New k-dist for ref points in the k-Neighborhood of the points b

        dists = Tab.dist[b, 1:k]# !! dist.ref[b, kNeighborhood.b], distances between b and points in his k_Neighborhood

        sum.lrd.neighborhood <- sum.lrd.neighborhood + k / sum(pmax.int(new.k.dist, dists)) # add the lrd of point b

      }
      else # The new point is in the k_Neighborhood of point b
      {

        if (k > 1) {
          old.kNeighborhood.b <- Tab.id[b, seq_len(k-1)] # !! which(dist.ref[b,] <= k.dist.ref[k-1,b]) # old (k-1)-neighborhood of point b
          new.k.dist = sapply(old.kNeighborhood.b, FUN = kdist.new.local) # New k-dist for ref points in the old (k-1)-Neighborhood of the points b
          dists = Tab.dist[b, seq_len(k-1)] # !! dist.ref[b, old.kNeighborhood.b] # distances between b and points in his old (k-1)_Neighborhood
        } else {
          dists <- 0
          new.k.dist <- 0
        }

        sum.lrd.neighborhood <- sum.lrd.neighborhood + k / (sum(pmax.int(new.k.dist, dists)) + pmax(dist.new.ref.sorted[k], dist.new.ref[b]))
      }
    }

    ### Ratio :

    result[i] <- inv.lrd.new * sum.lrd.neighborhood / k # LOF_k for the new point

    # with more than k duplicates lrd can become infinity
    # we define them not to be outliers
    # See dbscan::lof
    result[is.nan(result)] <- 1

  }
  result
}

#' @title New k-dist
#'
#' @description
#' Compute the new k-dist for a point in the reference data set. It uses the old k-distances and the distance to the new point.
#'
#' @param old.kdist k-dist computed without the new points.
#' @param old.belowkdist (k-1)-dist computed without the new points.
#' @param dist.new distance to the new point
#'
#' @return new k-dist including the new point.
#' @keywords internal
#'
kdist.new <- function(old.kdist, old.belowkdist, dist.new)
{
  if (old.kdist < dist.new) #the new point is not in the k-neighborhood, nothing changes
  {
    return(old.kdist)
  }
  # the new point is in the k-neighborhood, so k-dist becomes :
  return(max(old.belowkdist, dist.new))
}

#' @title Mean kNN for one new point
#'
#' @description
#' Compute the mean kNN for the new point with respect to the data set of reference for various values of k.
#'
#' @inheritParams lof.internal
#' @param ... further arguments to be passed to the function (ignored).
#'
#' @return vector of mean kNN values for the new point with the respect to the reference data set, for each k in k_all.
#' @keywords internal
#'
meankNN <- function(new, set.ref, k_all, ...) {
  result = numeric(length(k_all)) # will contain mean kNN values
  dist.new.ref = calDIST(new, set.ref) # distances between the new point and particles in ref
  dist.new.ref.sorted = sort(dist.new.ref) # increasing order.
  return(sapply(k_all, function(kk) mean(dist.new.ref.sorted[1:kk])))
}
