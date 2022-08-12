




JaccardBase <- function(clabs1, clabs2, return_dist=FALSE) {
  # Return: Jaccard index or distance.
  # https://github.com/cran/clusteval/blob/master/src/rcpp_comembership_table.cpp
  # https://en.wikipedia.org/wiki/Jaccard_index
  n1 <- length(clabs1)
  stopifnot(n1 == length(clabs2))
  
  mat1 <- outer(clabs1, clabs1, FUN=`==`) & upper.tri(diag(n1))
  mat2 <- outer(clabs2, clabs2, FUN=`==`) & upper.tri(diag(n1))
  
  # The counts of comembership pairs.
  # n_11: the number of comemberships in both partitions
  # n_10: the number of comemberships in clustering 1 but not in clustering 2
  # n_01: the number of comemberships in clustering 2 but not in clustering 1
  # n_00: the number of non-comemberships in both partitions
  n_11 <- sum(mat1 & mat2)
  n_10 <- sum(mat1 & !mat2)
  n_01 <- sum(!mat1 & mat2)
  
  jsc <- n_11 / (n_01 + n_10 + n_11)  
  ifelse(return_dist, 1-jsc, jsc)
}

JaccardMat <- function(clabsmat, clabsref, return_dist=FALSE, do_parallel=F) {
  # second argument is n-length vector.
  # first argument should be (nd, n) matrix.
  # Return: nd-length vector of Jaccard indices or distances.
  stopifnot(is.matrix(clabsmat))
  stopifnot(is.vector(clabsref))
  if (ncol(clabsmat) != length(clabsref)) {
    print("ncol(clabsmat) != length(clabsref).  Attempting to use t(clabsmat).")
    stopifnot(nrow(clabsmat) == length(clabsref))
    clabsmat <- t(clabsmat)
  }
  if (!do_parallel) {
    return( apply(clabsmat, 1, JaccardBase, clabs2=clabsref, return_dist=return_dist) )
  } else {
    library(parallel)
    numCores <- detectCores()
    nr <- nrow(clabsmat)
    jm_res <- mclapply(1:nr, function(i) JaccardBase(clabsmat[i, ], clabs2=clabsref, return_dist=return_dist))
    return( unlist(jm_res) )
  }
}

if (FALSE) {  # for testing purposes
  (tmpmat <- rbind(c(1,1,2,3), c(4,4,5,5), c(1,1,1,7)))
  (tmpref <- c(1,1,2,2))
  JaccardMat(tmpmat, tmpref, T)
}

