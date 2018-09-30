# A collection of convenience functions mainly to augment the magrittr-package

#' Inverse of matrix using singular value decomposition
#'
#' @param M Matrix
#' @param tol not used right now...
#'
#' @export
#'
svdinv <- function(M, tol = .Machine$double.eps*1e2) {
  mysvd <- svd(M)
  return( mysvd$v %*% diag(1/mysvd$d) %*% t(mysvd$u))
}



