
#' Inverse of matrix using singular value decomposition
#'
#' @param M Matrix
#' @param tol not used right now...
#'
#' @return
#' @export
#'
#' @examples
svdinv <- function(M, tol = .Machine$double.eps*1e2) {
  mysvd <- svd(M)
  return( mysvd$v %*% diag(1/mysvd$d) %*% t(mysvd$u))
}
