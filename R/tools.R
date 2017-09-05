# A collection of convenience functions mainly to augment the magrittr-package

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



#' str with max.level = 0
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
str0 <- function(...) {
  str(..., max.level = 0)
}

#' str with max.level = 1
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
str1 <- function(...) {
  str(..., max.level = 1)
}

#' str with max.level = 2
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
str2 <- function(...) {
  str(..., max.level = 2)
}


#' magrittr-style two sided comparison
#'
#' @param x
#' @param lower_lim
#' @param upper_lim
#'
#' @return
#' @export
#'
#' @examples
is_in_between <- function(x, lower_lim, upper_lim) {
  (x > lower_lim) & (x < upper_lim)
}


#' magrittr-Style Apply an expression that's more complicated
#'
#' @param x
#' @param expr
#'
#' @return
#' @export
#'
#' @examples
apply_expression <- function(x, ...) {
  myenv <- environment()
  eval(as.expression(substitute(...)), envir = myenv)
}

#' magrittr-Style set attributes
#'
#' @param x
#' @param expr
#'
#' @return
#' @export
#'
#' @examples
set_attributes <- function(x, attribute, value) {
  for(i in 1:length(attribute)) {
    attr(x, attribute[i]) <- value[i]
  }
  return(x)
}





#' The runtime of some code
#'
#' This function can be used, when you want to benchmark, but also access the results.
#' @param ... Some code
#'
#' @return The output of the code with an attribute "runtime"
#' @export
#'
#' @examples
runtime <- function( ... ) {
  pt <- proc.time()
  myenv <- environment()
  out <- eval(as.expression(substitute(...)), envir = myenv)
  pt <- proc.time()-pt
  attr(out, "runtime") <- pt
  return(out)
}
