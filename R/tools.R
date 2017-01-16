#' Automatically generaty parameter trafos to quickly test different parameter values as conditions
#'
#' This function can be used, if one wants quickly wants to check the behaviour of the model with different parameters.
#' The named list ic_list contains all the parameter values that one wants to try out.
#' The function then returns a set of parameter transformations p_ic which can be put into g*x*p_ic to have a look at the model prediction.
#'
#' @param ic_list Named list of parameter values to be set.
#' @param pars The original set of parameters to be replaced.
#'
#' @return A "parfn" of all conditions.
#' @export
#'
#' @examples
p_ic_fun <- function(ic_list, pars = pars_0) {

  mydf <- expand.grid(ic_list)

  p_ic_list <- lapply(1:nrow(mydf), function(i) {
    mytrafo <- structure(names(pars), names = names(pars))
    myics <- mydf[i,] %>% as.character %>% structure(names = names(ic_list))
    mytrafo[names(myics)] <- myics
    p <- P(mytrafo, condition = paste0(names(myics), myics, collapse = "_"))
  })

  p_ic <- NULL
  for(i in 1:length(p_ic_list)) p_ic <- p_ic + p_ic_list[[i]]

  return(p_ic)
}

#' Generate a list of parameter trafos.
#'
#' This function is similar to p_ic_fun, but it returns a list of parfns.
#' `p_ic_list` needs to be used like `lapply(p_ic_list, function(p_ic) {...(g*x*p_ic)})`.
#' This allows to do optimization over more than one condition (e.g. replace the dots by `obj`).
#'
#' @param ic_list Named list of parameter values to be set.
#' @param pars The original set of parameters to be replaced.
#'
#' @return
#' @export
#'
#' @examples
p_ic_list_fun <- function(ic_list, pars = pars_0) {

  mydf <- expand.grid(ic_list)

  p_ic_list <- lapply(1:nrow(mydf), function(i) {
    mytrafo <- structure(names(pars), names = names(pars))
    myics <- mydf[i,] %>% as.character %>% structure(names = names(ic_list))
    mytrafo[names(myics)] <- myics
    p <- P(mytrafo, condition = paste0(names(myics), myics, collapse = "_"))
  })

  names(p_ic_list) <- lapply(p_ic_list, function(p_ic) attr(p_ic, "conditions")) %>% do.call(c,.)

  return(p_ic_list)
}




#' The runtime of some code
#'
#' @param code
#'
#' @return The result of the code
#' @export
#'
#' @examples
runtime <- function( code ) {

  pt <- proc.time()
  out <- code
  pt <- proc.time()-pt
  attr(out, "runtime") <- pt
  return(out)
}


