#' Objective function to minimize certain elements of the local response matrix
#'
#' This function squares all elements that flip sign once the complexes are added to the observation.
#' This can then be used to minimize the values.
#'
#' @param pars
#' @param predfun
#' @param pars_optimization_0
#' @param r_fun
#' @param which_pars_perturbed
#'
#' @return The parameter values for which the sign-changing elements are minimal.
#' @export
#'
#' @examples
obj_frobenius <- function(pars_optimization = pars_opt_0,
                                      condition = 1,
                                      r_kept = list(rep(1,length(modules)^2)), #dirty way to keep all elements in condition 1
                                      pars = pars_0,
                                      predfun = (g*xs),
                                      which_pars_perturbed = which_pars_perturbed_0,
                                      r_fun = r_numDeriv_fun) {


  elements_fun  <- function(pars_opt) r_fun(pars = pars,
                                            predfun = predfun,
                                            which_pars_perturbed = which_pars_perturbed,
                                            pars_optimization = pars_opt) [[ condition ]] [ r_kept[[condition]] ]


  elements      <- elements_fun(pars_optimization)
  elements_jac  <- numDeriv::jacobian(elements_fun, pars_optimization)

  value   <- (elements^2) %>% sum()
  grad    <- (2*(t(elements)%*%elements_jac)) %>% as.numeric %>% structure(names = names(pars_optimization))
  hessian <- (2*(t(elements_jac) %*% elements_jac)) %>% structure(dimnames = list(names(pars_optimization),names(pars_optimization)))
  # hessian <- matrix(0, nrow = length(pars_optimization), ncol = length(pars_optimization)) %>% structure(dimnames = list(names(pars_optimization),names(pars_optimization)))

  return(list(value = value, gradient = grad, hessian = hessian))
}



