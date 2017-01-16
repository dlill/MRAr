#' Objective function to minimize certain elements of the local response matrix
#'
#' This function squares all elements that flip sign once the complexes are added to the observation.
#' This can then be used to minimize the values.
#'
#' @param pars
#' @param time_max
#' @param predfun
#' @param which_pars_perturbed
#' @param pars_optimization The parameters which are needed to optimize the matrix
#'
#' @return The parameter values for which the sign-changing elements are minimal.
#' @export
#'
#' @examples
obj_frobenius <- function(pars_optimization = pars_opt_0,
                          modelpars = pars_0,
                          time_max = 12000,
                          predfun = (g*x),
                          deparsed_predfun = NULL,
                          which_pars_perturbed = which_pars_perturbed_0,
                          r_fun = r_numDeriv_fun) {

  if(is.null(deparsed_predfun)) {deparsed_predfun <- deparse(substitute(predfun))}

  conditions <- attr(predfun, "conditions")
  if(is.null(conditions)) conditions <- 1;

  # unfortunately I have to abandon the conditions here :( Only one condition allowed from here

  r_sign_changed <- r_sign_changed_fun(pars = modelpars,
                                       time_max = time_max,
                                       predfun = predfun,
                                       deparsed_predfun = deparsed_predfun,
                                       which_pars_perturbed = which_pars_perturbed,
                                       pars_optimization = pars_optimization,
                                       r_fun = r_fun)[[1]]

  elements_fun  <- function(pars_opt) r_fun(pars = modelpars,
                                            time_max = time_max,
                                            predfun = predfun,
                                            deparsed_predfun = deparsed_predfun,
                                            which_pars_perturbed = which_pars_perturbed,
                                            pars_optimization = pars_opt)[[1]][r_sign_changed]

  elements      <- elements_fun(pars_optimization)
  elements_jac  <- numDeriv::jacobian(elements_fun, pars_optimization)

  value   <- (elements^2) %>% sum()
  grad    <- (2*(t(elements)%*%elements_jac)) %>% as.numeric %>% structure(names = names(pars_optimization))
  hessian <- 2*(t(elements_jac) %*% elements_jac) %>% structure(dimnames = list(names(pars_optimization),names(pars_optimization)))
  # hessian <- matrix(0, nrow = length(pars_optimization), ncol = length(pars_optimization)) %>% structure(dimnames = list(names(pars_optimization),names(pars_optimization)))

  return(list(value = value, gradient = grad, hessian = hessian))
}


#' Title
#'
#' @param fit_args
#'
#' @return
#' @export
#'
#' @examples
fit_fun <- function(fit_args) {

}

#' Title
#'
#' @param fit
#' @param fit_args
#'
#' @return
#' @export
#'
#' @examples
frobenius_fit_analysis <- function(fit, fit_args) {
  g <- fit_args$g
  p <- fit_args$p
  pars <- fit_args$pars
  deriv <- fit_args$deriv
  p_ic_list <- fit_args$p_ic_list

  pars_opt <- fit$argument

  r_fun <- if(deriv == "sens") r_sens_fun else r_numDeriv_fun

  lapply(p_ic_list, function(p_ic) {
    myname <- names(p_ic)
  })


}
