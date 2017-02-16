#' Objective function to minimize certain elements of the local response matrix
#'
#' This function squares all elements that flip sign once the complexes are added to the observation.
#' This can then be used to minimize the values.
#'
#' @param pars
#' @param predfun
#' @param pars_opt_0
#' @param r_fun
#' @param which_pars_perturbed
#'
#' @return The parameter values for which the sign-changing elements are minimal.
#' @export
#'
#' @examples
obj_frobenius <- function(pars_opt = pars_opt_0,
                          r_kept = list(rep(1,length(modules)^2)), #dirty way to keep all elements in condition 1
                          pars = pars_0,
                          predfun = (g*xs),
                          which_pars_perturbed = which_pars_perturbed_0,
                          R_fun = R_numDeriv_fun) {

  if(length(attr(predfun, "conditions"))>1) warning("Prediction function has more than one condition.")

  elements_fun  <- function(pars_optimization) local_response_matrix_eq10(R_fun(pars = pars,
                                            predfun = predfun,
                                            which_pars_perturbed = which_pars_perturbed,
                                            pars_opt = pars_optimization)) [[ 1 ]] [ r_kept [[1]] ]


  elements      <- elements_fun(pars_opt)
  elements_jac  <- numDeriv::jacobian(elements_fun, pars_opt)

  value   <- (elements^2) %>% sum()
  grad    <- (2*(t(elements)%*%elements_jac)) %>% as.numeric %>% structure(names = names(pars_opt))
  hessian <- (2*(t(elements_jac) %*% elements_jac)) %>% structure(dimnames = list(names(pars_opt),names(pars_opt)))
  # hessian <- matrix(0, nrow = length(pars_opt), ncol = length(pars_opt)) %>% structure(dimnames = list(names(pars_opt),names(pars_opt)))

  return(list(value = value, gradient = grad, hessian = hessian))
}


#' Objective function to minimize certain elements of the local response matrix
#'
#' This function squares all elements that flip sign once the complexes are added to the observation.
#' This can then be used to minimize the values.
#'
#' @param pars_opt
#' @param perturbation_prediction
#' @param r_kept
#' @param obs_fun
#' @param p_fun
#' @param pars
#'
#' @return The parameter values for which the sign-changing elements are minimal.
#' @export
#'
#' @examples
obj_alpha <- function(pars_opt = pars_opt_0,
                      perturbation_prediction = perturbation_prediction_0,
                      r_kept = r_kept_0,
                      obs_fun = g,
                      p_fun = p_pert,
                      mypars = pars_0,
                      fixed = NULL) {



  elements_fun  <- function(pars_optimization) local_response_matrix_eq10(R_fun(pars_opt = pars_optimization,
                                                                                perturbation_prediction = perturbation_prediction,
                                                                                obs_fun = obs_fun,
                                                                                p_fun = p_fun,
                                                                                pars = mypars))[r_kept]

  elements      <- elements_fun(pars_opt)
  elements_jac  <- numDeriv::jacobian(elements_fun, pars_opt)

  value   <- (elements^2) %>% sum()
  grad    <- (2*(t(elements)%*%elements_jac)) %>% as.numeric %>% structure(names = names(pars_opt))
  hessian <- (2*(t(elements_jac) %*% elements_jac)) %>% structure(dimnames = list(names(pars_opt),names(pars_opt)))

  return(list(value = value, gradient = grad, hessian = hessian))
}
