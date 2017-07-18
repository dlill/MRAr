#' Objective function to minimize certain elements of the local response matrix
#'
#' This function squares all elements that flip sign once the complexes are added to the observation.
#' This can then be used to minimize the values.
#'
#' @param pars_opt
#' @param perturbation_prediction
#' @param r_kept
#' @param obs_fun
#' @param mypars
#' @param fixed
#' @param method Method to take numerical derivatives
#' @param p_fun
#'
#' @return The parameter values for which the sign-changing elements are minimal.
#' @export
#'
#' @examples
obj_alpha <- function(pars = pars_opt_0,
                      fixed = NULL,
                      perturbation_prediction = perturbation_prediction_0,
                      r_kept = r_kept_0,
                      obs_fun = g,
                      p_fun = (p_log * p_pert),
                      mypars = pars_0,
                      method = "simple") {

  # implement functionality of "fixed" when it is needed...
  if (!is.null(fixed)) {
    pars <- c(pars, fixed)
  }

  elements_fun  <- function(pars_optimization) local_response_matrix_eq10(R_fun(pars_opt = pars_optimization,
                                                                                perturbation_prediction = perturbation_prediction,
                                                                                obs_fun = obs_fun,
                                                                                p_fun = p_fun,
                                                                                pars = mypars))[r_kept]

  elements      <- elements_fun(pars)
  elements_jac  <- numDeriv::jacobian(elements_fun, pars, method = method)

  value   <- (elements^2) %>% sum()
  grad    <- (2*(t(elements)%*%elements_jac)) %>% as.numeric %>% structure(names = names(pars))
  hessian <- (2*(t(elements_jac) %*% elements_jac)) %>% structure(dimnames = list(names(pars),names(pars)))

  if(!is.null(fixed)) {
    entries <- !(names(grad)%in%names(fixed))
    grad <- grad[entries]
    hessian <- hessian[entries, entries, drop = F]
  }

  return(list(value = value, gradient = grad, hessian = hessian))
}





#' Objective function which doesn't compute derivatives
#'
#' @param pars_opt
#' @param perturbation_prediction
#' @param r_kept
#' @param obs_fun
#' @param p_fun
#' @param mypars
#'
#' @return
#' @export
#'
#' @examples
obj_alpha_no_derivs <- function(pars_opt = pars_opt_0,
                       perturbation_prediction = perturbation_prediction_0,
                       r_kept = r_kept_0,
                       obs_fun = g,
                       p_fun = (p_log * p_pert),
                       mypars = pars_0
) {



  local_response_matrix_eq10(R_fun(pars_opt = pars_opt,
                                   perturbation_prediction = perturbation_prediction,
                                   obs_fun = obs_fun,
                                   p_fun = p_fun,
                                   pars = mypars)) %>%
    extract(r_kept) %>%
    raise_to_power(2) %>%
    sum()

}

