#' Objective function to minimize certain elements of the local response matrix
#'
#' This function squares all elements that flip sign once the complexes are added to the observation.
#' This can then be used to minimize the values.
#'
#'
#' @param pars Free parameters a_ij to control the influence of the complexes
#' @param fixed Fixed parameters are parameters for which no derivative is calculated
#' @param perturbation_prediction Simulated steady state perturbation data
#' @param r_kept elements kept in the minimization (=elements that describe sequestration connections)
#' @param obs_fun Observation function which returns the communicating species
#' @param p_fun Parameter transformation which returns different conditions
#' @param mypars Parameters of the ODE model
#' @param method method for taking the derivative wrt pars_opt.
#' @param ... not needed except for throwing away an unnecessary input when applying the new mstrust
#'
#'
#' @export
obj_alpha <- function(pars,
                      fixed = NULL,
                      perturbation_prediction,
                      r_kept,
                      obs_fun,
                      p_fun,
                      mypars,
                      method = "simple",
                      ...) {

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
#' I used this function to compare different optimizers. trust() was as fast as optim, but more reliable.
#' All other methods like simulated annealing were not as reliable.
#'
#' @inheritParams obj_alpha
#'
#' @export
#'
#'
obj_alpha_no_derivs <- function(pars_opt,
                       perturbation_prediction,
                       r_kept,
                       obs_fun,
                       p_fun,
                       mypars
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

