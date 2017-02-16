#' Calculate the local response matrix r from the global responses
#'
#' @param R The matrix of global responses
#'
#' @details
#'
#' @return r The local response matrix, with the attribute which_pars_perturbed
#' @export
#'
#' @examples
local_response_matrix_eq10 <- function(R) {

  R_svd <- svd(R)

  condition_number <- R_svd$d[1]/R_svd$d[length(R_svd$d)]
  if(abs(condition_number) > 1e12) warning("R is ill conditioned. Condition number is ", condition_number)

  R_inv <- R_svd$v %*% diag(1/R_svd$d) %*% t(R_svd$u)

  R_inv_diag <- diag(R_inv)

  r <-  - R_inv / R_inv_diag
  dimnames(r) <- list(rownames(R), rownames(R))
  return(r)
}


#' Compute R from fractional changes of steady states
#'
#' @param pars_opt
#' @param perturbation_prediction
#' @param obs_fun
#' @param p_fun
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
R_fun <- function(pars_opt = pars_opt_0,
                  perturbation_prediction = perturbation_prediction_0,
                  obs_fun = g,
                  p_fun = p_pert,
                  pars = pars_0) {

  mypars <- pars
  mypars[names(pars_opt)] <- pars_opt
  mypars <- p_fun(mypars)
  # Apply the observation function to add the states

  my_perturbation_prediction <- do.call(c, lapply(1:length(perturbation_prediction), function(i) obs_fun(out = perturbation_prediction[[i]], pars = mypars[[i]])))
  names(my_perturbation_prediction) <- attr(p_fun, "conditions")


  # R
  sapply(2:length(my_perturbation_prediction), function(i) {
    Ctr <- my_perturbation_prediction[[1]]
    perturbed <- my_perturbation_prediction[[i]]
    return((2*(perturbed-Ctr)/(Ctr+perturbed))[,-1])
  })

}


#' Which elements to keep in the optimization procedure?
#'
#' Elements are kept if they are either small or if they change sign when switching on the complexes.
#'
#' To do: have a look at what happens, if pars_opt = .5 or combinations of pars_opt (eg c(0,1), c(1,0) ...)
#'
#' @param pars_opt
#' @param pars
#' @param time_max
#' @param predfun
#' @param which_pars_perturbed
#'
#' @return
#' @export
#'
#' @examples
r_kept_fun <- function(pars_opt = pars_opt_0,
                       perturbation_prediction = perturbation_prediction_0,
                       obs_fun = g,
                       p_fun = p_pert,
                       pars = pars_0) {

  # set Pars_opt to zero
  pars_opt <- structure(rep(0, length(pars_opt)), names = names(pars_opt))

  r_0 <- R_fun(pars_opt = pars_opt,
               perturbation_prediction = perturbation_prediction,
               obs_fun = obs_fun,
               p_fun = p_fun,
               pars = pars) %>% local_response_matrix_eq10()

  r_1 <- R_fun(pars_opt = pars_opt+1,
               perturbation_prediction = perturbation_prediction,
               obs_fun = obs_fun,
               p_fun = p_fun,
               pars = pars) %>% local_response_matrix_eq10()

  # Which elements are kept in the optimization process?

    signs <- (sign(r_0) - sign(r_1)) %>% as.logical() %>% matrix(nrow = nrow(r_0))
    abs_vals <- (abs(r_0) < 1e-3) | (abs(r_1) < 1e-3)

    keep <- (abs_vals | signs)

  return(keep)
}
