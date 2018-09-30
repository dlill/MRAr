#' Calculate the local response matrix r from the global responses
#'
#' @param R The (square) matrix of global responses, in which column j means that module j only was perturbed
#'
#' @details
#'
#' @return r The local response matrix, with the attribute which_pars_perturbed
#' @export
#'
#' @examples
local_response_matrix_eq10 <- function(R) {

  R_svd <- svd(R)

  condition_number <- R_svd$d[length(R_svd$d)]/R_svd$d[1]
  if(abs(condition_number) < .Machine$double.eps) warning("R is ill conditioned. Condition number is ", condition_number)

  R_inv <- R_svd$v %*% diag(1/R_svd$d) %*% t(R_svd$u)

  R_inv_diag <- diag(R_inv)

  r <-  - R_inv / R_inv_diag
  dimnames(r) <- list(rownames(R), rownames(R))
  return(r)
}


#' Compute R from fractional changes of steady states
#'
#' @param pars_opt Free parameters a_ij to control the influence of the complexes
#' @param perturbation_prediction Simulated steady state perturbation data
#' @param obs_fun Observation function which returns the communicating species
#' @param p_fun Parameter transformation which returns different conditions
#' @param pars Parameters of the ODE model
#'
#' @return
#' @export
#'
#' @examples
R_fun <- function(pars_opt,
                  perturbation_prediction,
                  obs_fun,
                  p_fun,
                  pars ) {
  mypars <- pars
  mypars[names(pars_opt)] <- pars_opt
  mypars <- p_fun(mypars, deriv = F)

  my_perturbation_prediction <- lapply(1:length(perturbation_prediction),
                                       function(i) obs_fun(out = perturbation_prediction[[i]],
                                                           pars = mypars[[i]],
                                                           deriv = F)) %>% do.call(c,.)
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
#' This function compares the connection coefficient in the local response matrices when either only free
#' active kinases are measured or when the complexes are switched on.
#' All free parameters a_ij are increased simultaneously
#'
#' @param pars_opt Free parameters a_ij to control the influence of the complexes
#' @param perturbation_prediction Simulated steady state perturbation data
#' @param obs_fun Observation function which returns the communicating species
#' @param p_fun Parameter transformation which returns different conditions
#' @param pars Parameters of the ODE model
#' @param alpha Compares r(pars_opt) with r(pars_opt + alpha).
#' Since pars_opt_0 defaults to log(0+.Machine$double.eps), I chose its negative value as the default.
#' This is basically the same as comparing r(0) with r(1)
#'
#' @return
#' @export
#'
#' @examples
r_kept_fun <- function(pars_opt,
                       perturbation_prediction,
                       obs_fun,
                       p_fun,
                       pars,
                       alpha = -log(0+.Machine$double.eps)) {

  r_0 <- R_fun(pars_opt = pars_opt,
               perturbation_prediction = perturbation_prediction,
               obs_fun = obs_fun,
               p_fun = p_fun,
               pars = pars) %>% local_response_matrix_eq10()

  r_1 <- R_fun(pars_opt = pars_opt + alpha, #yields asymptotically the same results as only complexes without free species.
               perturbation_prediction = perturbation_prediction,
               obs_fun = obs_fun,
               p_fun = p_fun,
               pars = pars) %>% local_response_matrix_eq10()

  # Which elements are kept in the optimization process?

  # sign switch
  sign_flips <- (sign(r_0) - sign(r_1)) %>% as.logical() %>% matrix(nrow = nrow(r_0))

  return(sign_flips)
}


#' Function to compare two local response matrices
#' in order to choose which elements to keep in the minimization process
#'
#' @param r_0
#' @param r_alpha
#'
#' @return
#' @export
#'
#' @examples
r_kept_fun2 <-function(r_0, r_alpha) {

  # Which elements are kept in the optimization process?
  sign_flips <- (sign(r_0) - sign(r_alpha)) %>% as.logical() %>% matrix(nrow = nrow(r_0))

  return(sign_flips)
}


#' Parametric local response matrix
#'
#' @param pars_opt Free parameters a_ij to control the influence of the complexes
#' @param perturbation_prediction Simulated steady state perturbation data
#' @param obs_fun Observation function which returns the communicating species
#' @param p_fun Parameter transformation which returns different conditions
#' @param pars Parameters of the ODE model
#'
#' @return
#' @export
#'
#' @examples
r_alpha_fun <- function(pars_opt,
                    perturbation_prediction,
                    obs_fun,
                    p_fun,
                    pars) {

  R_fun(pars_opt = pars_opt,
        perturbation_prediction = perturbation_prediction,
        obs_fun = obs_fun,
        p_fun = p_fun,
        pars = pars) %>% local_response_matrix_eq10()
}






#' Compute R from fractional changes of steady states,
#' but do the linear combination only after the derivatives have been taken
#'
#' @param pars_opt Free parameters a_ij to control the influence of the complexes
#' @param perturbation_prediction Simulated steady state perturbation data
#' @param obs_fun Observation function which returns the communicating species
#' @param p_fun Parameter transformation which returns different conditions
#' @param pars Parameters of the ODE model
#'
#' @return
#' @export
#'
#' @examples
R_fun_2 <- function(pars_opt,
                  perturbation_prediction,
                  obs_fun,
                  p_fun,
                  pars) {


  mypars <- pars
  mypars[names(pars_opt)] <- pars_opt
  mypars <- p_fun(mypars, deriv = F)

  myderiv <- lapply(2:length(perturbation_prediction), function(i) {
    Ctr <- perturbation_prediction[[1]]
    perturbed <- perturbation_prediction[[i]]
    return(2*(perturbed-Ctr)/(Ctr+perturbed))   # calculate the fractional derivative for each species individually
  }) %>%
    c(perturbation_prediction[1],.) %>%      # append the control condition, in order to be able to apply p_fun() conveniently
    set_names(attr(p_fun, "conditions"))

  myderiv <- lapply(seq_along(myderiv), function(i) {      # apply obs_fun()
    obs_fun(out = myderiv[[i]],
            pars = mypars[[i]],
            deriv = F)
  }) %>%
    extract(-1) %>%                             # kick out the control condition, which doesn't contain derivatives
    sapply(function(i) i[[1]][-1])                   # kick out "time", make matrix out of it and transpose it

  return(myderiv)
}
