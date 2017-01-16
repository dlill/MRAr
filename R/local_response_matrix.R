#' Print the local response matrix without R and steady_state attributes.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
print.r_small <- function(x) {
  attr(x, "R") <- NULL
  attr(x, "steady_state") <- NULL
  attr(x, "which_pars_perturbed") <- NULL
  attr(x, "class") <- "matrix"
  print(x)
}

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

  which_pars_perturbed <- colnames(R)

  condition_number <- R_svd$d[1]/R_svd$d[length(R_svd$d)]
  if(abs(condition_number) > 1e12) warning("R is ill conditioned. Condition number is ", condition_number)

  R_inv <- R_svd$v %*% diag(1/R_svd$d) %*% t(R_svd$u)

  R_inv_diag <- diag(R_inv)

  r <-  - R_inv / R_inv_diag
  dimnames(r) <- list(rownames(R), rownames(R))
  attr(r, "which_pars_perturbed") <- which_pars_perturbed
  class(r) <- c("r_small", "matrix")
  return(r)
}

#' Title
#'
#' This function only works for one condition
#'
#' @param pars
#' @param time_max
#' @param predfun
#' @param which_pars_perturbed
#' @param pars_optimization
#'
#' @return
#' @export
#'
#' @examples
r_numDeriv_fun <- function(pars = pars_0,
                           predfun = (g*x),
                           deparsed_predfun = NULL,
                           which_pars_perturbed = which_pars_perturbed_0,
                           pars_optimization = NULL) {

  if(!is.null(pars_optimization)) pars[names(pars_optimization)] <- pars_optimization

  conditions <- attr(predfun, "conditions")
  if(is.null(conditions)) conditions <- 1

  if(is.null(deparsed_predfun)) {deparsed_predfun <- deparse(substitute(predfun))}


  # define a function that puts out the steady states as a function of pars_perturbed

  reduced_steady_states <- function(pars_perturbed) {
    mypars <- pars
    mypars[which_pars_perturbed] <- pars_perturbed
    steady_states_numerical(pars = mypars, deparsed_predfun = deparsed_predfun)[[1]][[1]][,-1]  #first condition ([[1]]) without time column ([,-1])
  }

  pars_perturbed <- pars[which_pars_perturbed]

  myderiv <- numDeriv::jacobian(reduced_steady_states, pars_perturbed)

  dimnames(myderiv) = list(modules, which_pars_perturbed)

  r_sens <- local_response_matrix_eq10(myderiv)
  attr(r_sens, "R") <- myderiv
  attr(r_sens, "steady_state") <- steady_states_numerical(pars = pars, predfun = predfun, deparsed_predfun = deparsed_predfun)[[1]]

  out <- list(r_sens)
  names(out) <- conditions

  return(out)
}



#' Retrieve the local response matrix from the sensitivity equations
#'
#' @param pars The parameters for which r shall be computed
#' @param which_pars_perturbed The perturbations
#' @param time_max First guess for the time the steady state has been reached
#' @param predfun Prediction function
#' @param pars_optimization
#'
#' @details The prediction function can also incorporate more than one condition.
#'
#' @return A conditions-named list of local response matrices "r" with attributes "R" and "steady state".
#' @export
#'
#' @examples
r_sens_fun <- function(pars = pars_0,
                       predfun = (g*x),
                       deparsed_predfun = NULL,
                       which_pars_perturbed = which_pars_perturbed_0,
                       pars_optimization = NULL) {

  # time_max <- steady_states_numerical(pars = pars, time_max = time_max) %>% attr("counter_time_max") %>% extract(2)
  # times <- seq(1,time_max,by = 100)

  conditions <- attr(predfun, "conditions")
  if(is.null(conditions)) conditions <- 1

  if(is.null(deparsed_predfun)) {deparsed_predfun <- deparse(substitute(predfun))}
  steady_states <- steady_states_numerical(pars = pars, deparsed_predfun = deparsed_predfun)

  time_max <- lapply(conditions, function(cond) steady_states[[cond]][[1]][1]) %>% do.call(c,.) %>% max
  times <- seq(1,time_max, length.out = 200)

  if(!is.null(pars_optimization)) pars[names(pars_optimization)] <- pars_optimization

  mypred <- predfun(times, pars, deriv = T)

  out <- lapply(conditions, function(cond) {

    my_steady_state <- mypred[[cond]][length(times),]

    myderiv <- attr(mypred[[cond]], "deriv")[length(times),]
    derivnames <- outer(modules, which_pars_perturbed, paste, sep = ".")

    myderiv <- myderiv[derivnames] %>% matrix(nrow = length(modules), ncol = length(which_pars_perturbed))
    dimnames(myderiv) = list(modules, which_pars_perturbed)

    r_sens <- local_response_matrix_eq10(myderiv)
    attr(r_sens, "R") <- myderiv
    attr(r_sens, "steady_state") <- my_steady_state

    return(r_sens)
  })

  names(out) <- conditions

  return(out)

}







#' Determine which signs have flipped over by including the complexes
#'
#' To do: have a look at what happens, if pars_opt = .5 or combinations of pars_opt (eg c(0,1), c(1,0) ...)
#'
#' @param pars_optimization
#' @param pars
#' @param time_max
#' @param predfun
#' @param which_pars_perturbed
#'
#' @return
#' @export
#'
#' @examples
r_sign_changed_fun <- function(pars = pars_0,
                               time_max = 12000,
                               predfun = (g*x),
                               deparsed_predfun = NULL,
                               which_pars_perturbed = which_pars_perturbed_0,
                               pars_optimization = pars_opt_0,
                               r_fun = r_numDeriv_fun) {

  if(is.null(deparsed_predfun)) {deparsed_predfun <- deparse(substitute(predfun))}
  pars_optimization <- structure(rep(0, length(pars_optimization)), names = names(pars_optimization))

  r_0 <- r_fun(pars = pars,
               time_max = time_max,
               predfun = predfun,
               deparsed_predfun = deparsed_predfun,
               which_pars_perturbed = which_pars_perturbed,
               pars_optimization = pars_optimization)

  r_1 <- r_fun(pars = pars,
               time_max = time_max,
               predfun = predfun,
               deparsed_predfun = deparsed_predfun,
               which_pars_perturbed = which_pars_perturbed,
               pars_optimization = pars_optimization + 1)

  conditions <- attr(predfun, "conditions")
  if(is.null(conditions)) conditions <- 1;

  out <- lapply(conditions, function(cond) {
    (sign(r_0[[cond]]) - sign(r_1[[cond]])) %>% as.logical() %>% matrix(nrow = nrow(r_0[[cond]]))
  })

  names(out) <- conditions

  return(out)
}

