#' Simulate the System and take the last value of the simulation as the values for the steady states
#'
#' The system is assumed to have reached steady state,
#' if the relative changes of the observables are less than 1e-3
#'
#' @param pars The parameters needed for the prediction function
#' @param time_max Initial guess at which time steady state is reached
#' @param predfun The prediction function. Right now, it should contain only one condition.
#'
#' @details
#'
#' @return The steady state concentrations that were obtained during the simulation.
#' @export
#'
#' @examples
steady_states_numerical <- function(pars,  time_max = 12000, predfun = (g*x)) {
  flag <- TRUE
  counter <- 1
  while (flag) {
    times <- seq(1,time_max, by = 10)
    mypred <- predfun(times, pars, deriv = F)

    conditions <- attr(predfun, "conditions")
    if(is.null(conditions)) conditions <- 1;
    # quickly check, if steady state has been reached
    out <- lapply(conditions, function(cond) {
      tail_of_mypred <- mypred[[cond]][(length(times)-1) :length(seq(times)),-1]

      diffs_bigger_than_threshold <- ((tail_of_mypred[1,]/tail_of_mypred[2,]) - 1)  %>%
        abs() %>%
        is_greater_than(1e-3)

      if(any(diffs_bigger_than_threshold) & counter < 4) {
        time_max <<- 10*time_max
        counter <<- counter + 1
      } else {
        flag <<- FALSE
        return(tail_of_mypred[2,]) # steady state values
      }
    })
  }
  names(out) <- conditions
  attr(out, "counter_time_max") <- c(counter, time_max)
  return(out)
}


#' Global response matrix R as a function of the perturbed parameters
#'
#' @param pars_perturbed
#' @param pars
#' @param logderiv Logical. TRUE: Compute d ln(x_i) /dp_j, FALSE: d x_i /dp_j
#'
#' @return
#' @export
#'
#' @examples
global_response_matrix <- function(pars_perturbed, pars, time_max = 12000, logderiv = FALSE) {

  steady_states_0 <- steady_states_numerical(pars = pars, time_max = time_max)[[1]] #right now just for one condition

  sapply(names(pars_perturbed), function(i) {

    delta_p <- pars_perturbed[i]-pars[i]
    new_pars <- pars
    new_pars[i] <- pars_perturbed[i] # Apply the perturbation

    new_steady_states <- steady_states_numerical(pars = new_pars, time_max = time_max)[[1]]

    deriv <- (new_steady_states-steady_states_0)/delta_p

    if(logderiv) deriv <- deriv / steady_states_0

    return(deriv)
  })

}


#' Global response matrix R
#'
#' @param pars_perturbed
#' @param pars
#'
#' @description This version calculates d ln x_i/dp_j as suggested in the paper
#'
#' @return
#' @export
#'
#' @examples
global_response_matrix_fractional <- function(pars_perturbed, pars) {

  steady_states_0 <- steady_states_numerical(pars = pars)[[1]] #right now just for one condition

  sapply(names(pars_perturbed), function(i) {

    new_pars <- pars
    new_pars[i] <- pars_perturbed[i] # Apply the perturbation

    new_steady_states <- steady_states_numerical(pars = new_pars)[[1]]

    logderiv <- 2*(new_steady_states-steady_states_0)/(new_steady_states+steady_states_0)

    return(logderiv)
  })

}




#' Calculate the local response matrix r from the global responses
#'
#' @param R
#'
#' @details
#'
#' @return r The local response matrix
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



#' Retrieve the local response matrix from the sensitivity equations
#'
#' @param pars The parameters for which r shall be computed
#'
#' @return
#' @export
#'
#' @examples
r_sens_fun <- function(pars_perturbed, pars, time_max = 12000) {

  # time_max <- steady_states_numerical(pars) %>% attr("counter_time_max") %>% extract(2)
  times <- seq(1,time_max,by = 10)
  mypred <- (g*x)(times, pars_0, deriv = T)
  mysteady_state <- mypred[[1]][length(times),]

  myderiv <- attr(mypred[[1]], "deriv")[length(times),]
  derivnames <- outer(modules, names(pars_perturbed), paste, sep = ".")

  myderiv <- myderiv[derivnames] %>% matrix(nrow = length(modules), ncol = length(pars_perturbed))
  dimnames(myderiv) = list(modules, names(pars_perturbed))

  r_sens <- local_response_matrix_eq10(myderiv)
  attr(r_sens, "R") <- myderiv
  return(r_sens)
}



#' Local response matrix r as a function of the some explicitly specified parameters
#'
#' @description The local response matrix r is computed with the standard set of parameters "pars".
#' If a parameter is given in pars_optimization, this value is taken instead.
#'
#' @param pars_optimization
#'
#' @return
#' @export
#'
#' @examples
r_opt_fun <- function(pars_optimization,pars, time_max = 12000) {

  newpars <- pars
  newpars[names(pars_optimization)] = pars_optimization

  return(global_response_matrix(pars_perturbed = pars_perturbed, pars = newpars, time_max = time_max) %>% local_response_matrix_eq10)
}



#' Determine which signs have flipped over by including the complexes
#'
#' To do: have a look at what happens, if pars_opt = .5 or combinations of pars_opt (eg c(0,1), c(1,0) ...)
#'
#' @param pars_optimization
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
r_sign_changed_fun <- function(pars_optimization, pars) {
  pars_optimization <- structure(rep(0, length(pars_optimization)), names = names(pars_optimization))
  r_0 <- r_opt_fun(pars_optimization = pars_optimization, pars = pars)
  r_1 <- r_opt_fun(pars_optimization = pars_optimization+1, pars = pars)
  # r_05 <- r_opt_fun(pars_optimization = pars_optimization+0.5, pars = pars)

  out <- (sign(r_0) - sign(r_1)) %>% as.logical() %>% matrix(nrow = nrow(r_0))
  return(out)
}


#' Objective function to minimize certain elements of the local response matrix
#'
#' This function squares all elements that flip sign once the complexes are added to the observation.
#' This can then be used to minimize the values.
#'
#' @param pars_optimization The parameters which are needed to optimize the matrix
#' @param r_sign_changed A matrix of logicals specifying the elements to minimize.
#'
#' @return The parameter values for which the sign-changing elements are minimal.
#' @export
#'
#' @examples
obj_frobenius <- function(pars_optimization, obj_pars, r_sign_changed) {

  elements_fun  <- function(pars_optimization) r_opt_fun(pars_optimization = pars_optimization, pars = pars)[r_sign_changed]
  elements      <- elements_fun(pars_optimization)
  elements_jac  <- numDeriv::jacobian(elements_fun, pars_optimization)

  value   <- elements^2 %>% sum()
  grad    <- 2*(t(elements)%*%elements_jac) %>% as.numeric %>% structure(names = names(pars_optimization))
  # hessian <- 2*(t(elements_jac) %*% elements_jac) %>% structure(dimnames = list(names(pars_optimization),names(pars_optimization)))
  hessian <- matrix(0, nrow = length(pars_optimization), ncol = length(pars_optimization)) %>% structure(dimnames = list(names(pars_optimization),names(pars_optimization)))

  return(list(value = value, gradient = grad, hessian = hessian))
}


#' Compute dose-response curves
#'
#' Compute the steady state dose-response of a parameter.
#' Put the dose response out in a tidy data format.
#'
#' @param par_dose Character of length 1: Which par to perturb
#' @param dosages Numeric including the dosages eq(1,50,by = 1)
#' @param pars
#' @param predfun
#'
#' @return
#' @export
#'
#' @examples
dose_response_fun <- function(par_dose, dosages, pars, predfun = (g*x*p)) {

  d_r <- mclapply(dosages, function(i) {
    mypars <- pars
    mypars[par_dose] = i
    conditions <- attr(predfun, "conditions")
    my_steady_states <-
      data.frame( condition = conditions, par_dose = i, do.call(rbind, steady_states_numerical(mypars,time_max = 20000, predfun = predfun)))
  }, mc.cores = 4) %>% do.call(rbind, .)
  dose_response_colnames <- d_r %>% colnames
  dose_response_colnames <- dose_response_colnames[-c(1:2)]

  d_r <- gather_(data = d_r, key_col = "observable", value_col =  "concentration", gather_cols = dose_response_colnames, factor_key = T)
  return(d_r)

}


#' Automatically generaty parameter trafos to quickly test different parameter values as conditions
#'
#' This function can be used, if one wants quickly wants to check the behaviour of the model with different parameters.
#' The named list ic_list contains all the parameter values that one wants to try out.
#' The function then returns a set of parameter transformations p_ic which can be put into g*x*p_ic to have a look at the model prediction.
#'
#' @param ic_list Named list of parameter values to be set.
#' @param pars The original set of parameters to be replaced.
#'
#' @return
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







