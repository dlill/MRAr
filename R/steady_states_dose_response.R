#' Steady states out of a simulation of a system
#'
#' @description This function gives access to the runsteady-function of deSolve within the formalism of dMod.
#' One can pass a prediction function like (g*x*p), (g*x) or (x*p).
#' The steady states of the ODE are computed by deSolve::runsteady and g and p are executed if they exist.
#'
#' @param pars The parameters for which the steady states are to be computed.
#' @param predfun An expression like (g*x*p), (g*x) or (x*p). Brackets can be omitted.
#' @param odemodel The object in which the odemodel is stored.
#'
#' @return
#' @export
#'
#' @examples
steady_states_numerical <- function(pars = pars_0, predfun = (g*x), deparsed_predfun = NULL) {

  if(is.null(deparsed_predfun)) {deparsed_predfun <- deparse(substitute(predfun))}
  deparsed_predfun <- str_replace_all(deparsed_predfun, "[() ]", "")
  deparsed_predfun <- str_split(deparsed_predfun, "\\*")[[1]]
  g_ss <- tryNULL(get(deparsed_predfun[str_detect(deparsed_predfun, "g")]))
  p_ss <- tryNULL(get(deparsed_predfun[str_detect(deparsed_predfun, "p")]))

  # cat("\n deparsed_predfun: \n",
  #     deparsed_predfun, "\n \n")

  if(is.null(p_ss)) {
    pars_list <- list(pars)
    conditions <- 1
  } else {
    pars_list <- p_ss(pars)
    conditions <- attr(p_ss, "conditions")
  }


  steadys <- lapply(pars_list, function(pars){
    ic    <- structure(as.numeric(pars[names(ic_raw)]), names = names(ic_raw))
    parms <- structure(as.numeric(pars[names(pars_raw)]), names = names(pars_raw))

    # loadDLL(myodemodel$func)
    mysol  <- runsteady(y = ic, times = c(0,Inf), func = "derivs", parms = parms, dllname = myodemodel$func, cinitfunc = "initmod")

    # apply observation function to the steady state
    my_ss <- matrix(c(attr(mysol, "time"),mysol$y), nrow = 1, dimnames = list(NULL, c("time", names(ic_raw))))

    if(!is.null(g_ss)) {
      my_ss <- g_ss(my_ss, pars = pars, deriv = F)
    }
    return(my_ss)
  })
  names(steadys) <- conditions
  return(steadys)
}


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
#'  It is a condition-named list of steady_state concentrations with attributes counter_time_max =  c(counter, time_max)
#'  If count == 4, the predictions are also appended as attributes to look at the predictions, since in this case, the original time_max_guess has been already increased 1000 times, but steady state has not been reached.
#' @export
#'
#' @examples
steady_states_numerical_2 <- function(pars = pars_0,  time_max = 12000, predfun = (g*x), by = 10) {#   flag <- TRUE
  counter <- 1
  while (flag) {
    times <- seq(1,time_max, by = by)
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
        steady_state <- tail_of_mypred[2,]
        if(counter == 4 ) {attr(steady_state, "prediction") <- mypred[[cond]]
          warning("In condition", cond, "the steady state criterion has not been met. Look at the prediction (attribute) to see what's happening.")
          }
        return(steady_state) # steady state values
      }
    })
  }
  names(out) <- conditions
  attr(out, "counter_time_max") <- c(counter, time_max)
  return(out)
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
dose_response_fun <- function(par_dose, dosages, pars, predfun = (g*x), depased_predfun = NULL) {

  if(is.null(deparsed_predfun)) {deparsed_predfun <- deparse(substitute(predfun))}

  d_r <- mclapply(dosages, function(i) {
    mypars <- pars
    mypars[par_dose] = i
    conditions <- attr(predfun, "conditions")
    my_steady_states <-
      data.frame( condition = conditions, par_dose = i, do.call(rbind, steady_states_numerical(pars = mypars, deparsed_predfun = depased_predfun)))
  }, mc.cores = 4) %>% do.call(rbind, .)
  dose_response_colnames <- d_r %>% colnames
  dose_response_colnames <- dose_response_colnames[-c(1:2)]

  d_r <- gather_(data = d_r, key_col = "observable", value_col =  "concentration", gather_cols = dose_response_colnames, factor_key = T)
  return(d_r)

}








