# Some adaptions of dMod and cOde functions for steady states

#' Alternative version of expand.grid
#' @param seq1 Vector, numeric or character
#' @param seq2 Vector, numeric or character
#' @return Matrix ob combinations of elemens of \code{seq1} and \code{seq2}
expand.grid.alt <- function(seq1, seq2) {
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}


#' Interface to runsteady()
#'
#' @param y named vector of type numeric. Initial values for the integration
#' @param times vector of type numeric. Integration times
#' @param func return value from funC()
#' @param parms named vector of type numeric.
#' @param ... further arguments going to \code{ode()}
#' @details See deSolve-package for a full description of possible arguments
#' @return matrix with times and states
#' @export
runsteadyC <- function(y, times=c(0,Inf), func, parms, ...)
{
  nGridpoints <- attr(func, "nGridpoints")
  times.inner <- times
  which.times <- match(times, times.inner)
  yout <- c(attr(func, "forcings"), names(attr(func, "outputs")))
  y <- y[attr(func, "variables")]
  parms <- parms[attr(func, "parameters")]
  parms <- c(parms, rep(0, length(y)))
  arglist <- list(y = y, times = times.inner, func = paste0(func, "_derivs"),
                  parms = parms, dllname = func, initfunc = paste0(func, "_initmod"))
  if (attr(func, "jacobian") == "full")
    arglist <- c(arglist, list(jacfunc = paste0(func, "_jacobian")))
  if (attr(func, "jacobian") == "inz.lsodes") {
    inz <- attr(func, "inz")
    lrw <- 20 + 3 * dim(inz)[1] + 20 * length(y)
    arglist <- c(arglist, list(sparsetype = "sparseusr",
                               inz = inz, lrw = lrw))
  }
  if (attr(func, "jacobian") == "jacvec.lsodes") {
    inz <- attr(func, "inz")
    arglist <- c(arglist, list(sparsetype = "sparseusr",
                               jacvec = paste0(func, "_jacvec"), inz = inz))
  }
  if (!is.null(attr(func, "forcings")) & attr(func, "fcontrol") ==
      "nospline")
    arglist <- c(arglist, list(initforc = paste0(func, "_initforc")))
  if (!is.null(attr(func, "rootfunc")))
    arglist <- c(arglist, list(rootfunc = paste0(func, "_myroot"), nroot = length(attr(func,
                                                                        "rootfunc"))))
  if (!is.null(yout)) {
    arglist <- c(arglist, list(nout = length(yout), outnames = yout))
  }
  moreargs <- list(...)
  if (any(names(moreargs) == "forcings") & attr(func, "fcontrol") ==
      "einspline")
    moreargs <- moreargs[-which(names(moreargs) == "forcings")]
  i <- match(names(moreargs), names(arglist))
  is.overlap <- which(!is.na(i))
  is.new <- which(is.na(i))
  arglist[i[is.overlap]] <- moreargs[is.overlap]
  arglist <- c(arglist, moreargs[is.new])
  out <- do.call(rootSolve::runsteady, arglist)
  out <- c(time = attr(out, "time"), out[[1]])
  out <- matrix(out, nrow = 1, dimnames = list(NULL, names(out)))
  return(out)
}



#' Steady state prediction function for ODE models.
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function \code{x(times, pars, deriv = TRUE)} returning ODE output and sensitivities.
#' @param odemodel object of class \link{odemodel}
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings.
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See \link[deSolve]{events}.
#' Within \code{Xs()} a \code{data.frame} of additional events is generated to
#' reset the sensitivities appropriately, depending on the event method.
#' ATTENTION: The addional events are not dynamically recalculated. If you call the prediction
#' function with alternative events, the prediction is fine but the sensitivities can be wrong.
#' @param names character vector with the states to be returned. If NULL, all states are returned.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @return Object of class \link{prdfn}. If the function is called with parameters that
#' result from a parameter transformation (see \link{P}), the Jacobian of the parameter transformation
#' and the sensitivities of the ODE are multiplied according to the chain rule for
#' differentiation. The result is saved in the attributed "deriv",
#' i.e. in this case the attibutes "deriv" and "sensitivities" do not coincide.
#' @export
#' @import deSolve
Xs_steady <- function(odemodel, forcings=NULL, events=NULL, names = NULL, condition = NULL, optionsOde=list(method = "lsoda"), optionsSens=list(method = "lsodes")) {

  func <- odemodel$func
  extended <- odemodel$extended
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")

  myforcings <- forcings
  myevents <- events

  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  forcnames <- attr(func, "forcings")

  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables")%in%variables]
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) v[2]))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)

  # Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar

  #Additional events for resetting the sensitivities when events are supplied
  myevents.addon <- NULL
  if(!is.null(myevents)) {
    myevents.addon <- lapply(1:nrow(myevents), function(i) {
      newevent <- with(myevents[i, ], {
        newvar <- sensvar[senssplit.1 == var]
        newtime <- time
        newvalue <- switch(as.character(method), replace = 0, add = 0, multiply = value)
        newmethod <- method
        if(length(newvar) > 0 && method != "add") {
          data.frame(var = newvar, time = newtime, value = newvalue, method = newmethod)
        } else {
          NULL
        }
      })

      return(newevent)

    })
    myevents.addon <- do.call(rbind, myevents.addon)
  }

  # Names for deriv output
  sensGrid <- expand.grid(variables, c(svariables, sparameters), stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")

  # Only a subset of all variables/forcings is returned
  if (is.null(names)) names <- c(variables, forcnames)

  # Update sensNames when names are set
  select <- sensGrid[, 1] %in% names
  sensNames <- paste(sensGrid[,1][select], sensGrid[,2][select], sep = ".")


  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    names = names,
    events.addon = myevents.addon,
    optionsOde = optionsOde,
    optionsSens = optionsSens
  )

  P2X <- function(times, pars, deriv=FALSE){


    yini <- pars[variables]
    mypars <- pars[parameters]

    events <- controls$events
    forcings <- controls$forcings
    myevents.addon <- controls$events.addon
    optionsOde <- controls$optionsOde
    optionsSens <- controls$optionsSens
    names <- controls$names


    # Add event time points (required by integrator)
    event.times <- unique(events$time)
    times <- sort(union(event.times, times))

    myderivs <- NULL
    mysensitivities <- NULL
    if (!deriv) {

      # Evaluate model without sensitivities
      # loadDLL(func)
      if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
      out <- do.call(runsteadyC, c(list(y=yini, times=times, func=func, parms=pars, forcings=forc,events = list(data = events)), optionsOde))
      out <- submatrix(out, cols = c("time", names))
      #out <- cbind(out, out.inputs)


    } else {

      # Evaluate extended model
      # loadDLL(extended)
      if (!is.null(forcings)) forc <- setForcings(extended, forcings) else forc <- NULL
      outSens <- do.call(runsteadyC, c(list(y = c(unclass(yini), yiniSens), times = times, func = extended, parms = mypars,
                                            forcings = forc,
                                            events = list(data = rbind(events, myevents.addon))), optionsSens))
      #out <- cbind(outSens[,c("time", variables)], out.inputs)
      out <- submatrix(outSens, cols = c("time", names))
      mysensitivities <- submatrix(outSens, cols = !colnames(outSens) %in% c(variables, forcnames))


      # Apply parameter transformation to the derivatives
      variables <- intersect(variables, names)
      sensLong <- matrix(outSens[,sensNames], nrow = dim(outSens)[1]*length(variables))
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        sensLong <- sensLong %*% submatrix(dP, rows = c(svariables, sparameters))
        sensGrid <- expand.grid.alt(variables, colnames(dP))
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep = ".")
      }
      myderivs <- matrix(0, nrow = nrow(outSens), ncol = 1 + length(sensNames), dimnames = list(NULL, c("time", sensNames)))
      myderivs[, 1] <- out[, 1]
      myderivs[, -1] <- sensLong

    }

    #prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = unique(sensGrid[,2]))
    prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = pars)

  }

  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]


  prdfn(P2X, c(variables, parameters), condition)


}


#' Augmentation of getDerivs to subset wrt to which_states.which_pars
#'
#' @param prediction a prdframe which contains derivatives as an attribute
#' @param which_states which states shall be taken the derivative of ...
#' @param which_pars ... and with regard to which pars?
#'
#' @return
#' @export
#'
#' @examples
extract_derivs <- function(prediction, which_states = NULL, which_pars = NULL) {
  prediction %>%
    getDerivs() %>%
    lapply(function(i) {
      mynames <-  colnames(i)
      state_columns <- par_columns <- rep(T, length(mynames))
      if(!is.null(which_states)) {
        state_columns <- sapply(which_states, function(ws) {str_detect(mynames, paste0("^", ws, "\\."))}) %>% matrix(ncol = length(which_states)) %>% apply(1, any)
      }
      if(!is.null(which_pars)) {
        par_columns <- sapply(which_pars, function(wp) str_detect(mynames, paste0("\\.", wp, "$"))) %>% matrix(ncol = length(which_pars)) %>% apply(1, any)
      }

      cols <- state_columns & par_columns
      cols[1] <- TRUE
      return(i[,cols])
    })
}

