
#' Convenient simulation of perturbation experiments
#'
#' A parameter transformation function to conveniently generate conditions, in which on parameter was perturbed
#'
#' @param pars_perturbed A named vector with the changes on the log scale of the perturbed parameters.
#' @param pars All parameters of the model.
#'
#' @return A parameter trafo function of class parfn, which goes from pars_outer (log scale) to pars_outer (log scale)
#' @export
#'
#' @examples
p_pert_fun <- function(pars_perturbed = pars_perturbed_0,
                       pars = pars_0) {


  p_control <- P(structure(names(pars), names = names(pars)),
                 condition = "Ctr",
                 compile = T,
                 modelname = "p_ctr")

  p_pert_list <- lapply(seq_along(pars_perturbed), function(i) {
    mytrafo <- structure(names(pars), names = names(pars))

    mytrafo[names(pars_perturbed)[i]] <- paste0(names(pars_perturbed)[i], " + 1 * (", pars_perturbed[i], ")")

    p <- P(mytrafo,
           condition = paste0(names(pars_perturbed)[i], " + ", round(pars_perturbed[i],1)),
           compile = T,
           modelname = paste0("p_pert",i))
  })

  p_pert <- NULL
  p_pert <- p_pert + p_control
  for(i in 1:length(p_pert_list)) p_pert <- p_pert + p_pert_list[[i]]

  return(p_pert)
}


#' Reduce ic_list to include only pars with more than one value
#'
#' This can be used to clean the output, such that the condition description is not as long anymore.
#'
#' @param ic_list
#'
#' @return list(ic_list, modified pars)
#' @export
#'
#' @examples
#'
#' pars <- c(a=5,b=1)
#' ic_list <- c(a=1,b=1:2)
#' pars[names(reduce_ic_list(ic_list)[[2]])] <- reduce_ic_list(ic_list)[[2]]
#' ic_list <- reduce_ic_list(ic_list)[[1]]
reduce_ic_list <- function(ic_list) {
  lengths <- sapply(ic_list, length)
  myic_list <- ic_list[lengths != 1]
  mypars <- ic_list[lengths == 1] %>% unlist
  return(list(myic_list,mypars))
}

#' Automatically generaty parameter trafos to quickly test different parameter values as conditions
#'
#' This function can be used, if one wants quickly wants to check the behaviour of the model with different parameters.
#' The named list ic_list contains all the parameter values that one wants to try out.
#' The function then returns a set of parameter transformations p_ic which can be put into g*x*p_ic to have a look at the model prediction.
#'
#' @param ic_list Named list of parameter values to be set.
#' @param pars The original set of parameters to be replaced.
#' @param mixed Should all combinations of the pars be gone through (TRUE)
#' or should the other pars be set to their first value supplied, when scanning through a certain parameter?
#'
#' @return A "parfn" of all conditions.
#' @export
#'
#' @examples
p_ic_fun <- function(ic_list, pars = pars_0, mixed = FALSE) {


  if(mixed) {
    mydf <- expand.grid(ic_list)
  } else {
    standardvalues <- lapply(ic_list, "[", 1)
    mydf <- lapply(seq_along(ic_list), function(i) {
      myic_list <- c(standardvalues[-i], ic_list[i])
      myic_list <-  myic_list[names(ic_list)]
      mygrid <- expand.grid(myic_list)
      if(i>1) mygrid <- mygrid[-1,]
      return(mygrid)
        }) %>% do.call(rbind,.)
  }

  p_ic_list <- lapply(1:nrow(mydf), function(i) {
    mytrafo <- structure(names(pars), names = names(pars))
    myics <- mydf[i,] %>% as.character %>% structure(names = names(ic_list))
    mytrafo[names(myics)] <- myics
    p <- P(mytrafo,
           condition = paste0(names(myics), "=", myics, collapse = "_"),
           compile = T,
           modelname = paste0("p_", sample(c(letters,0:9))))
  })

  p_ic <- NULL
  for(i in 1:length(p_ic_list)) p_ic <- p_ic + p_ic_list[[i]]

  return(p_ic)
}

#' Generate a list of parameter trafos.
#'
#' This function is similar to p_ic_fun, but it returns a list of parfns.
#' `p_ic_list` needs to be used like `lapply(p_ic_list, function(p_ic) {...(g*x*p_ic)})`.
#' This allows to do optimization over more than one condition (e.g. replace the dots by `obj`).
#'
#' @param ic_list Named list of parameter values to be set.
#' @param pars The original set of parameters to be replaced.
#' @param mixed Should all combinations of the pars be gone through (TRUE)
#' or should the other pars be set to their first value supplied, when scanning through a certain parameter?
#'
#' @return
#' @export
#'
#' @examples
p_ic_list_fun <- function(ic_list, pars = pars_0, mixed = FALSE) {

  if(mixed) {
    mydf <- expand.grid(ic_list)
  } else {
    standardvalues <- lapply(ic_list, "[", 1)
    mydf <- lapply(seq_along(ic_list), function(i) {
      myic_list <- c(standardvalues[-i], ic_list[i])
      myic_list <-  myic_list[names(ic_list)]
      mygrid <- expand.grid(myic_list)
      if(i>1) mygrid <- mygrid[-1,]
      return(mygrid)
    }) %>% do.call(rbind,.)
  }

  p_ic_list <- lapply(1:nrow(mydf), function(i) {
    mytrafo <- structure(names(pars), names = names(pars))
    myics <- mydf[i,] %>% as.character %>% structure(names = names(ic_list))
    mytrafo[names(myics)] <- myics
    p <- P(mytrafo)
    condition <- paste0(names(myics), "=", myics, collapse = "_")

    return(structure(list(p), names = condition))
  }) %>% do.call(c,.)

  return(p_ic_list)
}


#' Generate a log-trafo
#'
#' Just a convenience-function to quickly generate a log-trafo.
#'
#' @param parameters
#' @param condition
#' @param attach.input
#' @param keep.root
#' @param compile
#' @param modelname
#' @param method
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
P_log <- function(parameters = NULL, condition = NULL,
                  attach.input = FALSE, keep.root = TRUE, compile = FALSE,
                  modelname = NULL, method = c("explicit", "implicit"), verbose = FALSE) {

  trafo <- paste0("exp(log(", parameters, "))") %>% set_names(parameters)
  return(P(trafo, parameters = NULL, condition = NULL,
           attach.input = FALSE, keep.root = TRUE, compile = FALSE,
           modelname = NULL, method = c("explicit", "implicit"), verbose = FALSE)
  )
}
