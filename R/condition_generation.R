
#' Convenient simulation of perturbation experiments
#'
#' A parameter transformation function to conveniently generate conditions, in which on parameter was perturbed
#'
#' @param which_pars_perturbed The names of the parameters to be perturbed.
#' @param perturbation_strength The fold-change of the parameter on linear scale.
#' Either one value to set all of them to the same strength or a numeric vector with the individual perturbation strengths.
#' @param pars All parameters of the model.
#'
#' @return A parameter trafo function of class parfn.
#' @export
#'
#' @examples
p_pert_fun <- function(which_pars_perturbed = which_pars_perturbed_0,
                       perturbation_strength = 1.1,
                       pars = pars_0) {

  if(length(perturbation_strength) == 1)
    perturbation_strength <- rep_len(perturbation_strength, length.out = length(which_pars_perturbed))
  if(length(perturbation_strength) < length(which_pars_perturbed))
    stop("Number of perturbation_strength strengths not equal to number of perturbed parameters.")


  p_control <- P(structure(names(pars), names = names(pars)), condition = "Ctr")

  p_pert_list <- lapply(seq_along(which_pars_perturbed), function(i) {
    mytrafo <- structure(names(pars), names = names(pars))

    mytrafo[which_pars_perturbed[i]] <- paste0(which_pars_perturbed[i], "*", perturbation_strength[i])

    p <- P(mytrafo, condition = paste0(which_pars_perturbed[i], "*", perturbation_strength[i]))
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
#'
#' @return A "parfn" of all conditions.
#' @export
#'
#' @examples
p_ic_fun <- function(ic_list, pars = pars_0) {

  mydf <- expand.grid(ic_list)

  p_ic_list <- lapply(1:nrow(mydf), function(i) {
    mytrafo <- structure(names(pars), names = names(pars))
    myics <- mydf[i,] %>% as.character %>% structure(names = names(ic_list))
    mytrafo[names(myics)] <- myics
    p <- P(mytrafo, condition = paste0(names(myics), "=", myics, collapse = "_"))
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
#'
#' @return
#' @export
#'
#' @examples
p_ic_list_fun <- function(ic_list, pars = pars_0) {

  mydf <- expand.grid(ic_list)

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
