#' Run the optimization procedure at many conditions
#'
#' Even though this function takes only one parameter as input, it depends on many objects to be defined.
#'
#' xs(times, c(pars_raw, ic_raw))
#' (g*xs)(times, c(pars_raw, ic_raw, pars_opt))
#' pars_0, pars_opt_0
#'
#' which_pars_perturbed_0, p_pert
#'
#'
#' @param p_ic_list
#'
#' @return
#' @export
#'
#' @examples
run_simulations <- function(ic_list = NULL,
                            nfits= 10,
                            pars_opt = pars_opt_0,
                            obs_fun = g,
                            pars = pars_0,
                            which_pars_perturbed = which_pars_perturbed_0,
                            fixed = NULL) {

  
  p_pert <- p_pert_fun(which_pars_perturbed = which_pars_perturbed,
                       perturbation_strength = 1.1,
                       pars = pars)
  
  if(is.null(ic_list)) p_ic_list <- list(P(structure(names(pars), names = names(pars)))) else p_ic_list <- p_ic_list_fun(ic_list)
  

  out <- lapply(seq_along(p_ic_list), function(ic) {

    
    # Changed parameter trafo function
    cat("Condition ", ic, " of ", length(p_ic_list), "\n")
    p_fun <- (p_pert * p_ic_list[[ic]])
    p_fun(pars)
    # Perturbation Data generation,
    perturbation_prediction <-  (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)
    r_kept <- r_kept_fun(pars_opt = pars_opt,
                           perturbation_prediction = perturbation_prediction,
                           obs_fun = obs_fun,
                           p_fun = p_fun,
                           pars = pars)

    # Fitting
    myfits <- mstrust(obj_alpha,
                      center =  pars_opt+0.5,
                      studyname = "Fits",
                      cores = 3,
                      fits = nfits,
                      sd = 0.3,

                      perturbation_prediction = perturbation_prediction,
                      r_kept = r_kept,
                      obs_fun = obs_fun,
                      p_fun = p_fun,
                      mypars = pars,
                      fixed = NULL)

    # print(myfits)
    # Analysis of fits and output preparation
    best_fit <- try({
      bf <- myfits %>% as.parframe()
      bf[abs(bf$a_1) < 1 & abs(bf$a_2) < 1, ] %>% as.parvec()
      })

    # If no fits are in parameter range <1, take the best fit
    if(inherits(best_fit, "try-error")) 
      best_fit <- try({myfits %>% as.parframe() %>% as.parvec()})
    
    # if there are still problems, return just the fitlist and the condition
    if(inherits(best_fit, "try-error")) {return(list(all_fits = myfits %>% as.parframe,
                                                     best_fit = 0,
                                                     r_0 = matrix(0),
                                                     r_kept = matrix(0),
                                                     r_best_fit = matrix(0),
                                                     steady_states = matrix(0),
                                                     prediction = matrix(0),
                                                     condition = names(p_ic_list)[ic],
                                                     which_pars_perturbed = which_pars_perturbed))
    } else
    r_0 <- R_fun(pars_opt = pars_opt,
                 perturbation_prediction = perturbation_prediction,
                 obs_fun = obs_fun,
                 p_fun = p_fun,
                 pars = pars) %>% local_response_matrix_eq10()

    r_best_fit <- R_fun(pars_opt = best_fit,
                        perturbation_prediction = perturbation_prediction,
                        obs_fun = obs_fun,
                        p_fun = p_fun,
                        pars = pars) %>% local_response_matrix_eq10()

    steady_states <- (xs*p_ic_list[[ic]])(times = c(0,Inf), pars = pars)

    prediction <- (x*p_ic_list[[ic]])(times = seq(0,steady_states[[1]][1], length.out = 50), pars = pars)

    return(list(all_fits = myfits %>% as.parframe,
                best_fit = best_fit,
                r_0 = r_0,
                r_kept = r_kept,
                r_best_fit = r_best_fit,
                steady_states = steady_states,
                prediction = prediction,
                condition = names(p_ic_list)[ic],
                which_pars_perturbed = which_pars_perturbed))
  })

  names(out) <- names(p_ic_list)

  return(out)
}
