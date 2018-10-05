#' Run the algorithm
#'
#' This function has a lot of arguments, to ensure that each object is the correct one
#'
#' @param pars pars. The parameters used to sim
#' @param which_pars_perturbed character. which pars are used to perturb the system?
#' @param which_alpha_pars names of the alpha pars used, i.e. which complexes are measured?
#' @param alpha_pars the parameters which are used for weighting the complexes
#' @param alpha numeric 1L. The value of the free parameters, at which the two matrices are compared to see if a sign switch occured.
#' @param g observation function
#' @param xs steady state prediction function
#' @param p_log parameter trafo
#'
#' @return a tibble with r0, r1, rkep and ropt
#' @export
algo <- function(pars,
                 which_pars_perturbed,
                 which_alpha_pars,
                 alpha_pars,
                 alpha,
                 g,
                 xs,
                 p_log) {
  if (any(!(which_alpha_pars %in% names(alpha_pars))))
    stop("some alpha_pars in alpha_pars_settings which don't exist in alpha_pars")

  myperturbation <- rep(log(0.9), length(which_pars_perturbed)) %>% `names<-`(which_pars_perturbed)
  p_pert <- p_pert_fun(myperturbation, pars = pars, modelname = NULL)
  perturbation_prediction <- (xs*p_log*p_pert)(c(0,Inf), pars, deriv = F)

  out <- tibble(r_0 = list(NULL),
                r_1 = list(NULL),
                r_kept = list(NULL),
                parframes = list(NULL),
                r_opt = list(NULL))


  r_0 <- R_fun(pars_opt = alpha_pars[which_alpha_pars] * 0 + log(.Machine$double.eps),
               perturbation_prediction = perturbation_prediction,
               obs_fun = g,
               p_fun = (p_log*p_pert),
               pars = pars) %>% local_response_matrix_eq10()
  out$r_0 <- list(r_0 %>% round(2))

  r_1 <- R_fun(pars_opt = alpha_pars[which_alpha_pars] * 0 + alpha,
               perturbation_prediction = perturbation_prediction,
               obs_fun = g,
               p_fun = (p_log*p_pert),
               pars = pars) %>% local_response_matrix_eq10()
  out$r_1 <- list(r_1 %>% round(2))

  r_kept <- r_kept_fun2(r_0, r_1)
  out$r_kept <- list(r_kept)

  myfits <- mstrust(obj_alpha,
                    center =  structure(rep(0, length(which_alpha_pars)), names = which_alpha_pars),
                    studyname = "Fits",
                    cores = 1,
                    fits = 3,
                    sd = 1,
                    mypars = pars,
                    perturbation_prediction = perturbation_prediction,
                    r_kept = r_kept,
                    p_fun = (p_log * p_pert),
                    obs_fun = g)
  myparframe <- try(myfits %>% as.parframe())
  out$parframes <- list(myparframe)
  if (inherits(myparframe, "try-error")) {
    out$parframes <- list(myfits)
    return(out)
  }

  r_opt <- r_alpha_fun(pars_opt = myfits %>% as.parframe() %>% as.parvec() %>% unclass() ,
                       pars = pars,
                       perturbation_prediction = perturbation_prediction,
                       p_fun = (p_log*p_pert),
                       obs_fun = g)


  out$r_opt <- list(r_opt %>% round(2))
  return(out)
}
