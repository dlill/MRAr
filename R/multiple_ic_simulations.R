#' Run the optimization procedure at many conditions
#'
#' Even though this function takes only one parameter as input, it depends on many objects to be defined.
#'
#' xs(times, c(pars_raw, ic_raw))
#' (g*xs)(times, c(pars_raw, ic_raw, pars_opt))
#' pars_0, pars_opt_0
#'
#' pars_perturbed_0, p_pert
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
                            cores = 3,
                            pars_opt = pars_opt_0,
                            obs_fun = g,
                            pars = pars_0,
                            pars_perturbed = pars_perturbed_0,
                            fixed = NULL,
                            alpha = 1,
                            mixed = FALSE) {


  p_pert <- p_pert_fun(pars_perturbed = pars_perturbed,
                       pars = pars)

  if(is.null(ic_list)) p_ic_list <- list(P(structure(names(pars), names = names(pars)))) else p_ic_list <- p_ic_list_fun(ic_list, mixed = mixed)


  out <- lapply(seq_along(p_ic_list), function(ic) {


    # Changed parameter trafo function
    cat("Condition ", ic, " of ", length(p_ic_list), "\n")
    p_fun <- (p_pert * p_ic_list[[ic]])
    # p_fun(pars)
    # Perturbation Data generation,
    perturbation_prediction <-  (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)
    r_kept <- r_kept_fun(pars_opt = pars_opt,
                         perturbation_prediction = perturbation_prediction,
                         obs_fun = obs_fun,
                         p_fun = p_fun,
                         pars = pars,
                         alpha = alpha)

    # Fitting
    myfits <- mstrust(obj_alpha,
                      center =  pars_opt+0.5,
                      studyname = "Fits",
                      cores = cores,
                      fits = nfits,
                      sd = 0.3,
                      perturbation_prediction = perturbation_prediction,
                      r_kept = r_kept,
                      obs_fun = obs_fun,
                      p_fun = p_fun,
                      mypars = pars,
                      fixed = NULL,
                      iterlim = 50)

    # print(myfits)
    # Analysis of fits and output preparation

    # select rows, in which the optimiziation pars are smaller than 1
    best_fit <- try({
      bf <- myfits %>% as.parframe()
      mynames <- bf %>% as.parvec() %>% names()

      if(length(mynames) == 1) {
        myrows <- bf[,mynames] > 0
      } else {
        myrows <- apply(bf[,mynames], 1, function(i) {
          i <- i > 0
          do.call("&", as.list(i)) # connect all alpha_j with "&"
        })
      }

      if(all(!myrows)) myrows <- 1 # If there exists no fit with alpha > 0

      bf[myrows, ]%>% as.parvec()
      })

    # If no fits are in parameter range <1, take the best fit
    if(inherits(best_fit, "try-error"))
      best_fit <- try({myfits %>% as.parframe() %>% as.parvec()})

    # if there are still problems, return just the fitlist and the condition
    if(inherits(best_fit, "try-error")) {return(list(all_fits = myfits %>% as.parframe,
                                                     best_fit = 0,
                                                     r_0 = matrix(0),
                                                     r_alpha = matrix(0),
                                                     r_kept = matrix(0),
                                                     r_best_fit = matrix(0),
                                                     steady_states = matrix(0),
                                                     prediction = matrix(0),
                                                     condition = names(p_ic_list)[ic],
                                                     pars_perturbed = pars_perturbed))
    } else

    # r_0 before optimization
    r_0 <- R_fun(pars_opt = structure(rep(0, length(pars_opt)), names = names(pars_opt)),
                 perturbation_prediction = perturbation_prediction,
                 obs_fun = obs_fun,
                 p_fun = p_fun,
                 pars = pars) %>% local_response_matrix_eq10()
    # r_0 before optimization with different observation function
    r_alpha <- R_fun(pars_opt = structure(rep(alpha, length(pars_opt)), names = names(pars_opt)),
                                      perturbation_prediction = perturbation_prediction,
                                      p_fun = p_fun,
                                      pars = pars) %>% local_response_matrix_eq10()


    r_best_fit <- R_fun(pars_opt = best_fit,
                        perturbation_prediction = perturbation_prediction,
                        obs_fun = obs_fun,
                        p_fun = p_fun,
                        pars = pars) %>% local_response_matrix_eq10()

    steady_states <- (xs*p_ic_list[[ic]])(times = c(0,Inf), pars = pars, deriv = F)

    prediction <- (x*p_ic_list[[ic]])(times = seq(0,steady_states[[1]][1], length.out = 50), pars = pars, deriv = F)

    return(list(all_fits = myfits %>% as.parframe,
                best_fit = best_fit,
                r_0 = r_0,
                r_alpha = r_alpha,
                r_kept = r_kept,
                r_best_fit = r_best_fit,
                steady_states = steady_states,
                prediction = prediction,
                condition = names(p_ic_list)[ic],
                pars_perturbed = pars_perturbed))
  })

  names(out) <- names(p_ic_list)

  return(out)
}








#' Run simulations with noise
#'
#' @param ic_list
#' @param nfits
#' @param cores
#' @param pars_opt
#' @param obs_fun
#' @param pars
#' @param pars_perturbed
#' @param fixed
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
run_simulations_with_noise <- function(ic_list = NULL,
                            # nfits= 10,
                            n_samples = 100,
                            sdlog = 0.01,
                            cores = 3,
                            pars_opt = pars_opt_0,
                            obs_fun = g,
                            pars = pars_0,
                            pars_perturbed = pars_perturbed_0,
                            fixed = NULL,
                            alpha = 1,
                            mixed = FALSE) {


  p_pert <- p_pert_fun(pars_perturbed = pars_perturbed,
                       pars = pars)

  if(is.null(ic_list)) p_ic_list <- list(P(structure(names(pars), names = names(pars)))) else p_ic_list <- p_ic_list_fun(ic_list, mixed = mixed)


  out <- lapply(seq_along(p_ic_list), function(ic) {

    # Changed parameter trafo function
    cat("Condition ", ic, " of ", length(p_ic_list), "\n")
    p_fun <- (p_pert * p_ic_list[[ic]])
    # p_fun(pars)
    # Perturbation Data generation,
    perturbation_prediction <-  (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)


    # Add Log-Normal noise to it and simulate many many times,
    noisysimulation <- mclapply(1:n_samples, function(i) {

      perturbation_prediction <- lapply(perturbation_prediction, function(j) j*rlnorm(length(j), sdlog = sdlog))

      r_kept <- r_kept_fun(pars_opt = pars_opt,
                           perturbation_prediction = perturbation_prediction,
                           obs_fun = obs_fun,
                           p_fun = p_fun,
                           pars = pars,
                           alpha = alpha)

      # Fitting
      myfits <- mstrust(obj_alpha,
                        center =  pars_opt+0.5,
                        studyname = "Fits",
                        cores = cores,
                        fits = 1,
                        sd = 0.01,
                        perturbation_prediction = perturbation_prediction,
                        r_kept = r_kept,
                        obs_fun = obs_fun,
                        p_fun = p_fun,
                        mypars = pars,
                        fixed = NULL,
                        iterlim = 50)

      # print(myfits)
      # Analysis of fits and output preparation

      # select rows, in which the optimiziation pars are smaller than 1
      best_fit <- try({
        bf <- myfits %>% as.parframe()
        mynames <- bf %>% as.parvec() %>% names()

        if(length(mynames) == 1) {
          myrows <- bf[,mynames] > 0
        } else {
          myrows <- apply(bf[,mynames], 1, function(i) {
            i <- i > 0
            do.call("&", as.list(i)) # connect all alpha_j with "&"
          })
        }

        if(all(!myrows)) myrows <- 1 # If there exists no fit with alpha > 0

        bf[myrows, ]%>% as.parvec()
      })

      # If no fits are in parameter range <1, take the best fit
      if(inherits(best_fit, "try-error"))
        best_fit <- try({myfits %>% as.parframe() %>% as.parvec()})

      # if there are still problems, return just the fitlist and the condition
      if(inherits(best_fit, "try-error")) {return(list(best_fit = pars_opt,
                                                       r_0 = matrix(0, nrow = length(pars_perturbed), ncol = length(pars_perturbed)),
                                                       r_kept = matrix(FALSE, nrow = length(pars_perturbed), ncol = length(pars_perturbed)),
                                                       r_best_fit = matrix(0, nrow = length(pars_perturbed), ncol = length(pars_perturbed))
                                                       ))
      }

      # r_0 before optimization
      r_0 <- R_fun(pars_opt = structure(rep(0, length(pars_opt)), names = names(pars_opt)),
                   perturbation_prediction = perturbation_prediction,
                   obs_fun = obs_fun,
                   p_fun = p_fun,
                   pars = pars) %>% local_response_matrix_eq10()

      r_best_fit <- R_fun(pars_opt = best_fit,
                          perturbation_prediction = perturbation_prediction,
                          obs_fun = obs_fun,
                          p_fun = p_fun,
                          pars = pars) %>% local_response_matrix_eq10()


      return(list(best_fit = best_fit,
                  r_0 = r_0,
                  r_kept = r_kept,
                  r_best_fit = r_best_fit))
    }, mc.cores = cores)


    # Bring all objects in data.frame-format so that the mc simulations can be plotted in histograms.
    noisysimulation <- lapply(seq_along(noisysimulation[[1]]), function(i) lapply(seq_along(noisysimulation), function(j) noisysimulation[[j]][[i]] %>% matrix(nrow = 1) %>% as.data.frame)  %>%  do.call(rbind,.))
    # naming of the columns of the data.frames
    names(noisysimulation[[1]]) <- names(pars_opt)
    for(i in 2:4) names(noisysimulation[[i]]) <- paste0("r", outer(1:(length(pars_perturbed)),1:length(pars_perturbed), FUN = paste0))


    steady_states <- (xs*p_ic_list[[ic]])(times = c(0,Inf), pars = pars, deriv = F)
    prediction <- (x*p_ic_list[[ic]])(times = seq(0,steady_states[[1]][1], length.out = 50), pars = pars, deriv = F)



    return(list(best_fits = noisysimulation[[1]],
                r_0 = noisysimulation[[2]],
                r_kept = noisysimulation[[3]],
                r_best_fit = noisysimulation[[4]],
                steady_states = steady_states,
                prediction = prediction,
                condition = names(p_ic_list)[ic],
                pars_perturbed = pars_perturbed,
                sdlog = sdlog))
  })

  names(out) <- names(p_ic_list)

  return(out)
}
