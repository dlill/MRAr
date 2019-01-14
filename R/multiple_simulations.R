#' Run the optimization procedure at many conditions (old, not used for results presented in the paper)
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
#'
#' @export
#'
#'
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


  r_names_0 <- paste0("r", outer(as.character(1:length(modules0)),as.character(1:length(modules0)), paste0))

  p_pert <- p_pert_fun(pars_perturbed = pars_perturbed,
                       pars = pars)

  if(is.null(ic_list)) p_ic_list <- list(P(structure(names(pars), names = names(pars)))) else p_ic_list <- p_ic_list_fun(ic_list, mixed = mixed)


  out <- lapply(seq_along(p_ic_list), function(ic) {


    # Changed parameter trafo function
    cat("Condition ", ic, " of ", length(p_ic_list), "\n")
    p_fun <- (p_log*p_pert*p_ic_list[[ic]])
    # p_fun(pars)
    # Perturbation Data generation,
    perturbation_prediction <-  (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)




    alpha_scan_range <- seq(-9,12,length.out = 50)

    # r_alpha to plot and to compute r_kept
    r_alpha_df <- sapply(alpha_scan_range, function(a) {
        # print(a)
        myr <- r_alpha_fun(pars_opt = pars_opt_0 - pars_opt_0 + a, # center around logalpha = 0
                           perturbation_prediction = perturbation_prediction,
                           obs_fun = obs_fun,
                           p_fun = p_fun,
                           pars = pars)
        # print(myr)
        return(myr)
      }) %>%
        t %>%
        structure(. , dimnames = list(NULL, r_names_0)) %>%
        cbind(alpha = alpha_scan_range, .) %>%
        as.data.frame(stringsAsFactors =  F) %>%
        tidyr::gather_(key = "r_element", value = "r_alpha", gather_cols = r_names_0)


    r_0 <- r_alpha_fun(pars_opt = pars_opt_0-pars_opt_0 - 9,
                perturbation_prediction = perturbation_prediction,
                obs_fun = obs_fun,
                p_fun = p_fun,
                pars = pars)

    r_kept_df <- sapply(alpha_scan_range, function(a) {
        myr <- r_alpha_df %>%
          dplyr::filter(alpha == a) %>%
          extract2("r_alpha") %>%
          matrix(nrow = length(modules))
        # print(myr)
        return(r_kept_fun2(r_0,myr))
      }) %>%
        add(seq(0,0.05,length.out = 9)) %>%
        t %>%
        structure(. , dimnames = list(NULL, r_names_0)) %>%
        cbind(alpha = alpha_scan_range, .) %>%
        as.data.frame(stringsAsFactors =  F) %>%
        tidyr::gather_(key = "r_element", value = "r_kept", gather_cols = r_names_0)


    # determine the alpha values at which a sign flip occured.
    flipped <- r_kept_df %>%
      tidyr::spread(r_element, value) %>%
      select(-alpha) %>%
      apply(2, diff) %>%
      rbind(0,.) %>%
      apply(1, function(i) as.logical(i) %>% any)

    # when no sign flip occurs between two alphas, one can be sure that the alpha_opt, that is found will be found again.
    # No fitting is necessary in this case. Initialize an object (booh) to intermediately store the r_opt to be able to access it within the loop over the alphas
    r_opt_prev <- r_0
    best_fits <- list()
    fits <- list()

    r_opt_df <- sapply(seq_along(alpha_scan_range), function(i) {

      if(flipped[i]) { # if a sign flip occured in any of the matrix elements, optimize

        r_kept <- r_kept_df %>%
          dplyr::filter(alpha == alpha_scan_range[i] ) %>%
          extract2("value") %>%
          is_greater_than(0.5) %>%
          matrix(nrow = length(modules))

          myfits <- mstrust(obj_alpha,
                            center =  pars_opt_0-pars_opt_0, # center around log(pars_opt_0) = 0
                            studyname = "pars_controlled",
                            resultPath = ".pars_controlled/",
                            cores = 3,
                            fits = 6,
                            sd = 3,
                            perturbation_prediction = perturbation_prediction,
                            r_kept = r_kept,
                            obs_fun = obs_fun,
                            p_fun = p_fun,
                            mypars = pars,
                            iterlim = 50)

          mybest_fit <- myfits %>% as.parframe %>% as.parvec

          # "prev" = previous wrt alpha, see further up, where r_opt_prev is constructed
          r_opt <- r_alpha_fun(pars_opt = mybest_fit,
                                              perturbation_prediction = perturbation_prediction,
                                              obs_fun = obs_fun,
                                              p_fun = p_fun,
                                              pars = pars)

          # save r_opt and bestfit outside of this environment
          r_opt_prev <<- r_opt
          fits <<- c(fits, myfits %>% as.parframe)
          best_fits <<- c(best_fits, i = mybest_fit)

      } else {
        r_opt <- r_opt_prev
        # print(i)
      }
      return(r_opt)
    }) %>%
      t %>%
      structure(. , dimnames = list(NULL, r_names_0)) %>%
      cbind(alpha = alpha_scan_range, .) %>%
        as.data.frame(stringsAsFactors =  F) %>%
        tidyr::gather_(key = "r_element", value = "r_opt", gather_cols = r_names_0)



    r_df <- join(r_alpha_df, r_kept_df) %>% join(r_opt_df)



    steady_states <- (xs*p_log*p_ic_list[[ic]])(times = c(0,Inf), pars = pars, deriv = F)

    prediction <- (x*p_log*p_ic_list[[ic]])(times = seq(0,steady_states[[1]][1], length.out = 50), pars = pars, deriv = F)

    return(list(fits = fits,
                best_fits = best_fits,
                r_df = r_df,
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
#' @inheritParams run_simulations
#'
#'
#' @export
#'
#'
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
    p_fun <- (p_log * p_pert * p_ic_list[[ic]])
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
      r_0 <- R_fun(pars_opt = structure(rep(-36, length(pars_opt)), names = names(pars_opt)),
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
