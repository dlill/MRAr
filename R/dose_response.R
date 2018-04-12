#' Compute dose-response curves
#'
#' Compute the steady state dose-response wrt a parameter.
#' Return the dose response in a tidy data format.
#'
#' @param which_par Character of length 1: Which par to perturb
#' @param dosages Numeric including the dosages eq(1,50,by = 1)
#' @param pars dynamical parameters
#' @param predfun prediction function, can be concatenated with parameter trafos and observation functions
#'
#' @return
#' @export
#'
#' @examples
dose_response_fun <- function(which_par, dosages, pars = pars_0, predfun = (xs*p_log)) {

  d_r <- mclapply(dosages, function(dose) {

    mypars <- pars
    mypars[which_par] = dose

    conditions <- attr(predfun, "conditions")
    if(is.null(conditions)) conditions <- "1"

    my_steady_states <-  data.frame(condition = conditions,
                                    dose = dose,
                                    which_par = which_par,
                                    do.call(rbind, predfun(times = c(0,Inf),  pars = mypars, deriv = F)))
  }, mc.cores = 4) %>% do.call(rbind, .)

  # exclude time
  d_r <- d_r[,!(colnames(d_r) %in% "time")]

  observed_names <- colnames(d_r)[-c(1:3)]

  d_r <- gather(data = d_r, key = "name", value =  "value", rlang::UQS(rlang::syms(observed_names)), factor_key = T)
  return(d_r)

}

#' Plot dose_response curves
#'
#' @param d_r
#'
#' @return
#' @export
#'
#' @examples
plot_dose_response <- function(d_r) {
  ggplot(d_r, aes(x = log(which_par), y = log(value), color = condition)) + facet_wrap(facets = ~name) + geom_line()+theme_dMod()
}



#' Combine a list of dose_response data.frames with different pars_opt combinations
#'
#' @param dr_list
#'
#' @return
#' @export
#'
#' @examples
combine_dr_list <- function(dr_list) {
  dr_list_names <- dr_list %>%
    lapply(names) %>%
    do.call(c,.) %>%
    unique()

  dr_list %>%
    seq_along() %>%
    lapply(function(i) {
      not_in_i <- dr_list_names[! dr_list_names %in% names(dr_list[[i]])]
      if(length(not_in_i)==0) {
        return(data.frame(dr_list[[i]], par_opt_setting = i))
      }
      else {
        return(data.frame(dr_list[[i]],
                          matrix(rep(NA, length(not_in_i)*nrow(dr_list[[i]])),
                                 ncol=length(not_in_i)) %>% set_colnames(not_in_i),
                          par_opt_setting = i))
      }
    }) %>%
    do.call(rbind,.)
}



#' Compute dose responses for one parameter
#'
#' Scan free parameters separately
#'
#' @param which_par Which dynamical parameter shall be scanned?
#' @param dosages Scan range of this parameter
#' @param alpha_scan_range Which values of free parameters a shall be scanned?
#' @param pars_opt Free parameters a_ij to control the influence of the complexes
#' @param obs_fun Observation function which returns the communicating species
#' @param p_fun Parameter transformation which returns different conditions
#' @param pars Parameters of the ODE model
#'
#' @return
#' @export
#'
#' @examples
r_opt_dose_response <- function(which_par,
                                dosages,
                                alpha_scan_range = seq(-9,9,by = 2),
                                pars_opt = pars_opt_0,
                                obs_fun = g,
                                p_fun = (p_log * p_pert),
                                pars = pars_0,
                                pars_opt_default_value = 0,
                                cores = detectFreeCores()) {


  r_names_0 <- paste0("r", outer(as.character(1:length(modules0)),as.character(1:length(modules0)), paste0))



  out <- lapply(dosages, function(dose) {

    pars[which_par] <- structure(dose, names = which_par)

    perturbation_prediction <-  (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)

    # r_alpha to plot and to compute r_kept


    r_alpha_df <- sapply(alpha_scan_range, function(a) {
      # print(a)
      myr <- r_alpha_fun(pars_opt = structure(rep(a, length(pars_opt)), names = names(pars_opt)), # center around logalpha = 0
                         perturbation_prediction = perturbation_prediction,
                         obs_fun = obs_fun,
                         p_fun = p_fun,
                         pars = pars)
      return(myr)
    }) %>%
      t %>%
      structure(. , dimnames = list(NULL, r_names_0)) %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      tidyr::gather(key = "r_element", value = "r_alpha", rlang::UQS(rlang::syms(r_names_0)))


    # r_0
    r_0 <- r_alpha_fun(pars_opt = structure(rep(log(0+.Machine$double.eps), length(pars_opt)), names = names(pars_opt)),
                       perturbation_prediction = perturbation_prediction,
                       obs_fun = obs_fun,
                       p_fun = p_fun,
                       pars = pars)

    # r_kept
    r_kept_df <- sapply(alpha_scan_range, function(a) {
      myr <- r_alpha_df %>%
        dplyr::filter(alpha == a) %>%
        .[["r_alpha"]] %>%
        matrix(nrow = length(modules0))
      # print(myr)
      return(r_kept_fun2(r_0,myr))
    }) %>%
      {.+(seq(0,0.05,length.out = length(r_names_0)))} %>%  # könnte man evtl auch mit "jitter" lösen
      t %>%
      structure(. , dimnames = list(NULL, r_names_0)) %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      gather(key = "r_element", value = "r_kept", rlang::UQS(rlang::syms(r_names_0)))

    # r_kept <- r_kept_df  %>%  # könnte man evtl auch mit "jitter" lösen
    #   extract2(2) %>%
    #   is_greater_than(0.5) %>%
    #   matrix(ncol = (obs_fun %>% attr("ma") %>% extract2(1) %>% attr("eq") %>% length)) %>%
    #   t()

    # determine the alpha values at which a sign flip occured.
    flipped <- r_kept_df %>%
      tidyr::spread(r_element, r_kept) %>%
      select(-alpha) %>%
      apply(2, diff) %>%
      rbind(0,.) %>%
      apply(1, function(i) as.logical(i) %>% any)

    if(any(r_kept_df[r_kept_df[["alpha"]] == min(alpha_scan_range), "r_kept"]>0.5)) {
      flipped[1] <- TRUE
    }
    # when no sign flip occurs between two alphas, one can be sure that the alpha_opt, that is found will be found again.
    # No fitting is necessary in this case. Initialize an object (booh) to intermediately store the r_opt to be able to access it within the loop over the alphas
    r_opt_prev <- r_0
    best_fit_prev <- structure(rep(min(alpha_scan_range), length(pars_opt)), names = names(pars_opt))

    # r_opt
    r_opt_df <- sapply(seq_along(alpha_scan_range), function(i) {

      if(flipped[i]) { # if a sign flip occured in any of the matrix elements, optimize

        r_kept <- r_kept_df %>%
          dplyr::filter(alpha == alpha_scan_range[i] ) %>%
          .[["r_kept"]] %>%
          {.>0.5} %>%
          matrix(nrow = length(modules0))

        pars_opt[!which_pars_opt(pars_opt,r_kept,modules)] <- pars_opt_default_value
        pars[names(pars_opt)] <- pars_opt


        assign("obj", obj_alpha, envir = .GlobalEnv)


        myfits <- mstrust(obj,
                          center =  structure(rep(-1, sum(which_pars_opt(pars_opt,r_kept,modules))), names = names(pars_opt[which_pars_opt(pars_opt,r_kept,modules)])), # center around log(pars_opt) = -1
                          studyname = "pars_controlled",
                          resultPath = ".pars_controlled/",
                          cores = round(cores),
                          fits = round(cores),
                          sd = 2,
                          perturbation_prediction = perturbation_prediction,
                          r_kept = r_kept,
                          obs_fun = obs_fun,
                          p_fun = p_fun,
                          mypars = pars,
                          iterlim = 50)

        mybest_fit <- try(myfits %>% as.parframe %>% as.parvec)

        #

        if(inherits(mybest_fit, "try-error")) {
          r_opt <- r_opt_prev
          mybest_fit <- best_fit_prev
        } else {
          # "prev" = previous wrt alpha, see further up, where r_opt_prev is constructed
          r_opt <- r_alpha_fun(pars_opt = mybest_fit,
                               perturbation_prediction = perturbation_prediction,
                               obs_fun = obs_fun,
                               p_fun = p_fun,
                               pars = pars)

          # save r_opt and bestfit outside of this environment
          r_opt_prev <<- r_opt
          mybest_fit <- c(mybest_fit, pars_opt[!which_pars_opt(pars_opt,r_kept,modules)])[names(pars_opt)]
          best_fit_prev <<- mybest_fit
        }
      } else {
        r_opt <- r_opt_prev
        mybest_fit <- best_fit_prev
        # print(i)
      }



      r_opt <- c(r_opt, mybest_fit) %>% set_names(c(r_names_0, names(pars_opt)))

      return(r_opt)
    }) %>%
      t %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      tidyr::gather(key = "r_element", value = "r_opt", rlang::UQS(rlang::syms(r_names_0)))


    steady_states <- (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)
    if(any(steady_states[[1]][-1] > 95)) warning("Warning in dose level ", dose, ". Saturation might have been reached. \n", steady_states[[1]][-1], "\n")


    r_df <- plyr::join(r_alpha_df, r_kept_df) %>% plyr::join(r_opt_df)
    r_df <- data.frame(r_df, which_par = which_par, dose = dose)

    return(r_df)
  })


  out <- out %>%
    do.call(rbind,.) %>%
    gather(key = "matrix", value  = "value", r_alpha, r_kept, r_opt) %>%
    # gather_(key_col = "par_opt", value_col  = "par_opt_value", gather_cols = names(pars_opt))
    gather(key = "par_opt", value  = "par_opt_value", rlang::UQS(rlang::syms(names(pars_opt))))

}



#' Compute dose responses for one parameter, scan pars_opt separately
#'
#' Scan each free parameter a_ij individually while keeping others at fixed values.
#' Optimize whenever a sign flip occurs
#'
#' @param which_par Which dynamical parameter shall be scanned?
#' @param dosages Scan range of this parameter
#' @param alpha_scan_range Which values of free parameters a shall be scanned?
#' @param pars_opt Free parameters a_ij to control the influence of the complexes
#' @param obs_fun Observation function which returns the communicating species
#' @param p_fun Parameter transformation which returns different conditions
#' @param pars Parameters of the ODE model
#'
#' @return
#'
#' @examples
r_opt_dose_response_separate_scanning <- function(which_par,
                                                  dosages,
                                                  alpha_scan_range = seq(-9,9,by = 2),
                                                  which_par_opt = "loga",
                                                  pars_opt = pars_opt_0,
                                                  obs_fun = g,
                                                  p_fun = (p_log * p_pert),
                                                  pars = pars_0) {


  r_names_0 <- paste0("r", outer(as.character(1:length(modules0)),as.character(1:length(modules0)), paste0))


  out <- lapply(dosages, function(dose) {

    pars[which_par] <- structure(dose, names = which_par)

    perturbation_prediction <-  (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)

    # r_alpha to plot and to compute r_kept

    r_alpha_df <- sapply(alpha_scan_range, function(a) {
      # print(a)
      mypars_opt <- pars_opt
      mypars_opt[which_par_opt] <- a
      myr <- r_alpha_fun(pars_opt = mypars_opt, # center around logalpha = 0
                         perturbation_prediction = perturbation_prediction,
                         obs_fun = obs_fun,
                         p_fun = p_fun,
                         pars = pars)
      return(myr)
    }) %>%
      t %>%
      structure(. , dimnames = list(NULL, r_names_0)) %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      tidyr::gather(key = "r_element", value = "r_alpha", rlang::UQS(rlang::syms(r_names_0)))


    # r_0
    r_0 <- r_alpha_fun(pars_opt = structure(rep(log(0+.Machine$double.eps), length(pars_opt)), names = names(pars_opt)),
                       perturbation_prediction = perturbation_prediction,
                       obs_fun = obs_fun,
                       p_fun = p_fun,
                       pars = pars)

    # r_kept
    r_kept_df <- sapply(alpha_scan_range, function(a) {
      myr <- r_alpha_df %>%
        dplyr::filter(alpha == a) %>%
        .[["r_alpha"]] %>%
        matrix(nrow = length(modules0))
      # print(myr)
      return(r_kept_fun2(r_0,myr))
    }) %>%
      {.+(seq(0,0.05,length.out = length(r_names_0)))} %>%  # könnte man evtl auch mit "jitter" lösen
      t %>%
      structure(. , dimnames = list(NULL, r_names_0)) %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      tidyr::gather(key = "r_element", value = "r_kept", rlang::UQS(rlang::syms(r_names_0)))

    # r_kept <- r_kept_df  %>%  # könnte man evtl auch mit "jitter" lösen
    #   extract2(2) %>%
    #   is_greater_than(0.5) %>%
    #   matrix(ncol = (obs_fun %>% attr("ma") %>% extract2(1) %>% attr("eq") %>% length)) %>%
    #   t()

    # determine the alpha values at which a sign flip occured.
    flipped <- r_kept_df %>%
      tidyr::spread(r_element, r_kept) %>%
      select(-alpha) %>%
      apply(2, diff) %>%
      rbind(0,.) %>%
      apply(1, function(i) as.logical(i) %>% any)

    if(any(r_kept_df[r_kept_df[["alpha"]] == min(alpha_scan_range), "r_kept"]>0.5)) {
      flipped[1] <- TRUE
    }
    # when no sign flip occurs between two alphas, one can be sure that the alpha_opt, that is found will be found again.
    # No fitting is necessary in this case. Initialize an object (booh) to intermediately store the r_opt to be able to access it within the loop over the alphas
    r_opt_prev <- r_0
    best_fit_prev <- structure(rep(min(alpha_scan_range), length(pars_opt)), names = names(pars_opt))

    # r_opt
    r_opt_df <- sapply(seq_along(alpha_scan_range), function(i) {

      if(flipped[i]) { # if a sign flip occured in any of the matrix elements, optimize

        r_kept <- r_kept_df %>%
          dplyr::filter(alpha == alpha_scan_range[i] ) %>%
          .[["r_kept"]] %>%
          {.>0.5} %>%
          matrix(nrow = length(modules0))

        assign("obj", obj_alpha, envir = .GlobalEnv)

        myfits <- mstrust(obj,
                          center =  structure(rep(-1, length(pars_opt)), names = names(pars_opt)), # center around log(pars_opt) = -1
                          studyname = "pars_controlled",
                          resultPath = ".pars_controlled/",
                          cores = 3,
                          fits = 9,
                          sd = 1,
                          perturbation_prediction = perturbation_prediction,
                          r_kept = r_kept,
                          obs_fun = obs_fun,
                          p_fun = p_fun,
                          mypars = pars,
                          iterlim = 50)

        mybest_fit <- try(myfits %>% as.parframe %>% as.parvec)

        #

        if(inherits(mybest_fit, "try-error")) {
          r_opt <- r_opt_prev
          mybest_fit <- best_fit_prev
        } else {
          # "prev" = previous wrt alpha, see further up, where r_opt_prev is constructed
          r_opt <- r_alpha_fun(pars_opt = mybest_fit,
                               perturbation_prediction = perturbation_prediction,
                               obs_fun = obs_fun,
                               p_fun = p_fun,
                               pars = pars)

          # save r_opt and bestfit outside of this environment
          r_opt_prev <<- r_opt
          best_fit_prev <<- mybest_fit
        }
      } else {
        r_opt <- r_opt_prev
        mybest_fit <- best_fit_prev
        # print(i)
      }

      r_opt <- c(r_opt, mybest_fit) %>% set_names(c(r_names_0, names(pars_opt)))

      return(r_opt)
    }) %>%
      t %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      tidyr::gather(key = "r_element", value = "r_opt", rlang::UQS(rlang::syms(r_names_0)))


    steady_states <- (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)
    if(any(steady_states[[1]][-1] > 95)) warning("Warning in dose level ", dose, ". Saturation might have been reached. \n", steady_states[[1]][-1], "\n")


    r_df <- plyr::join(r_alpha_df, r_kept_df) %>% plyr::join(r_opt_df)
    r_df <- data.frame(r_df, which_par = which_par, dose = dose)

    return(r_df)
  }) %>%
    do.call(rbind,.) %>%
    gather(key = "matrix", value  = "value", r_alpha, r_kept, r_opt) %>%
    gather(key = "par_opt", value  = "par_opt_value", rlang::UQS(rlang::syms(names(pars_opt))))

}



#' Title
#'
#' @param which_par
#' @param dosages
#' @param alpha_scan_range
#' @param pars_opt
#' @param obs_fun
#' @param p_fun
#' @param pars
#' @param grid_points_fine
#' @param grid_points_raw
#'
#' @return
#' @export
#'
#' @examples
scan_grid <- function(which_par,
                      dosages,
                      alpha_scan_range = seq(-9,9,by = 2),
                      pars_opt = pars_opt_0,
                      obs_fun = g,
                      p_fun = (p_log * p_pert),
                      pars = pars_0,
                      grid_points_fine = alpha_scan_range,
                      grid_points_raw = c(-9, 0, 2)) {

  grid_raw <- lapply(1:(length(pars_opt)-1), function(x) return(grid_points_raw)) %>%
    do.call(expand.grid,.)


  permutations <- sapply(1:length(pars_opt), function(i) {
    rep(1:length(pars_opt),2)[i:(i+length(pars_opt)-1)]
  })


  apply(permutations, 1, function(i) {
    which_par_opt <- names(pars_opt)[i==1]
    not_which_par_opt <- names(pars_opt)[i!=1]

    one_dim <- apply(grid_raw, 1, function(j) {
      pars_opt[not_which_par_opt] <- j
      mydose <- r_opt_dose_response_separate_scanning(which_par = which_par,
                                                      dosages = dosages,
                                                      alpha_scan_range = alpha_scan_range,
                                                      which_par_opt = which_par_opt,
                                                      pars_opt = pars_opt,
                                                      obs_fun = obs_fun,
                                                      p_fun = p_fun,
                                                      pars = pars)

      pars_opt_df <- data.frame(mydose[["alpha"]], data.frame(t(pars_opt[not_which_par_opt]))) %>%
        set_names(c(which_par_opt, not_which_par_opt))

      mydose %>%
        select(-alpha) %>% cbind(pars_opt_df)
    }) %>% do.call(rbind,.)

  }) %>% do.call(rbind,.)

}

