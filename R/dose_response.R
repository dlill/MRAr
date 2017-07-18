#' Compute dose-response curves
#'
#' Compute the steady state dose-response of a parameter.
#' Put the dose response out in a tidy data format.
#'
#' @param which_par Character of length 1: Which par to perturb
#' @param dosages Numeric including the dosages eq(1,50,by = 1)
#' @param pars
#' @param predfun
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

  d_r <- gather_(data = d_r, key_col = "name", value_col =  "value", gather_cols = observed_names, factor_key = T)
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
#' Optimize whenever a sign flip occurs
#'
#' @param which_par
#' @param dosages
#' @param alpha_scan_range
#' @param pars_opt
#' @param obs_fun
#' @param p_fun
#' @param pars
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
                                pars = pars_0) {


  r_names_0 <- paste0("r", outer(as.character(1:length(modules0)),as.character(1:length(modules0)), paste0))

  out <- lapply(dosages, function(dose) {


    pars[which_par] <- structure(dose, names = which_par)

    perturbation_prediction <-  (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)

    # r_alpha to plot and to compute r_kept
    r_alpha_df <- sapply(alpha_scan_range, function(a) {
      # print(a)
      myr <- r_alpha_fun(pars_opt = pars_opt - pars_opt + a, # center around logalpha = 0
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
      tidyr::gather_(key = "r_element", value = "r_alpha", gather_cols = r_names_0)


    # r_0
    r_0 <- r_alpha_fun(pars_opt = pars_opt-pars_opt - 16,
                       perturbation_prediction = perturbation_prediction,
                       obs_fun = obs_fun,
                       p_fun = p_fun,
                       pars = pars)

    # r_kept
    r_kept_df <- sapply(alpha_scan_range, function(a) {
      myr <- r_alpha_df %>%
        dplyr::filter(alpha == a) %>%
        extract2("r_alpha") %>%
        matrix(nrow = length(modules0))
      # print(myr)
      return(r_kept_fun2(r_0,myr))
    }) %>%
      add(seq(0,0.05,length.out = length(r_names_0))) %>%  # könnte man evtl auch mit "jitter" lösen
      t %>%
      structure(. , dimnames = list(NULL, r_names_0)) %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      tidyr::gather_(key = "r_element", value = "r_kept", gather_cols = r_names_0)

    r_kept <- r_kept_df  %>%  # könnte man evtl auch mit "jitter" lösen
      extract2(2) %>%
      is_greater_than(0.5) %>%
      matrix(ncol = (obs_fun %>% attr("ma") %>% extract2(1) %>% attr("eq") %>% length)) %>%
      t()

    # determine the alpha values at which a sign flip occured.
    flipped <- r_kept_df %>%
      tidyr::spread(r_element, r_kept) %>%
      select(-alpha) %>%
      apply(2, diff) %>%
      rbind(0,.) %>%
      apply(1, function(i) as.logical(i) %>% any)

    if(any(r_kept_df[["r_kept"]]>0.5)&(!any(flipped))) {
      flipped[1] <- TRUE
    }
    # when no sign flip occurs between two alphas, one can be sure that the alpha_opt, that is found will be found again.
    # No fitting is necessary in this case. Initialize an object (booh) to intermediately store the r_opt to be able to access it within the loop over the alphas
    r_opt_prev <- r_0
    best_fit_prev <- pars_opt

    # r_opt
    r_opt_df <- sapply(seq_along(alpha_scan_range), function(i) {

      if(flipped[i]) { # if a sign flip occured in any of the matrix elements, optimize

        r_kept <- r_kept_df %>%
          dplyr::filter(alpha == alpha_scan_range[i] ) %>%
          extract2("r_kept") %>%
          is_greater_than(0.5) %>%
          matrix(nrow = length(modules0))

        obj <<- obj_alpha

        myfits <- mstrust(obj,
                          center =  pars_opt-pars_opt - 2, # center around log(pars_opt) = 0
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

        mybest_fit <- myfits %>% as.parframe %>% as.parvec

        # "prev" = previous wrt alpha, see further up, where r_opt_prev is constructed
        r_opt <- r_alpha_fun(pars_opt = mybest_fit,
                             perturbation_prediction = perturbation_prediction,
                             obs_fun = obs_fun,
                             p_fun = p_fun,
                             pars = pars)

        # save r_opt and bestfit outside of this environment
        r_opt_prev <<- r_opt
        best_fit_prev <<- mybest_fit

      } else {
        r_opt <- r_opt_prev
        mybest_fit <<- best_fit_prev
        # print(i)
      }

      r_opt <- c(r_opt, mybest_fit) %>% set_names(c(r_names_0, names(pars_opt)))

      return(r_opt)
    }) %>%
      t %>%
      cbind(alpha = alpha_scan_range, .) %>%
      as.data.frame(stringsAsFactors =  F) %>%
      tidyr::gather_(key = "r_element", value = "r_opt", gather_cols = r_names_0)


    steady_states <- (xs*p_fun)(times = c(0,Inf), pars = pars, deriv = F)
    if(any(steady_states[[1]][-1] > 95)) warning("Warning in dose level ", dose, ". Saturation might have been reached. \n", steady_states[[1]][-1], "\n")


    r_df <- plyr::join(r_alpha_df, r_kept_df) %>% plyr::join(r_opt_df)
    r_df <- data.frame(r_df, which_par = which_par, dose = dose)

    return(r_df)
  }) %>%
    do.call(rbind,.) %>%
    gather(key = "matrix", value  = "value", r_alpha, r_kept, r_opt) %>%
    gather_(key_col = "par_opt", value_col  = "par_opt_value", gather_cols = names(pars_opt))

}
