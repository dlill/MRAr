#' Run the simulations and add noise to the data
#'
#' @inheritParams run_different_perturbations
#'
#' @param errorpars list of vectors with absolute and relative errors: list(setting1 = c(sabs = 0, srel = 0.1))
#' @param N number of realizations for each setting
#' @param nreplicates number of replicates. normally distributed noise is assumed, replicates are averaged before starteing the algorithm
#'
#' @return tibble
#' @export
#'
#' @seealso [run_different_perturbations()]
run_different_perturbations_noisy_with_replicates <- function(dose_pars, pars0, alpha_pars, alphas, perturbations, xs, p_log, p_pert, g, cores, alpha_par_settings,
                                                              errorpars,
                                                              N = 100,
                                                              nreplicates = 3
) {

  # ---------------------------------------------------------- #
  # helper function ----
  # ---------------------------------------------------------- #

  algo <- function(which_alpha_pars, alpha_pars, alpha, perturbation_prediction, g, p_log, p_pert, pars, r_0, obj_alpha, prefix = NULL) {
    if (any(!(which_alpha_pars %in% names(alpha_pars))))
      stop("some alpha_pars in alpha_pars_settings which don't exist in alpha_pars")

    out <- tibble(r_1 = list(NULL), r_kept = list(NULL), parframes = list(NULL), r_opt = list(NULL))

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
    names(out) <- paste0(names(out), prefix)
    return(out)
  }

  # ----------------------------------------------- #
  # .. 1 dose_par ----
  # have one parameter for which one makes a dose_response
  # ----------------------------------------------- #
  map(seq_along(dose_pars), function(dose_par) {
    dose_par <- dose_pars[dose_par]
    pars <- pars0
    # browser()
    if (!(names(dose_par) %in% names(pars)))
      stop("Dose parameter not in pars")

    pars[names(dose_par)] <- dose_par

    # ----------------------------------------------- #
    # .. 2 perturbations ----
    # different perturbations
    # ----------------------------------------------- #
    map(perturbations, function(myperturbation) {
      p_pert <- p_pert_fun(myperturbation, pars = pars, modelname = NULL)
      perturbation_prediction <- (xs*p_log*p_pert)(c(0,Inf), pars, deriv = F)
      # ----------------------------------------------- #
      # .. 3 alphas ----
      # define different alpha values at which the algorithm is run
      # ----------------------------------------------- #
      map(alphas, function(alpha) {


        # ----------------------------------------------- #
        # .. 4 different errorpar settings ----
        # ----------------------------------------------- #
        algo_results <- map(errorpars, function(err) {
          sabs <- err[1]
          srel <- err[2]

          # ----------------------------------------------- #
          # .. 5 different random realizations of data ----
          # ----------------------------------------------- #
          mclapply(seq_len(N), function(i) {
            set.seed(i)
            perturbation_prediction <- map(perturbation_prediction, function(.x) .x + rnorm(length(.x), 0, sabs + srel * .x))
            replicates <- map(seq_len(nreplicates), function(repl) {
              perturbation_prediction <- map(perturbation_prediction, function(.x) .x + rnorm(length(.x), 0, sabs + srel * .x))
            })

            perturbation_prediction <- replicates %>% transpose() %>% map(~do.call(rbind, .x) %>% colMeans %>% t)



            try({r_0 <- R_fun(pars_opt = alpha_pars,
                              perturbation_prediction = perturbation_prediction,
                              obs_fun = g,
                              p_fun = (p_log*p_pert),
                              pars = pars) %>% local_response_matrix_eq10()})
            if(inherits(r_0, "try-error")) return(NULL)

            # ----------------------------------------------- #
            # .. 6alpha_par_settings ----
            # ----------------------------------------------- #
            algo_results.i <- imap(alpha_par_settings, ~try(algo(which_alpha_pars = .x,
                                                                 alpha_pars=alpha_pars, alpha=alpha,
                                                                 perturbation_prediction=perturbation_prediction,
                                                                 g=g, p_log = p_log,
                                                                 p_pert = p_pert, pars = pars,
                                                                 r_0 = r_0, obj_alpha = obj_alpha, prefix = .y))) %>% bind_cols
            algo_results.i <- bind_cols(tibble(r_0 = list(r_0 %>% round(2)),
                                               pars_perturbed = list(myperturbation),
                                               alpha = alpha,
                                               dose_name = names(dose_par),
                                               dose_par = dose_par,
                                               realization = i, srel = srel, sabs = sabs,
                                               replicates = list(replicates),
                                               perturbation_prediction = list(perturbation_prediction)),
                                        algo_results.i)
            return(algo_results.i)

          }, mc.cores = cores) %>% bind_rows
        }) %>% bind_rows()
      }) %>%
        bind_rows
    }) %>%
      bind_rows
  }) %>%
    bind_rows
}


#' @inherit algo
#' @seealso [algo()]
algo_noisy <- function(pars,
                       which_pars_perturbed,
                       which_alpha_pars,
                       alpha_pars,
                       alpha,
                       g,
                       xs,
                       p_log,

                       srel = 0.1,
                       sabs = 0.1,
                       N = 1000,

                       perturbation_effect = 0.9
) {
  if (any(!(which_alpha_pars %in% names(alpha_pars))))
    stop("some alpha_pars in alpha_pars_settings which don't exist in alpha_pars")

  myperturbation <- rep(log(perturbation_effect), length(which_pars_perturbed)) %>% `names<-`(which_pars_perturbed)
  p_pert <- p_pert_fun(myperturbation, pars = pars, modelname = NULL)
  perturbation_prediction <- (xs*p_log*p_pert)(c(0,Inf), pars, deriv = F)

  out <- map(seq_len(N), function(i) {
    set.seed(i)
    perturbation_prediction <- map(perturbation_prediction, function(.x) .x + rnorm(length(.x), 0, sabs + srel * .x))
    out.i <- tibble(r_0 = list(NULL),
                    r_1 = list(NULL),
                    r_kept = list(NULL),
                    parframes = list(NULL),
                    r_opt = list(NULL),
                    perturbation_predictioni = list(perturbation_prediction))


    r_0 <- R_fun(pars_opt = alpha_pars[which_alpha_pars] * 0 + log(.Machine$double.eps),
                 perturbation_prediction = perturbation_prediction,
                 obs_fun = g,
                 p_fun = (p_log*p_pert),
                 pars = pars) %>% local_response_matrix_eq10()
    out.i$r_0 <- list(r_0 %>% round(2))

    r_1 <- R_fun(pars_opt = alpha_pars[which_alpha_pars] * 0 + alpha,
                 perturbation_prediction = perturbation_prediction,
                 obs_fun = g,
                 p_fun = (p_log*p_pert),
                 pars = pars) %>% local_response_matrix_eq10()
    out.i$r_1 <- list(r_1 %>% round(2))

    r_kept <- r_kept_fun2(r_0, r_1)
    out.i$r_kept <- list(r_kept)

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
    out.i$parframes <- list(myparframe)
    if (inherits(myparframe, "try-error")) {
      out.i$parframes <- list(myfits)
      return(out.i)
    }

    r_opt <- r_alpha_fun(pars_opt = myfits %>% as.parframe() %>% as.parvec() %>% unclass() ,
                         pars = pars,
                         perturbation_prediction = perturbation_prediction,
                         p_fun = (p_log*p_pert),
                         obs_fun = g)

    out.i$r_opt <- list(r_opt %>% round(2))
    return(out.i)
  })
  out <- do.call(rbind, out)
  # out <- bind_cols(out, perturbation_prediction0 = map(seq_len(N),~perturbation_prediction))
  return(out)
}
