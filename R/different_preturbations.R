#' Run a set of different perturbations over varying alphas and a parameter which can be varied
#'
#' @param dose_pars dosages eg  setNames(1:10, "x1")
#' @param pars0 the "default" pars for the model
#' @param alpha_pars The default alpha parameters e.g. c(a_tox_Cxy = log(.Machine$double.eps))
#' @param alphas vector of values for the free pars a to compare r_0 against r_a
#' @param perturbations list of sets of perturbations eg list(c(x1 = log(0.9), y1 = log(0.9)))
#' @param xs xs
#' @param p_log p_log
#' @param p_pert p_pert
#' @param g g
#' @param cores detectFreeCores()
#' @param alpha_par_settings  named list. Entries are vectors of the names of aplha_pars you want to use for the algorithm e.g. list(upstream = c("a_tox_Cxy", "a_toy_Cyz"))
#'
#' @return a tibble
#' @export
run_different_perturbations <- function(dose_pars, pars0, alpha_pars, alphas, perturbations, xs, p_log, p_pert, g, cores, alpha_par_settings) {

  map(seq_along(dose_pars), function(dose_par) {
    dose_par <- dose_pars[dose_par]
    pars <- pars0
# browser()
    if (!(names(dose_par) %in% names(pars)))
      stop("Dose parameter not in pars")

    pars[names(dose_par)] <- dose_par

    map(alphas, function(alpha) {

      mclapply(perturbations, function(myperturbation) {
        p_pert <- p_pert_fun(myperturbation, pars = pars, modelname = NULL)

        perturbation_prediction <- (xs*p_log*p_pert)(c(0,Inf), pars, deriv = F)

        try({r_0 <- R_fun(pars_opt = alpha_pars,
                          perturbation_prediction = perturbation_prediction,
                          obs_fun = g,
                          p_fun = (p_log*p_pert),
                          pars = pars) %>% local_response_matrix_eq10()})
        if(inherits(r_0, "try-error")) return(NULL)

        algo <- function(which_alpha_pars, prefix = NULL) {
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

        algo_results <- imap(alpha_par_settings, ~try(algo(.x, .y))) %>% bind_cols

        bind_cols(tibble(r_0 = list(r_0 %>% round(2)),
                         pars_perturbed = list(myperturbation),
                         alpha = alpha,
                         dose_name = names(dose_par),
                         dose_par = dose_par),
                         algo_results)

      }, mc.cores = cores) %>%
        bind_rows
    }) %>%
      bind_rows
  }) %>%
    bind_rows
}



#' clustering and print function
#'
#' @param results_unnested_matrices result from above
#' @param a alpha, the value at which the local response matrices are compared
#' @param col_start regex like "up_" or "both_"
#' @param nclust nclusters
#' @param which_clust print perturbations of which cluster?
#' @export
cluster_matrices_and_print_perturbations <- function(results_unnested_matrices, a, col_start, nclust, which_clust) {
  kclusts <- results_unnested_matrices %>%
    filter(near(alpha, a)) %>%
    select(starts_with(col_start)) %>%
    {.} %T>% {set.seed(1)} %>%
    kmeans(nclust) %>%
    {.} %T>% print(digits = 2) %>%
    {.}

  kclusts$centers %>%
    round(2) %>%
    apply(1, . %>% matrix( nrow = 3) %>%
    {.} %T>% assign(".mat", ., pos = .GlobalEnv) %>%
      print %>% {NULL})

  augmented <- broom::augment(kclusts, results_unnested_matrices %>% filter(near(alpha, a)))

  augmented %>%
    filter(.cluster %in% which_clust) %>%
    select(X,Y,Z) %>% as.matrix %>% print

  out <- augmented
}

#' Get the relevant matrices for a perturbation and alpha-setting
#'
#' @param results as is from runbg
#' @param a alpha value
#' @param pp names of perturbed parameters
#'
#' @return list of matrices
#' @export
get_r0_r1_ropt_rkept <- function(results, a, pp) {
  myresults <- filter(results, alpha == a, map_lgl(pars_perturbed, . %>% names %>% identical(pp)))
  list(perturbed_parameters = t(pp),
       r_0 = myresults$r_0[[1]],
       r_1 =  myresults$r_1upstream[[1]] %>% round(2),
       r_opt = myresults$r_optupstream[[1]] %>% round(2),
       r_kept = myresults$r_keptupstream[[1]])
}
