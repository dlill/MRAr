#' Run a set of dirreent perturbations over varying alphas and a parameter which can be varied
#'
#' @param pars0 the "default" pars for the model
#' @param alphas vector of values for the free pars a to compare r_0 against r_a
#' @param perturbations list of sets of perturbations eg list(c(x1 = log(0.9), y1 = log(0.9)))
#' @param xs xs
#' @param p_log p_log
#' @param g g
#' @param cores detectFreeCores()
#' @param p_pert p_pert
#' @param dose_pars dosages eg  setNames(1:10, "x1")
#' @param alpha_par_settings  named list e.g. list(upstream = c("a_tox_Cxy", "a_toy_Cyz"))
#'
#' @return a tibble
#' @export
run_different_perturbations <- function(dose_pars, pars0, alphas, perturbations, xs, p_log, p_pert, g, cores, alpha_par_settings ) {

  map(seq_along(dose_pars), function(dose_par) {
    dose_par <- dose_pars[dose_par]
    pars <- pars0
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

          out <- tibble(r_kept = list(NULL), parframes = list(NULL), r_opt = list(NULL))

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

        algo_results <- imap(alpha_par_settings, ~algo(.x, .y)) %>% bind_cols

        # to_upstream_module <- algo(c("a_tox_Cxy", "a_toy_Cyz", "a_toz_Czx"), "upstream")
        # to_both_modules <- algo(c("a_Cxy", "a_Cyz", "a_Czx"), "both")

        bind_cols(tibble(r_0 = list(r_0 %>% round(2)),
                         pars_perturbed = list(myperturbation),
                         alpha = alpha,
                         dose_name = names(dose_par),
                         dose_par = dose_par),
                  algo_results)
      }, mc.cores = cores) %>%
        do.call(rbind,.)
    }) %>%
      do.call(rbind,.)
  }) %>%
    do.call(rbind,.)
}



#' clustering and print function
#'
#' @param results_unnested_matrices result from above
#' @param a alpha, the value at which the local response matrices are compared
#' @param col_start regex like "up_" or "both_"
#' @param nclust nclusters
#' @param which_clust print perturbations of which cluster?
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
#' @param pars,xs,obs_fun,p_log objects that were used to generate the result
#'
#' @return list of matrices
get_r0_r1_ropt_rkept <- function(results, a, pp, pars, xs, obs_fun, p_log) {
  myresults <- filter(results, alpha == a, map_lgl(pars_perturbed, . %>% names %>% identical(pp)))

  p_pert <- p_pert_fun(myresults$pars_perturbed[[1]], pars, NULL)
  perturbation_prediction <- (xs*p_log*p_pert)(c(0,Inf), pars, deriv = F)
  r_1 <-  r_alpha_fun(pars_opt = myresults$parframes[[1]] %>% as.parvec %>% unclass %>%  `*`(0),
                      pars = pars,
                      perturbation_prediction = perturbation_prediction,
                      p_fun = (p_log*p_pert),
                      obs_fun = g) %>% round(2)

  list(perturbed_parameters = t(pp), r_0 = myresults$r_0[[1]], r_1 = r_1, r_opt = myresults$r_opt[[1]] %>% round(2), r_kept = myresults$r_kept[[1]])
}
