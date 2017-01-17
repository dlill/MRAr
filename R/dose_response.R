




#' Compute dose-response curves
#'
#' Compute the steady state dose-response of a parameter.
#' Put the dose response out in a tidy data format.
#'
#' @param par_dose Character of length 1: Which par to perturb
#' @param dosages Numeric including the dosages eq(1,50,by = 1)
#' @param pars
#' @param predfun
#'
#' @return
#' @export
#'
#' @examples
dose_response_fun <- function(par_dose, dosages, pars, predfun = (g*xs)) {

  d_r <- mclapply(dosages, function(i) {
    mypars <- pars
    mypars[par_dose] = i
    conditions <- attr(predfun, "conditions")
    my_steady_states <-  data.frame(condition = conditions,
                                    par_dose = i,
                                    do.call(rbind, predfun(times = c(0,Inf),  pars = mypars, deriv = F)))
  }, mc.cores = 4) %>% do.call(rbind, .)

  # exclude time
  d_r <- d_r[,!(colnames(d_r) %in% "time")]

  dose_response_colnames <- colnames(d_r)[-c(1:2)]

  d_r <- gather_(data = d_r, key_col = "observable", value_col =  "concentration", gather_cols = dose_response_colnames, factor_key = T)
  return(d_r)

}








