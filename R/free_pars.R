

#' Which pars_opt are associated with an element kept in the optimization procedure
#'
#' @param pars_opt
#' @param r_kept
#' @param modules
#'
#' @return
#' @export
#'
#' @examples
which_pars_opt <- function(pars_opt, r_kept, modules) {

  module_grid <- expand.grid(modules, modules, stringsAsFactors = F)

  out <- lapply(names(pars_opt), function(p_o) {
    sapply(1:nrow(module_grid), function(i) {

      str_detect(p_o, module_grid[i,] %>% unlist %>% stringr::regex(ignore_case = T)) %>% # detect module names in pars_opt_name
        all()

    })
  } ) %>%
    lapply("&", r_kept) %>% # which par_opts are associated with a kept element?
    sapply(any)

  return(out)
}



