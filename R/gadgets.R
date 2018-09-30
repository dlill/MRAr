#' Shiny Gadget for quick exploration of Complexes and alpha
#'
#' @param which_pars_opt All possible pars opt that you might want to include (you can kick them out later)
#' @param dMod.frame the dMod.frame with additional columns xs, p_pert, pars and pars_opt
#' @param hypothesis indec to subset the dMod.frame
#'
#' @return NULL
#' @export
ropt_gadget <- function(which_pars_opt, dMod.frame, hypothesis = 1) {

  # ------------ ui ------------
  ui <- miniPage(
    gadgetTitleBar("Explore the algorithm's behaviour depending on complexes and alphas"),
    miniContentPanel(
      # The brush="brush" argument means we can listen for
      # brush events on the plot using input$brush.
      tableOutput("r0"),
      tableOutput("rkept"),
      tableOutput("ropt"),
      selectInput("whichParOpt", "which pars to compare", which_pars_opt, selected = which_pars_opt, multiple = T),
      sliderInput("value", "value", -4,4, -10, step = 0.001),
      actionButton("save", "Save")
    )
  )

  # ------------ Preparation for server --------
  # get objects from dMod.frame
  g <- dMod.frame$g[[hypothesis]]
  p_pert <- dMod.frame$p_pert[[hypothesis]]
  p <- dMod.frame$p[[hypothesis]]
  xs <- dMod.frame$xs[[hypothesis]]
  pars <- dMod.frame$pars[[hypothesis]]
  pars_opt <- dMod.frame$pars_opt[[hypothesis]]


  perturbation_prediction <- (xs*p*p_pert)(c(0,Inf), pars, deriv = F)

  # ------------ Server  ------------
  server <- function(input, output, session) {

    r0 <- R_fun(pars_opt = pars_opt,
                perturbation_prediction = perturbation_prediction,
                obs_fun = g,
                p_fun = (p*p_pert),
                pars = pars) %>% local_response_matrix_eq10()
    output$r0 <- renderTable(r0)

    rkept <- reactive(r_kept_fun(pars_opt = pars_opt[input$whichParOpt],
                                 perturbation_prediction = perturbation_prediction,
                                 obs_fun = g,
                                 p_fun = (p * p_pert),
                                 pars = pars,
                                 alpha =
                                   -log(.Machine$double.eps) + input$value))
    output$rkept <- renderTable(rkept() %>% {mymat <- .;mymat[!mymat] <- mymat[!mymat] %>% str_replace("FALSE", "0"); mymat})

    ropt <- reactive({
      myfits <- mstrust(obj_alpha,
                        center =  structure(rep(0, length(input$whichParOpt)), names = input$whichParOpt),
                        studyname = "Fits",
                        cores = 3,
                        fits = 3,
                        sd = 1,
                        mypars = pars,
                        perturbation_prediction = perturbation_prediction,
                        r_kept = rkept(),
                        p_fun = (p * p_pert),
                        obs_fun = g)

      try(myfits %>% as.parframe()) %>%
      {e <- .; if (inherits(e, "try-error")) print(myfits)}


      out <- r_alpha_fun(pars_opt = myfits %>% as.parframe() %>% as.parvec() %>% unclass() ,
                         pars = pars,
                         perturbation_prediction = perturbation_prediction,
                         p_fun = (p*p_pert),
                         obs_fun = g)
      attr(out, "pars_opt") <- myfits %>% as.parframe() %>% as.parvec() %>% unclass()

      out
    })

    output$ropt <- renderTable(ropt())

    observeEvent(input$save, {
      tpaste0 <- function(...) {
        paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
               , "_" , ...)
      }
      filename <- paste0("~/Promotion/Projects/MRA/Simulations/Models/CascadeComprehensiveFeedbacks/Exploration/", tpaste0("result"))
      out <- rbind(r0 %>% round(2), rkept() %>% round(2),ropt() %>% round(2)) %>% as.data.frame %>%
        rbind(0) %>% rbind(0) %>% rbind(0)
      pars_opt <- ropt() %>% attr("pars_opt") %>% {foo <- .; attr(foo, "alpha") <- input$value; foo} %>% {foo <- .; c(foo, rep(0, nrow(out)-length(foo)))}
      write.csv(cbind(out, pars_opt, names(pars_opt)), file = paste0(filename, ".csv"))
      saveRDS(ropt() %>% attr("pars_opt") %>% {foo <- .; attr(foo, "alpha") <- input$value; foo}, paste0(filename, "_pars.rds"))
    })

    observeEvent(input$done, stopApp(NULL))
  }

  runGadget(ui, server)
}
