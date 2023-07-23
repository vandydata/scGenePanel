#' single cell genepanel shiny app
#' @export
#'
genepanel_shiny <- function() {
  shiny::runApp(appDir = system.file("shiny", "genepanel", package = "scGenePanel"))
  #shiny::shinyAppDir(system.file("shiny", package = "genepanel"))
  #shiny::shinyApp(ui=ui,server=server)

}
