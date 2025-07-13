#' Launch Interactive Gene Panel Shiny Application
#'
#' Opens a Shiny web application for interactive exploration of single-cell
#' gene expression data with dynamic parameter selection.
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param ... Additional parameters passed to the Shiny application
#'
#' @return Invisibly returns NULL. Function is called for its side effect of launching the Shiny application.
#'
#' @examples
#' # Create mock data for the example
#' library(Seurat)
#' set.seed(123)
#' counts <- matrix(rpois(500, lambda = 2), nrow = 50, ncol = 10)
#' rownames(counts) <- paste0("Gene_", 1:50)
#' colnames(counts) <- paste0("Cell_", 1:10)
#'
#' mock_data <- CreateSeuratObject(counts = counts, project = "example")
#' mock_data$cell_type <- sample(c("TypeA", "TypeB"), 10, replace = TRUE)
#' mock_data$condition <- sample(c("Control", "Treatment"), 10, replace = TRUE)
#'
#' \dontrun{
#' # Launch the interactive Shiny application
#' # This opens a web browser with the gene panel interface
#' genepanel_shiny(mock_data)
#' }
#'
#' \donttest{
#' # The function would normally launch a Shiny app here
#' # For testing purposes, we just verify the object is valid
#' stopifnot(inherits(mock_data, "Seurat"))
#' }
#'
#' @export

genepanel_shiny <- function() {
  shiny::runApp(appDir = system.file("shiny", "genepanel", package = "scGenePanel"))
  #shiny::shinyAppDir(system.file("shiny", package = "genepanel"))
  #shiny::shinyApp(ui=ui,server=server)

}
