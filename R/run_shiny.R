
#' Launch Interactive Gene Panel Shiny Application
#'
#' Opens a Shiny web application for interactive exploration of single-cell
#' gene expression data with dynamic parameter selection.
#'
#' @param object A Seurat or SingleCellExperiment object. If NULL, loads example data
#' @param ... Additional parameters passed to the Shiny application
#'
#' @return Invisibly returns NULL. Function is called for its side effect of launching the Shiny application.
#'
#' @examples
#' \dontrun{
#' # Launch with example data
#' genepanel_shiny()
#' 
#' # Launch with your own data
#' my_data <- readRDS("path/to/your_data.rds")
#' genepanel_shiny(my_data)
#' }
#'
#' @export

genepanel_shiny <- function(object = NULL) {
  
  # Load/prepare data
  if (!is.null(object)) {
    # Validate the object
    if (!inherits(object, c("Seurat", "SingleCellExperiment"))) {
      stop("object must be a Seurat or SingleCellExperiment object")
    }
    # Convert to Seurat if needed
    data <- .make_seurat(object)
    message("Using provided data object")
  } else {
    # Load example data
    data <- readRDS(system.file("extdata", "human_panc_islets.Rds", package = "scGenePanel"))
    message("Using example pancreatic islet data")
  }
  
  # Create UI and server
  ui <- create_genepanel_ui(data)
  server <- create_genepanel_server(data)
  
  # Launch the app
  shiny::shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}

