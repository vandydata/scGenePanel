
#' Launch Interactive Gene Panel Shiny Application
#'
#' `scgenepanel_shiny` launches a Shiny web application for interactive
#' exploration of single-cell gene expression data with dynamic parameter
#' selection.
#'
#' @param object A Seurat or SingleCellExperiment object. If NULL, loads example data
#' @return Invisibly returns NULL. Function is called for its side effect of launching the Shiny application.
#' @examples
#' # Example 1: Launch with example data
#' if (interactive()){
#'   scgenepanel_shiny()
#' }
#'
#' # Example 2: Launch with seurat object
#' data <- readRDS(system.file("extdata", "human_panc_islets.Rds", package = "scGenePanel"))
#' if (interactive()){
#'   scgenepanel_shiny(data)
#' }
#'
#' # Example 3: Launch with SCE object
#' data <- readRDS(system.file("extdata", "sce_object.Rds", package = "scGenePanel"))
#' if (interactive()){
#'   scgenepanel_shiny(my_data)
#' }
#'
#' @export

scgenepanel_shiny <- function(object = NULL) {

  # Capture the filename of the object
  object_name <- deparse(substitute(object))

  # Load/prepare data
  if (!is.null(object)) {
    # Validate the object
    if (!inherits(object, c("Seurat", "SingleCellExperiment"))) {
      stop("object must be a Seurat or SingleCellExperiment object")
    }
    # Convert to Seurat if needed
    data <- .make_seurat(object)

    # Store the object name for use in the app
    attr(data, "object_name") <- object_name
    message(paste("Using provided data object:", object_name))
  } else {
    # Load example data
    data <- readRDS(system.file("extdata", "human_panc_islets.Rds", package = "scGenePanel"))
    attr(data, "object_name") <- "test_data"
    message("Using example pancreatic islet data")
  }

  # Create UI and server
  ui <- create_genepanel_ui(data)
  server <- create_genepanel_server(data)

  # Launch the app
  shiny::shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}
