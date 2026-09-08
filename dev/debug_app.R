
# Debug version of the Shiny app to identify the issue

library(shiny)
library(shinydashboard)
library(shinyjs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)

# Load data and functions
source("R/validation.R")
source("R/colors.R")
source("R/panel_umap.R")
source("R/panel_violin.R")

data <- readRDS("inst/extdata/human_panc_islets.Rds")
cat("Debug app data loaded: ", ncol(data), " cells\n")

# Simple UI
ui <- dashboardPage(
  dashboardHeader(title = "Debug App"),
  dashboardSidebar(
    actionButton("test_umap", "Test UMAP"),
    actionButton("test_violin", "Test Violin"),
    br(), br(),
    textOutput("debug_info")
  ),
  dashboardBody(
    useShinyjs(),
    h3("Debug Output"),
    plotOutput("debug_plot", height = "400px")
  )
)

# Simple server
server <- function(input, output, session) {
  
  output$debug_info <- renderText({
    paste("Data available:", ncol(data), "cells")
  })
  
  observeEvent(input$test_umap, {
    cat("UMAP button clicked\n")
    output$debug_plot <- renderPlot({
      cat("Rendering UMAP plot\n")
      
      tryCatch({
        levels_idents <- unique(data[["Source"]][, 1])
        levels_idents <- as.character(levels_idents)
        
        result <- umap_panel(seurat_obj = data,
                           cell_type_colname = "CellTypes",
                           cell_type_name = "Beta",
                           meta_group = "Source",
                           gene = "INS",
                           levels_idents = levels_idents)
        
        cat("UMAP plot generated successfully\n")
        return(result)
        
      }, error = function(e) {
        cat("UMAP error:", e$message, "\n")
        plot(1, main = paste("Error:", e$message))
      })
    })
  })
  
  observeEvent(input$test_violin, {
    cat("Violin button clicked\n")
    output$debug_plot <- renderPlot({
      cat("Rendering violin plot\n")
      
      tryCatch({
        levels_idents <- unique(data[["Source"]][, 1])
        levels_idents <- as.character(levels_idents)
        
        result <- violin_sig_panel(seurat_obj = data,
                                 cell_type_colname = "CellTypes",
                                 cell_type_name = "Beta",
                                 meta_group = "Source",
                                 gene = "INS",
                                 col_palette = "Set3",
                                 levels_idents = levels_idents)
        
        cat("Violin plot generated successfully\n")
        return(result)
        
      }, error = function(e) {
        cat("Violin error:", e$message, "\n")
        plot(1, main = paste("Error:", e$message))
      })
    })
  })
}

cat("Starting debug app...\n")
shinyApp(ui = ui, server = server)

