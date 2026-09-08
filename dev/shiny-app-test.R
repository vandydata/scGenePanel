# Direct run script to bypass package corruption
# Run this in your scGenePanel directory

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

# Load the example data directly
shiny_data <- readRDS("inst/extdata/human_panc_islets.Rds")

# Make it available globally for the app
assign("shiny_data", shiny_data, envir = .GlobalEnv)

# Print data info
cat("=== Data loaded successfully ===\n")
cat("Cells:", ncol(shiny_data), "\n")
cat("Genes:", nrow(shiny_data), "\n")
cat("Cell types:", paste(unique(shiny_data$CellTypes), collapse = ", "), "\n")
cat("Meta groups:", paste(unique(shiny_data$Source), collapse = ", "), "\n\n")

cat("=== Running Shiny app ===\n")
cat("The app should now work without corruption errors.\n")
cat("Try these parameters:\n")
cat("- Gene: INS\n")
cat("- Cell type name: Beta\n")
cat("- Cell type column: CellTypes\n")
cat("- Meta group: Source\n")
cat("- Color palette: Set3\n\n")

# Run the app directly from the inst directory

runApp("standalone_shiny_app.R", launch.browser = TRUE)

source('test_r_approach.R')

# **After rebuilding package, users will run:**
library(scGenePanel)
scGenePanelShiny()           # Example data
scGenePanelShiny(my_data)    # Custom data


source('test_improved_app.R')
scGenePanelShiny()
