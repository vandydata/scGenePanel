
# Test Shiny app using local functions (after revert)

# Load required libraries
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(magrittr)
library(flextable)
library(scales)
library(grid)

# Source all local R functions
source("R/validation.R")
source("R/colors.R")
source("R/panel_umap.R")
source("R/panel_violin.R")
source("R/panel_cellfreqtable.R")
source("R/shiny_ui.R")
source("R/shiny_server.R")
source("R/run_shiny.R")

cat("=== Testing Shiny App from Local Functions ===\n")
cat("✅ All R functions sourced\n")
cat("✅ Libraries loaded\n\n")

cat("**Starting the app...**\n")
cat("This will launch the Shiny app using your local functions\n")
cat("without needing to install the package.\n\n")

# Run the app
scGenePanelShiny()
