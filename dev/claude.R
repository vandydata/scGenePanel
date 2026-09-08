library(ClaudeR)
claudeAddin()


#

library(Seurat)
data <- readRDS("d:/data/2021-172__KimT-scRNA-Seq-Leukemia/manuscript/Fig4-Entinostat/2022-237-Kim_T-scRNAseq-Leukemia-Entinostat-AntiPD1/04-data_annotated.rds")
data
data@reductions[["X3dtsne_"]] <- NULL
data@reductions[["UMAP3D_"]] <- NULL
data@reductions[["tSNE3D_"]] <- NULL
data@reductions[["tsne3d_"]] <- NULL
data@reductions[["umap3d_"]] <- NULL
data@reductions[["tsne3d2_"]] <- NULL

data <- UpdateSeuratObject(data)

colnames(data@meta.data)


scGenePanel(data,
                  gene="VSIR",
                  meta_group = "condition",
                  cell_type_name = "Monocyte",
                  cell_type_colname = "celltypes")

scGenePanel(data,
                  gene="VSIR",
                  meta_group = "condition",
                  cell_type_name = "Monocyte",
                  cell_type_colname = "celltypes",  # Note the capitals
                  col_palette = "Set3",
                  output_dir = output_dir)



## Shiny app
library(shiny)
library(shinyWidgets)
setwd('D:/data/scGenePanel')  # Make sure you're in the right directory
runApp('inst/shiny/genepanel')



source('test_r_approach.R')


🎉 COMPLETE! Here's what you now have:

📂 **R/ Directory Structure:**
   R/run_shiny.R      ← Main scGenePanelShiny() function
   R/shiny_ui.R       ← UI as a function
   R/shiny_server.R   ← Server as a function
   R/GenePanel.R      ← Your existing functions
   R/panel_*.R        ← Your panel functions

🧪 **To test immediately:**
   source('test_r_approach.R')

📦 **After rebuilding package, users will run:**
   library(scGenePanel)
   scGenePanelShiny()           # Example data
   scGenePanelShiny(my_data)    # Custom data
