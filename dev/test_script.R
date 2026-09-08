library(Seurat)
library(scGenePanel)

devtools::load_all("d:/data/scGenePanel")


data <- readRDS(system.file("extdata", "human_panc_islets.Rds", package="scGenePanel"))

scGenePanel::scGenePanel(
  data,
  gene = "INS",
  meta_group = "Age",
  cell_type_name = "Beta",
  cell_type_colname = "CellTypes",
  group_order = NULL,
  output_dir = "dev/"
)



# for manual testing inside of GenePanel.R
object <- data
gene <- "INS"
meta_group <- "Age"
cell_type_name = "Beta"
cell_type_colname <- "CellTypes"
col_palette <- "Dark2"
group_order <- NULL
output_dir <- "dev/"

