library(Seurat)

data <- readRDS(system.file("extdata", "human_panc_islets.Rds", package="scGenePanel"))

scGenePanel::create_gene_panel(
  data,
  gene = "INS",
  meta_group = "Age",
  cell_type_name = "Beta",
  cell_type_colname = "CellTypes",
  group_order = NULL,
  output_dir = "."
)
