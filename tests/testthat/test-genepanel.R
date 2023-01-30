test_that("scGenePanel works with seurat object input", {
  library(scGenePanel)
  data <- readRDS(system.file("extdata", "seurat_object.Rds", package="scGenePanel"))

  gene <- 'INS'
  meta_group <- 'Age'
  cell_type_name <- 'Beta'
  cell_type_colname <- 'CellTypes'
  col.palette <- 'Dark2'

  create_gene_panel(
    data,
    gene,
    meta_group,
    cell_type_name,
    cell_type_colname,
    col.palette,
    group_order = NULL,
    output_dir = "."
  )
})


test_that("scGenePanel works with sce object input", {
  library(scGenePanel)
  data <- readRDS(system.file("extdata", "sce_object.Rds", package="scGenePanel"))

  gene <- 'INS'
  meta_group <- 'Age'
  cell_type_name <- 'Beta'
  cell_type_colname <- 'CellTypes'
  col.palette <- 'Dark2'

  create_gene_panel(
    data,
    gene,
    meta_group,
    cell_type_name,
    cell_type_colname,
    col.palette,
    group_order = NULL,
    output_dir = "."
  )
})

