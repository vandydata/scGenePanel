# Load test data
test_data <- readRDS(system.file("extdata", "human_panc_islets.Rds", package="scGenePanel"))

# Basic functionality tests
test_that("create_gene_panel works with valid Seurat input", {
  library(scGenePanel)
  
  result <- create_gene_panel(
    object = test_data,
    gene = 'INS',
    meta_group = 'Age',
    cell_type_name = 'Beta',
    cell_type_colname = 'CellTypes',
    col_palette = 'Dark2',
    output_dir = tempdir()
  )
  
  expect_s3_class(result, 'ggplot')
})

test_that("create_gene_panel works with SingleCellExperiment input", {
  library(scGenePanel)
  sce_data <- readRDS(system.file("extdata", "sce_object.Rds", package="scGenePanel"))
  
  result <- create_gene_panel(
    object = sce_data,
    gene = 'INS',
    meta_group = 'Age',
    cell_type_name = 'Beta',
    cell_type_colname = 'CellTypes',
    col_palette = 'Dark2',
    output_dir = tempdir()
  )
  
  expect_s3_class(result, 'ggplot')
})

# Error handling tests
test_that("create_gene_panel handles missing gene error", {
  expect_error(
    create_gene_panel(
      object = test_data,
      gene = 'NONEXISTENT_GENE',
      meta_group = 'Age',
      cell_type_name = 'Beta',
      cell_type_colname = 'CellTypes'
    ),
    "not found"
  )
})

test_that("create_gene_panel handles invalid cell type column", {
  expect_error(
    create_gene_panel(
      object = test_data,
      gene = 'INS',
      meta_group = 'Age',
      cell_type_name = 'Beta',
      cell_type_colname = 'INVALID_COLUMN'
    ),
    "Cell type column.*not found in metadata"
  )
})

test_that("create_gene_panel handles invalid cell type name", {
  expect_error(
    create_gene_panel(
      object = test_data,
      gene = 'INS',
      meta_group = 'Age',
      cell_type_name = 'INVALID_CELL_TYPE',
      cell_type_colname = 'CellTypes'
    ),
    "Cell type.*not found in column"
  )
})

test_that("create_gene_panel handles invalid meta group", {
  expect_error(
    create_gene_panel(
      object = test_data,
      gene = 'INS',
      meta_group = 'INVALID_GROUP',
      cell_type_name = 'Beta',
      cell_type_colname = 'CellTypes'
    ),
    "Metadata group.*not found in object metadata"
  )
})

test_that("create_gene_panel handles invalid object type", {
  expect_error(
    create_gene_panel(
      object = data.frame(x = 1:10, y = 1:10),  # Invalid object type
      gene = 'INS',
      meta_group = 'Age',
      cell_type_name = 'Beta',
      cell_type_colname = 'CellTypes'
    ),
    "not a Seurat or SingleCellExperiment"
  )
})

# Parameter validation tests
test_that("create_gene_panel validates required parameters", {
  expect_error(
    create_gene_panel(
      object = test_data
      # Missing all other required parameters
    ),
    "argument.*missing"
  )
})

# Edge case tests
test_that("create_gene_panel works with different color palettes", {
  palettes_to_test <- c('Set1', 'Set2', 'Set3', 'Dark2', 'Paired')
  
  for(pal in palettes_to_test) {
    result <- create_gene_panel(
      object = test_data,
      gene = 'INS',
      meta_group = 'Age',
      cell_type_name = 'Beta',
      cell_type_colname = 'CellTypes',
      col_palette = pal,
      output_dir = tempdir()
    )
    expect_s3_class(result, 'ggplot')
  }
})

test_that("create_gene_panel works with different genes", {
  genes_to_test <- c('INS', 'GCG', 'SST')  # Different pancreatic markers
  cell_types <- c('Beta', 'Alpha', 'Delta')
  
  for(i in seq_along(genes_to_test)) {
    result <- create_gene_panel(
      object = test_data,
      gene = genes_to_test[i],
      meta_group = 'Age',
      cell_type_name = cell_types[i],
      cell_type_colname = 'CellTypes',
      output_dir = tempdir()
    )
    expect_s3_class(result, 'ggplot')
  }
})

test_that("create_gene_panel handles custom group ordering", {
  result <- create_gene_panel(
    object = test_data,
    gene = 'INS',
    meta_group = 'Age',
    cell_type_name = 'Beta',
    cell_type_colname = 'CellTypes',
    group_order = c('50y', '39y', '14y'),  # Custom order
    output_dir = tempdir()
  )
  expect_s3_class(result, 'ggplot')
})
