test_that("scGenePanel works with valid Seurat input and writes to output_dir", {
  output_dir <- tempfile("scGenePanel-test-")
  dir.create(output_dir)

  result <- scGenePanel(
    object = seurat_data,
    gene = "INS",
    meta_group = "Age",
    cell_type_name = "Beta",
    cell_type_colname = "CellTypes",
    col_palette = "Dark2",
    group_order = c("14y", "39y", "66y"),
    output_dir = output_dir
  )

  expect_s3_class(result, "ggplot")
  expect_true(file.exists(file.path(
    output_dir, "scGenePanel__INS_Beta_Age_umap.png"
  )))
  expect_true(file.exists(file.path(
    output_dir, "scGenePanel__INS_Beta_Age_umap.pdf"
  )))
})

test_that("scGenePanel works with SingleCellExperiment input", {
  output_dir <- tempfile("scGenePanel-test-")
  dir.create(output_dir)

  result <- scGenePanel(
    object = sce_data,
    gene = "INS",
    meta_group = "Age",
    cell_type_name = "Beta",
    cell_type_colname = "CellTypes",
    output_dir = output_dir
  )

  expect_s3_class(result, "ggplot")
})

test_that("scGenePanel validates inputs", {
  expect_error(
    scGenePanel(seurat_data, "NONEXISTENT_GENE", "Age", "Beta", "CellTypes"),
    "not found"
  )
  expect_error(
    scGenePanel(seurat_data, "INS", "Age", "Beta", "INVALID_COLUMN"),
    "Cell type column.*not found in metadata"
  )
  expect_error(
    scGenePanel(seurat_data, "INS", "Age", "INVALID_CELL_TYPE", "CellTypes"),
    "Cell type.*not found in column"
  )
  expect_error(
    scGenePanel(seurat_data, "INS", "INVALID_GROUP", "Beta", "CellTypes"),
    "Metadata group.*not found in object metadata"
  )
  expect_error(
    scGenePanel(data.frame(x = 1:10, y = 1:10), "INS", "Age", "Beta", "CellTypes"),
    "not a Seurat or SingleCellExperiment"
  )
  expect_error(scGenePanel(seurat_data), "argument.*missing")
})
