test_that("discrete palettes provide usable colors and warn on fallbacks", {
  for (palette in c("Set1", "Set2", "Set3", "Dark2", "Paired")) {
    colors <- discrete_col_palette(num_colors = 3, palette = palette)
    expect_gte(length(colors), 3)
    expect_true(all(grepl("^#", colors)))
  }

  expect_warning(
    colors <- discrete_col_palette(num_colors = 3, palette = "not-a-palette"),
    "unknown palette"
  )
  expect_equal(length(colors), 3)
})

test_that("plot panels expose the expected plot structure", {
  groups <- as.character(unique(seurat_data$Age))

  umap <- scGenePanel:::umap_panel(
    seurat_data, "CellTypes", "Beta", "Age", "INS", "umap", groups
  )
  expect_true(inherits(umap, "patchwork"))

  violin <- scGenePanel:::violin_panel(
    seurat_data, "CellTypes", "Beta", "Age", "INS", "Dark2", groups
  )
  expect_s3_class(violin, "ggplot")
  expect_true(any(vapply(
    violin$layers,
    function(layer) inherits(layer$geom, "GeomViolin"),
    logical(1)
  )))
  expect_equal(violin$labels$y, "Expression of INS")
  expect_equal(violin$labels$fill, "Age")

  violin_sig <- scGenePanel:::violin_sig_panel(
    seurat_data, "CellTypes", "Beta", "Age", "INS", "Dark2", groups
  )
  expect_true(inherits(violin_sig, "patchwork"))
})

test_that("cell-frequency summary values and renderer agree", {
  summary <- scGenePanel:::.cellfreq_summary(
    seurat_data, "CellTypes", "Beta", "Age", "INS"
  )


  expect_named(summary, c(
    "Age", "n_cells", "n_expressing (>0.25)", "% expressed",
    "25%", "50%", "75%", "100%"
  ))

  beta_cells <- seurat_data[, seurat_data$CellTypes == "Beta"]
  expected_counts <- table(beta_cells$Age)
  expect_equal(summary$n_cells, unname(as.integer(expected_counts[summary$Age])))
  expect_true(all(summary$`n_expressing (>0.25)` <= summary$n_cells))
  expect_true(all(summary$`% expressed` >= 0 & summary$`% expressed` <= 100))

  table_grob <- scGenePanel:::cellfreq_panel(
    seurat_data, "CellTypes", "Beta", "Age", "INS", "Dark2"
  )
  expect_true(grid::is.grob(table_grob))
})
