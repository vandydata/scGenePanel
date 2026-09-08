compact_test_object <- function(object) {
  metadata <- if (inherits(object, "Seurat")) {
    object[[]]
  } else {
    as.data.frame(SummarizedExperiment::colData(object))
  }

  selected_cells <- rownames(metadata)[metadata$CellTypes %in% c("Alpha", "Beta", "Delta")]
  strata <- interaction(
    metadata[selected_cells, "CellTypes"],
    metadata[selected_cells, "Age"],
    drop = TRUE
  )
  cells <- unname(unlist(lapply(split(selected_cells, strata), head, n = 4)))
  compact_object <- object[, cells]

  if (inherits(compact_object, "Seurat")) {
    compact_object$CellTypes <- droplevels(compact_object$CellTypes)
    compact_object$Age <- droplevels(compact_object$Age)
  } else {
    metadata <- SummarizedExperiment::colData(compact_object)
    metadata$CellTypes <- droplevels(metadata$CellTypes)
    metadata$Age <- droplevels(metadata$Age)
    SummarizedExperiment::colData(compact_object) <- metadata
  }

  compact_object
}

seurat_data <- compact_test_object(readRDS(system.file(
  "extdata", "human_panc_islets.Rds", package = "scGenePanel"
)))

sce_data <- compact_test_object(readRDS(system.file(
  "extdata", "sce_object.Rds", package = "scGenePanel"
)))
