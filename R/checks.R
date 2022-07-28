

#' Check Seurat Object
#'
#' Checks if object is a Seurat object
#'
#' @param seurat_object Seurat object name.
#' @return returns an error message if the object is not seurat


Is_Seurat <- function(
  seurat_object
) {
  if (class(x = seurat_object)[[1]] != "Seurat") {
    stop("entered 'seurat_object' is not a seurat object")
  }
}

#
#' Check if cell type metadata exists
#'
#' Check to see if cell type metadata column name exist in the seurat object, if not returns an error
#'
#' @param seurat_object Seurat object name.
#' @param cell_type_colname cell type column name in metadata
#' @return returns an error message if cell_type_colname does not exist



Is_celltype_colname <- function(
  seurat_object
) {
  meta_col_names <- colnames(x = seurat_object@meta.data)
  if (cell_type_colname %in% meta_col_names == "FALSE") {
    stop("Entered 'cell_type_colname' not found in object metadata")
  }
}




#' Check if cell type name exists
#'
#' Check to see if cell type  name exist in the seurat object, if not returns an error
#'
#' @param seurat_object Seurat object name.
#' @param cell_type_name cell type name to highlight
#' @return returns an error message if cell type name does not exist



Is_cell_type_name <- function(
  seurat_object
) {
  meta <- seurat_object@meta.data
  cell_type_column <- meta[,cell_type_colname]
  if (cell_type_name %in% cell_type_column == "FALSE") {
    stop("Entered 'cell_type_name' not found in the 'cell_type_colname' of object metadata")
  }
}


#' Check meta_group
#'
#' Check if 'meta_group' exists
#'
#' @param seurat_object Seurat object name.
#' @param meta_group name of group to split visuals by
#' @return returns an error message if meta_group does not exist

Is_meta_group <- function(
  seurat_object
) {
  if (meta_group %in% meta_col_names == "FALSE") {
    stop("Entered 'meta_group' not found in object metadata")
  }
}

#
#' Check gene exists
#'
#' Check if gene included in object
#'
#' @param seurat_object Seurat object name.
#' @param gene name of gene
#' @return returns an error message if gene is not found in seurat object


Is_gene <- function(
  seurat_object
) {
  gene_names <- rownames(seurat_object)
  if (gene %in% gene_names == "FALSE") {
    stop("Entered 'gene' in object not found")
  }
}

