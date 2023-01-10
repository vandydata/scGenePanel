

#' Check Seurat Object
#'
#' Convert SingleCellExperiment object to Seurat object
#'
#' @param object A Seurat or SingleCellExperiment object
#' @return returns an error message if the object is not among the expected input object
#' @importFrom methods is


make_seurat <- function(object) {

    if (class(x = object)[[1]] == "Seurat") {
      return(object)

    } else if (is(object,"SingleCellExperiment")) {
      # Check if the count assay exists
      if(dim(SummarizedExperiment::assay(object, "counts"))[1] == 0 || dim(SummarizedExperiment::assay(object,"counts"))[2] == 0){
        stop(paste0(
          "counts assay is empty",
        ))
    }

    # Convert SingleCellExperiment object to Seurat object
    converted_obj <- Seurat::as.Seurat(object, counts = "counts", data = NULL)

    # Add normalized counts if it doesn't exist
    tryCatch({
      if("logcounts" %in% SummarizedExperiment::assayNames(object)){
        converted_obj <- Seurat::as.Seurat(object, counts = "counts", data = "logcounts")
      }
    }, error = function(e) {
      # if log counts does not exist, add log normalization
      converted_obj <- Seurat::NormalizeData(converted_obj)
      print("logcounts assay is empty. log normalization added")
    })

    # Add a compatible assay name "RNA"
    converted_obj <- suppressWarnings(RenameAssays(object = converted_obj, originalexp = 'RNA'))

    return(converted_obj)

    }

  stop("entered 'object' is not a Seurat or SingleCellExperiment")

}

#
#' Check if cell type metadata exists
#'
#' Check to see if cell type metadata column name exist in the seurat object, if not returns an error
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param cell_type_colname cell type column name in metadata
#' @return returns an error message if cell_type_colname does not exist



Is_celltype_colname <- function(
  object, cell_type_colname
) {
  seurat_obj <- make_seurat(object)
  meta_col_names <- colnames(x = seurat_obj@meta.data)
  if (cell_type_colname %in% meta_col_names == "FALSE") {
    stop("Entered 'cell_type_colname' not found in object metadata")
  }
}




#' Check if cell type name exists
#'
#' Check to see if cell type  name exist in the seurat object, if not returns an error
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param cell_type_name cell type name to highlight
#' @param cell_type_colname The metadata column name for the cell identity annotations
#' @return returns an error message if cell type name does not exist



Is_cell_type_name <- function(
  object, cell_type_name, cell_type_colname
) {
  seurat_obj <- make_seurat(object)
  meta <- seurat_obj@meta.data
  cell_type_column <- meta[,cell_type_colname]
  if (cell_type_name %in% cell_type_column == "FALSE") {
    stop("Entered 'cell_type_name' not found in the 'cell_type_colname' of object metadata")
  }
}


#' Check meta_group
#'
#' Check if 'meta_group' exists
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param meta_group name of group to split visuals by
#' @return returns an error message if meta_group does not exist

Is_meta_group <- function(
  object, meta_group
) {
  seurat_obj <- make_seurat(object)
  meta_col_names <- colnames(x = seurat_obj@meta.data)
  if (meta_group %in% meta_col_names == "FALSE") {
    stop("Entered 'meta_group' not found in object metadata")
  }
}

#
#' Check gene exists
#'
#' Check if gene included in object
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param gene name of gene
#' @return returns an error message if gene is not found in seurat object


Is_gene <- function(
  object, gene
) {
  seurat_obj <- make_seurat(object)
  gene_names <- rownames(seurat_obj)
  if (gene %in% gene_names == "FALSE") {
    stop("Entered 'gene' in object not found")
  }
}

