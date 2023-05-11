

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
      if(dim(SingleCellExperiment::counts(object))[1] == 0 || dim(SingleCellExperiment::counts(object))[2] == 0){
        stop(paste0(
          "counts assay is empty",
        ))
    }

    # Convert SingleCellExperiment object to Seurat object
    converted_obj <- Seurat::as.Seurat(object, counts = "counts", data = NULL)

    # Add normalized counts if it doesn't exist
    tryCatch({
      if("logcounts" %in% SingleCellExperiment::logcounts(object)){
        converted_obj <- Seurat::as.Seurat(object, counts = "counts", data = "logcounts")
      }
    }, error = function(e) {
      # if log counts does not exist, add log normalization
      converted_obj <- Seurat::NormalizeData(converted_obj)
      print("logcounts assay is empty. log normalization added")
    })

    # Add a compatible assay name "RNA"
    #converted_obj <- suppressWarnings(RenameAssays(object = converted_obj, originalexp = 'RNA')) #errors out
    converted_obj@assays$RNA <- converted_obj@assays$originalexp

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
#' @param cell_type_colname The metadata column name that contains the cell identity annotations
#' @return returns an error message if cell_type_colname does not exist



Is_celltype_colname <- function(
  object, cell_type_colname
) {
  seurat_obj <- make_seurat(object)
  meta_col_names <- colnames(x = seurat_obj@meta.data)
  if (cell_type_colname %in% meta_col_names == "FALSE") {
    error_message <- paste("Entered 'cell_type_colname' was not found in metadata. Enter a cell type identity containing column name among the following:",
                           paste(meta_col_names, collapse = ", "))
    stop(error_message)
  }
}




#' Check if cell type name exists
#'
#' Check to see if cell type  name exist in the seurat object, if not returns an error
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param cell_type_name The cell type identity to highlight in UMAP
#' @param cell_type_colname The metadata column name that contains the cell identity annotations
#' @return returns an error message if cell type name does not exist



Is_cell_type_name <- function(
  object, cell_type_name, cell_type_colname
) {
  seurat_obj <- make_seurat(object)
  meta <- seurat_obj@meta.data
  cell_type_column <- meta[,cell_type_colname]
  if (cell_type_name %in% cell_type_column == "FALSE") {
    error_message <- paste("Entered 'cell_type_name' was not found in metadata. Enter a cell type name among the following:",
                           paste(unique(meta[,cell_type_colname]), collapse = ", "))
    stop(error_message)
  }
}


#' Check meta_group
#'
#' Check if 'meta_group' exists
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param meta_group The metadata column name of the variable to split the UMAP, violinplot and cell frequency table by. for example, to split by disease condition
#' @return returns an error message if meta_group does not exist

Is_meta_group <- function(
  object, meta_group
) {
  seurat_obj <- make_seurat(object)
  meta_col_names <- colnames(x = seurat_obj@meta.data)
  if (meta_group %in% meta_col_names == "FALSE") {
    error_message <- paste("Entered 'meta_group' was not found in metadata. Enter a cell identity among the following:",
                           paste(meta_col_names, collapse = ", "))
    stop(error_message)
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

