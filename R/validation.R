#' Check and Convert Seurat Object with Automatic Version Updating
#'
#' Convert SingleCellExperiment object to Seurat object and automatically
#' update Seurat objects to current package version for compatibility
#'
#' @param object A Seurat or SingleCellExperiment object
#' @return returns an error message if the object is not among the expected
#'   input object types, or a properly versioned Seurat object
#' @importFrom methods is
#' @importFrom methods slotNames
#' @importFrom utils packageVersion combn
#' @importFrom Seurat UpdateSeuratObject
#' @noRd

.make_seurat <- function(object) {

  # Handle Seurat objects
  if (inherits(object, "Seurat")) {
    # Check if Seurat object needs version updating
    if ("version" %in% slotNames(object) && !is.null(object@version)) {
      object_version <- package_version(object@version)
      current_version <- packageVersion("Seurat")

      # If object version is older than current Seurat version, update it
      if (object_version < current_version) {
        message("Detected Seurat object version ", object_version,
                " (current: ", current_version, "). Updating object...")

        # Check if UpdateSeuratObject function exists
        if (exists("UpdateSeuratObject", where = asNamespace("Seurat"))) {
          tryCatch({
            object <- UpdateSeuratObject(object)
            message("Successfully updated Seurat object to version ", object@version)
          }, error = function(e) {
            warning("Could not update Seurat object: ", e$message,
                   "\nProceeding with original object, but compatibility issues may occur.")
          })
        } else {
          warning("UpdateSeuratObject function not available. ",
                 "Please update your Seurat installation.")
        }
      }
    }
    return(object)

  } else if (is(object, "SingleCellExperiment")) {
    # Check if the count assay exists
    if (dim(SingleCellExperiment::counts(object))[1] == 0 || dim(SingleCellExperiment::counts(object))[2] == 0) {
      stop(paste0(
        "counts assay is empty"
      ))
    }

    # Convert SingleCellExperiment object to Seurat object
    converted_obj <- Seurat::as.Seurat(object, counts = "counts", data = NULL)

    # Add normalized counts if it doesn't exist
    tryCatch({
      if ("logcounts" %in% SingleCellExperiment::assayNames(object)) {
        converted_obj <- Seurat::as.Seurat(object, counts = "counts", data = "logcounts")
      }
    }, error = function(e) {
      # if log counts does not exist, add log normalization
      converted_obj <- Seurat::NormalizeData(converted_obj)
      message("logcounts assay is empty. log normalization added")
    })

    # Add a compatible assay name "RNA"
    converted_obj@assays$RNA <- converted_obj@assays$originalexp

    return(converted_obj)
  }

  stop("entered 'object' is not a Seurat or SingleCellExperiment")

}

#
#' Check if cell type metadata exists
#'
#' Check to see if cell type metadata column name exist in the seurat object,
#'   if not returns an error
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param cell_type_colname The metadata column name that contains the cell
#' identity annotations
#' @return returns an error message if cell_type_colname does not exist
#' @noRd

.is_celltype_colname <- function(
  object, cell_type_colname
) {
  seurat_obj <- .make_seurat(object)
  meta_col_names <- colnames(x = seurat_obj@meta.data)

  if (!cell_type_colname %in% meta_col_names) {
    error_message <- paste0(
      "Cell type column '", cell_type_colname, "' not found in metadata.\n",
      "Available metadata columns: ", paste(meta_col_names, collapse = ", ")
    )
    stop(error_message, call. = FALSE)
  }
}


#' Check if cell type name exists
#'
#' Check to see if cell type  name exist in the seurat object, if not returns an
#'   error
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param cell_type_name The cell type identity to highlight in UMAP
#' @param cell_type_colname The metadata column name that contains the cell
#'   identity annotations
#' @return returns an error message if cell type name does not exist
#' @noRd

.is_cell_type_name <- function(
  object, cell_type_name, cell_type_colname
) {
  # Check for missing or empty cell_type_name
  if (missing(cell_type_name) || is.null(cell_type_name) || cell_type_name == "") {
    seurat_obj <- .make_seurat(object)
    available_types <- unique(seurat_obj@meta.data[[cell_type_colname]])
    stop(paste0("Please specify a cell_type_name. Available options: ",
                paste(available_types, collapse = ", ")), call. = FALSE)
  }

  seurat_obj <- .make_seurat(object)
  meta <- seurat_obj@meta.data
  cell_type_column <- meta[, cell_type_colname]

  if (!cell_type_name %in% cell_type_column) {
    available_types <- unique(meta[, cell_type_colname])
    error_message <- paste0(
      "Cell type '", cell_type_name, "' not found in column '", cell_type_colname, "'.\n",
      "Available cell types: ", paste(available_types, collapse = ", ")
    )
    stop(error_message, call. = FALSE)
  }
}

#' Check meta_group
#'
#' Check if 'meta_group' exists
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param meta_group The metadata column name of the variable to split the UMAP,
#' violinplot and cell frequency table by. for example, to split by disease
#' condition
#' @return returns an error message if meta_group does not exist
#' @noRd

.is_meta_group <- function(
  object, meta_group
) {
  seurat_obj <- .make_seurat(object)
  meta_col_names <- colnames(x = seurat_obj@meta.data)

  if (!meta_group %in% meta_col_names) {
    error_message <- paste0(
      "Metadata group '", meta_group, "' not found in object metadata.\n",
      "Available metadata columns: ", paste(meta_col_names, collapse = ", ")
    )
    stop(error_message, call. = FALSE)
  }

  # Additional check: warn if meta_group has only one unique value
  meta_values <- unique(seurat_obj@meta.data[[meta_group]])
  if (length(meta_values) == 1) {
    warning("Meta group '", meta_group, "' has only one unique value (",
            meta_values, "). Plots may not be meaningful.", call. = FALSE)
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
#' @noRd

.is_gene <- function(
  object, gene
) {
  seurat_obj <- .make_seurat(object)
  gene_names <- rownames(seurat_obj)
  if (gene %in% gene_names == "FALSE") {
    stop(paste0("Gene '", gene, "' not found in the object.\nUse rownames(your_object) to see all ", length(gene_names), " available genes."), call. = FALSE)
  }
}

#' Check cell count adequacy
#'
#' Check if there are enough cells in each group for meaningful analysis
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param meta_group The metadata column to group by
#' @param cell_type_colname The cell type column
#' @param min_cells Minimum number of cells required per group (default: 3)
#' @return warning if any groups have too few cells
#' @noRd

.check_cell_counts <- function(
  object, meta_group, cell_type_colname, min_cells = 3
) {
  seurat_obj <- .make_seurat(object)
  meta_data <- seurat_obj@meta.data

  # Create cross-tabulation of cell types and meta groups
  cell_counts <- table(meta_data[[cell_type_colname]], meta_data[[meta_group]])

  # Check for groups with too few cells
  low_count_groups <- which(cell_counts < min_cells, arr.ind = TRUE)

  if (nrow(low_count_groups) > 0) {
    low_groups <- paste(
      rownames(cell_counts)[low_count_groups[,1]],
      "in",
      colnames(cell_counts)[low_count_groups[,2]],
      collapse = ", "
    )
    warning(
      "Some groups have very few cells (< ", min_cells, "): ", low_groups,
      ". Results may not be meaningful.",
      call. = FALSE
    )
  }
}


#' Check if dimension reduction embedding exists
#'
#' Check if the specified dimension reduction embedding exists in the Seurat object
#' @param object A Seurat or SingleCellExperiment object
#' @param dim_red The name of the dimension reduction embedding to check
#' @return returns an error message if the specified dimension reduction does not exist
#' #' @noRd
.check_if_dim_red_exists <- function(
  object, dim_red
) {
  seurat_obj <- .make_seurat(object)

  if (!dim_red %in% names(seurat_obj@reductions)) {
    stop(paste0("Dimension reduction '", dim_red, "' not found in the object.\n",
                "Available reductions: ", paste(names(seurat_obj@reductions), collapse = ", ")), call. = FALSE)
  }
}