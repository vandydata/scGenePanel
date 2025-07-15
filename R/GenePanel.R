#' scGenePanel - A method for generating scRNAseq multi-panel gene expression
#' visuals that are focused on one gene, all cell types, and a chosen metadata
#' attribute of interest.
#'
#' `scgenepanel` exports a multi-panel scRNAseq gene expression plots of
#'   UMAP clusters highlighted with user defined cell types and split view on
#'   user defined sample groups, violin plot with user defined sample groups
#'   along with tabular cell counts, ratios and expression values per sample
#'   group, in one visual
#'
#'
#' @param object A Seurat or SingleCellExperiment object. Seurat objects
#'   will be automatically updated to the current package version if needed
#' @param gene Name of gene to explore gene expression in UMAP, violin plot, and
#'   cell frequency table
#' @param meta_group The metadata column name of the variable to split the UMAP,
#'   violinplot and cell frequency table by. for example, to split by disease
#'   condition
#' @param cell_type_name The cell type identity to highlight in UMAP
#' @param cell_type_colname The metadata column name that contains the cell
#'   identity annotations
#' @param dim_red Dimension reduction embedding to use
#' @param col_palette Color palettes to choose for violinplot panel. Options
#'   are "Tableau" or RColorBrewer qualitative variables like "Dark2", "Paired",
#'   "Set1", "Set2", "Set3", "Accent" etc.
#' @param group_order User defined order of the meta_groups to be displayed
#' @param output_dir Output directory where the image will be saved
#' @return A ggplot object containing a multi-panel visualization with UMAP
#' plots, violin plots, and cell frequency statistics for the specified gene
#' and cell type
#' @importFrom stats quantile
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @examples
#' # Load example pancreatic islet data from package
#' example_file <- system.file(
#'   "extdata",
#'   "human_panc_islets.Rds",
#'   package = "scGenePanel")
#' human_panc_islets <- readRDS(example_file)
#'
#' # Basic usage - create gene panel for insulin expression in beta cells
#' panel <- scgenepanel(
#'   object = human_panc_islets,
#'   gene = "INS",                      # Insulin gene
#'   meta_group = "Age",                # Group by age (14y, 39y, 50y, 59y, 66y)
#'   cell_type_name = "Beta",           # Highlight beta cells
#'   cell_type_colname = "CellTypes"    # Cell type column
#' )
#'
#' \donttest{
#' # Advanced usage with custom colors and ordering
#' panel_advanced <- scgenepanel(
#'   object = human_panc_islets,
#'   gene = "GCG",                      # Glucagon for alpha cells
#'   meta_group = "Sex",                # Group by sex (F/M)
#'   cell_type_name = "Alpha",          # Highlight alpha cells
#'   cell_type_colname = "CellTypes",
#'   col_palette = "Set2",              # Custom color palette
#'   group_order = c("F", "M"),         # Custom group order
#'   output_dir = tempdir()             # Save to temp directory
#' )
#'
#' # Example with different cell type and gene
#' somatostatin_panel <- scgenepanel(
#'   object = human_panc_islets,
#'   gene = "SST",                      # Somatostatin for delta cells
#'   meta_group = "Batch",              # Group by batch
#'   cell_type_name = "Delta",          # Highlight delta cells
#'   cell_type_colname = "CellTypes"
#' )
#' }
#' @export

scgenepanel <- function(object,
                              gene,
                              meta_group,
                              cell_type_name,
                              cell_type_colname,
                              dim_red = "umap",
                              col_palette = "Tableau",
                              group_order = NULL,
                              output_dir = getwd()) {

  # Convert to Seurat object
  seurat_obj <- suppressMessages(.make_seurat(object = object))

  # Check if gene included in object
  .is_gene(seurat_obj, gene = gene)

  # Check if 'cell_type_colname' exists
  .is_cell_type_colname(seurat_obj, cell_type_colname = cell_type_colname)

  # Check if 'cell_type_name' exists
  .is_cell_type_name(
    seurat_obj,
    cell_type_colname = cell_type_colname,
    cell_type_name = cell_type_name
  )

  # Check if 'meta_group' exists
  .is_meta_group(seurat_obj, meta_group = meta_group)

  # Check if dimension reduction embedding exists
  .check_if_dim_red_exists(seurat_obj, dim_red = dim_red)

    # Check cell count adequacy for meaningful analysis
  .check_cell_counts(
    seurat_obj,
    meta_group = meta_group,
    cell_type_colname = cell_type_colname
  )

  # Check for and clean NAs in key metadata columns
  seurat_obj <- .validate_and_clean_nas(
    seurat_obj,
    cell_type_colname = cell_type_colname,
    meta_group = meta_group
  )


  # Check if there are enough cells for meaningful analysis
  meta_data <- seurat_obj@meta.data
  cell_counts <- table(meta_data[[cell_type_colname]], meta_data[[meta_group]])
  if (any(cell_counts < 10)) {
    warning("Some groups have very few cells (< 10). Results may not be meaningful.",
            call. = FALSE)
  }

  # Retrieve group idents to visualize with
  Seurat::Idents(seurat_obj) <- meta_group
  levels_idents <- unique(seurat_obj[[meta_group]][, 1])
  levels_idents <- as.character(levels_idents)


  # UMAP plots
  message("Step 1 - UMAP plots")
  umap <- umap_panel(seurat_obj,
                         cell_type_colname,
                         cell_type_name,
                         meta_group,
                         gene,
                         dim_red,
                         levels_idents,
                         group_order
                         )
  # Violin plots
  message("Step 2 - Violin plots")
  violin <- violin_sig_panel(seurat_obj,
                               cell_type_colname,
                               cell_type_name,
                               meta_group,
                               gene,
                               col_palette,
                               levels_idents)

  # Cell frequency table
  message("Step 3 - Table")
  table <- cellfreq_panel(seurat_obj,
                                cell_type_colname,
                                cell_type_name,
                                meta_group,
                                gene,
                                col_palette)

  message("Step 4 - Assembling panel")

  # Combine figures
  figure  <- suppressWarnings(
    ggpubr::ggarrange(
      umap,
      NULL,
      violin,
      NULL,
      ggpubr::ggarrange(NULL, table, NULL, ncol=3, widths = c(1, 18, 1)),
      heights = c(1.8, 0.1, 3, 0.1, 1.5),
      ncol = 1,
      nrow = 5)
  )

  figure <- ggpubr::annotate_figure(
    figure,
    top = ggpubr::text_grob(
      label = paste0(gene, " expression in ", cell_type_name, " cells"),
      face = "bold",
      size = 20
      ),
    bottom = ggpubr::text_grob(
      label = "Generated by scGenePanel (Shrestha S, et al, 2024).",
      hjust = 1,
      x = 1,
      face = "italic",
      size = 12)
  )

  # Sanitize file name
  cell_type_name_clean <- cell_type_name

  if(grepl("[/?!]", cell_type_name_clean)) {
    cell_type_name_clean <- gsub("[#\\/?!{}<>*'\"|=@\`]", "-", cell_type_name_clean)
    warning("Cell type name contains disallowed characters. Replacing with '-'.")
  }

  if(grepl("[ \t\n]", cell_type_name_clean)) {
    cell_type_name_clean <- gsub("[ ]", "_", cell_type_name_clean)
    warning("Cell type name contains whitespace characters. Replacing with '_'.")
  }

  filename_partial <- paste0(output_dir, "/scGenePanel__", gene, "_", cell_type_name_clean, "_", meta_group, "_", dim_red)

  # Check the operating system and dir slashes
  if (.Platform$OS.type == "unix") {
    # Convert backslashes to forward slashes
    filename_partial <- gsub("\\\\", "/", filename_partial)
  } else if (.Platform$OS.type == "windows") {
    # Convert forward slashes to backslashes
    filename_partial <- gsub("/", "\\\\", filename_partial)
  } else {
    cat("Unsupported operating system, path slashed may be wrong")
  }

  figure <- figure + ggpubr::rremove("grid")

  # Export figure
  filename <- paste0(filename_partial , ".png")
  ggplot2::ggsave(
    filename,
    plot = figure,
    width = 12,
    height = 7,
    scale = 2,
    dpi = 600,
    units = "in",
    limitsize = FALSE,
    bg = "white")
  message(paste0("Step 5 - Plot saved: ", filename))

  filename <- paste0(filename_partial , ".pdf")
  ggplot2::ggsave(
    filename,
    plot = figure,
    width = 12,
    height = 7,
    scale = 2,
    dpi = 600,
    units = "in",
    limitsize = FALSE,
    bg = "white",
    device = cairo_pdf)
  message(paste0("Step 5 - Plot saved: ", filename))

  return(figure)
}
