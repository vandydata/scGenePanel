

#' scGenePanel - A method for generating scRNAseq multi-panel gene expression
#' visuals that are focused on one gene, all cell types, and a chosen metadata
#' attribute of interest.
#'
#' `create_gene_panel` exports a multi-panel scRNAseq gene expression plots of
#'   UMAP clusters highlighted with user defined cell types and split view on
#'   user defined sample groups, violin plot with user defined sample groups
#'   along with tabular cell counts, ratios and expression values per sample
#'   group, in one visual
#'
#'
#' @param object A Seurat or SingleCellExperiment
#' @param gene Name of gene to explore gene expression in UMAP, violin plot, and
#'   cell frequency table
#' @param meta_group The metadata column name of the variable to split the UMAP,
#'   violinplot and cell frequency table by. for example, to split by disease
#'   condition
#' @param cell_type_name The cell type identity to highlight in UMAP
#' @param cell_type_colname The metadata column name that contains the cell
#'   identity annotations
#' @param col_palette Color palettes to choose for violinplot panel. Options
#'   are "Tableau" or RColorBrewer qualitative variables like "Dark2", "Paired",
#'   "Set1", "Set2", "Set3", "Accent" etc.
#' @param group_order User defined order of the meta_groups to be displayed
#' @param output_dir Output directory where the image will be saved
#' @return Multi-panel plots
#' @importFrom stats quantile
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @export


create_gene_panel <- function(object,
                              gene,
                              meta_group,
                              cell_type_name,
                              cell_type_colname,
                              col_palette = "Tableau",
                              group_order = NULL,
                              output_dir = getwd()) {


  # Convert to Seurat object
  seurat_obj <- suppressMessages(.make_seurat(object = object))

  # Check if 'cell_type_colname' exists
  .is_celltype_colname(seurat_obj, cell_type_colname = cell_type_colname)

  # Check if 'cell_type_name' exists
  .is_cell_type_name(seurat_obj, cell_type_colname = cell_type_colname, cell_type_name = cell_type_name)

  # Check if 'meta_group' exists
  .is_meta_group(seurat_obj, meta_group = meta_group)

  # Check if gene included in object
  .is_gene(seurat_obj, gene = gene)

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
  table <- cellfreq_panel(object,
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

  figure <- ggpubr::annotate_figure(figure,
                                    top = ggpubr::text_grob(
                                      label = paste0(gene, " in ", cell_type_name, " cells"), face = "bold", size = 20
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
    warning(paste0("Cell type name contains disallowed characters. Replacing with '-'."))
  }

  if(grepl("[ \t\n]", cell_type_name_clean)) {
    cell_type_name_clean <- gsub("[ ]", "_", cell_type_name_clean)
    warning(paste0("Cell type name contains whitespace characters. Replacing with '_'."))
  }

  filename <- paste0(output_dir, "/scGenePanel__", gene, "_", cell_type_name_clean, "_", meta_group, ".png")

  # Check the operating system and dir slashes
  if (.Platform$OS.type == "unix") {
    # Convert backslashes to forward slashes
    filename <- gsub("\\\\", "/", filename)
  } else if (.Platform$OS.type == "windows") {
    # Convert forward slashes to backslashes
    filename <- gsub("/", "\\\\", filename)
  } else {
    # Handle other operating systems if needed
    cat("Unsupported operating system, path slashed may be wrong")
  }

  figure <- figure + ggpubr::rremove("grid")

  # Export figure
  ggplot2::ggsave(filename, plot = figure, width = 12, height = 7, scale = 2, dpi = 600, units = "in", limitsize = FALSE, bg = "white")
  message(paste0("Step 5 - Plot saved: ", filename))



  return(figure)

}
