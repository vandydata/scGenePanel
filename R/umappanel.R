#' scRNAseq multi-panel gene expression visual
#'
#' This functions exports a multi-panel scRNAseq gene expression plots of UMAP,
#'   violin plot with usen defined cell type and condition/groups along with
#'   tabular cell counts and ratios per chosen group in one visual
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
#' @param col_palette Color palettes to choose for violinplot panel. Options are
#'   "tableu" or RColorBrewer qualitative variables like "Dark2", "Paired",
#'   "Set1", "Set2", "Set3", "Accent" etc.
#' @return multi-panel plots
#' @importFrom stats quantile
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @export

umap_panel <- function(object,
                       cell_type_colname,
                       cell_type_name,
                       meta_group,
                       gene,
                       col_palette) {

  # Convert to Seurat object
  seurat_obj <- suppressMessages(make_seurat(object = object))

  # Check if 'cell_type_colname' exists
  Is_celltype_colname(seurat_obj, cell_type_colname = cell_type_colname)

  # Check if 'cell_type_name' exists
  Is_cell_type_name(seurat_obj, cell_type_colname = cell_type_colname, cell_type_name = cell_type_name)

  # Check if 'meta_group' exists
  Is_meta_group(seurat_obj, meta_group = meta_group)

  # Check if gene included in object
  Is_gene(seurat_obj, gene = gene)

  # retrieve group idents to visualize umaps/violinplots by
  Seurat::Idents(seurat_obj) <- meta_group
  levels_idents <- unique(seurat_obj[[meta_group]][, 1])
  levels_idents <- as.character(levels_idents)

  # loop through group idents for respective umaps
  loop_idents <- function(x) {
    obj_idents <- subset(seurat_obj, idents = x)
    p1.d <- suppressMessages(Seurat::FeaturePlot(obj_idents, features = gene,  pt.size = 1) +
                               viridis::scale_color_viridis(direction = -1) +
                               ggplot2::labs(title = paste(x, "_subset", sep = "")))
    p1.d[[1]]$layers[[1]]$aes_params$alpha <- .5
    p1.d.umap <- data.frame(obj_idents@reductions$umap@cell.embeddings)
    Seurat::Idents(obj_idents) <- cell_type_colname
    plot1 <- p1.d + ggforce::geom_mark_ellipse(ggplot2::aes(x = p1.d.umap$UMAP_1, y = p1.d.umap$UMAP_2,
                                                            fill = Seurat::Idents(obj_idents),
                                                            filter = Seurat::Idents(obj_idents) == cell_type_name),
                                               alpha = 0.1)
    suppressMessages(plot1)
  }

  panels <- lapply(levels_idents, loop_idents)
  panel_figure <- cowplot::plot_grid(plotlist = panels, ncol = 1)
  panel_figure

}
