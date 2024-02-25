#' scRNAseq multi-panel gene expression visual
#'
#' This functions exports a multi-panel scRNAseq gene expression plots of UMAP,
#'   violin plot with usen defined cell type and condition/groups along with
#'   tabular cell counts and ratios per chosen group in one visual
#'
#'
#' @param object A Seurat or SingleCellExperiment
#' @param gene Name of gene to explore gene expression in UMAP, violin plot,
#'   and cell frequency table
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

violin_panel <- function(object,
                         cell_type_colname,
                         cell_type_name,
                         meta_group,
                         gene,
                         col_palette,
                         levels_idents
                         ) {

  select_col <-  discrete_col_palette(num_colors = length(levels_idents), palette = col_palette)
  #max_value <- max(Seurat::GetAssayData(seurat_obj[["RNA"]])[gene, ], slot = "data")
  max_value <- Seurat::VlnPlot(object = seurat_obj, features = gene, group.by = cell_type_colname)
  max_value <- max_value$data
  max_value <- max(max_value[,1])
  plot2 <- suppressWarnings(Seurat::VlnPlot(object = seurat_obj, features = gene, group.by = cell_type_colname, split.by = meta_group, cols = select_col, pt.size = -1) +
                              #this is one-way ANOVA test, a significant p-value indicates that some of the group means are different, but we don’t know which pairs of groups are different.
                              #ggpubr::stat_compare_means(method = "anova", label.y = max_value, label = "p.signif", size = 8) + # Add global p-value #currently this does not work on all plot
                              # @JP yes, b/c it's not group-aware
                              ggplot2::geom_violin(trim = FALSE, alph = 0.5, scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
                              ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                                             axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
                                             axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
                                             axis.text = ggplot2::element_text(size = 20, face = "bold"),
                                             axis.text.x = ggplot2::element_text(size = 20, face = "bold"),
                                             axis.text.y = ggplot2::element_text(size = 20, face = "bold"),
                                             legend.text = ggplot2::element_text(size = 16, face = "bold"),
                                             legend.title = ggplot2::element_text(size = 16, face = "bold")
                              ) +
                              ggplot2::xlab(cell_type_colname))

  return(plot2)
}
