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

umap_panel <- function(seurat_obj,
                       cell_type_colname,
                       cell_type_name,
                       meta_group,
                       gene,
                       levels_idents,
                       group_order=NULL) {

  # loop through group idents for respective umaps
  loop_idents <- function(x) {
    #Idents(seurat_obj) <- meta_group
    obj_idents <- subset(seurat_obj, idents = x)
    p1.d <- suppressMessages(scCustomize::FeaturePlot_scCustom(obj_idents, features = gene,  pt.size = 1) +
                               ggplot2::labs(title = x))
    p1.d[[1]]$layers[[1]]$aes_params$alpha <- .5

    p1.d.umap <- data.frame(obj_idents@reductions$umap@cell.embeddings)
    Seurat::Idents(obj_idents) <- cell_type_colname
    plot1 <- p1.d + ggforce::geom_mark_ellipse(
      ggplot2::aes(x = p1.d.umap$UMAP_1, y = p1.d.umap$UMAP_2,
                   fill = Seurat::Idents(obj_idents),
                   filter = Seurat::Idents(obj_idents) == cell_type_name),
      alpha = 0.1)


    plot1[[1]]$labels$fill <- "Cell type"
    plot1[[1]]$labels$colour <- "Expression"
    plot1[[1]]$labels$x <- ''
    plot1[[1]]$labels$y <- ''
    plot1[[1]] + guides(fill = guide_legend(order = 1),
                        colour = guide_legend(order = 2))

    plot1 <- plot1 + ggplot2::theme(aspect.ratio=1)
    suppressMessages(plot1)
  }

  # V1
  #panels <- lapply(levels_idents, loop_idents)
  #panel_figure <- cowplot::plot_grid(plotlist = panels, nrow = 1)
  #return(panel_figure)




  # Generate panels
  if(!is.null(group_order)){
    levels_idents <- levels_idents[order(match(levels_idents, group_order))]
  }
  panels <- lapply(levels_idents, loop_idents)

  # Remove legend from all panels except the last one
  for(i in 1:length(panels)){
    if(i != length(levels_idents)){
      panels[[i]] <- panels[[i]] + theme(legend.position = "none")
    }
  }

  # Assemble multiplot panel
  panel_figure <- cowplot::plot_grid(plotlist = panels, ncol = length(levels_idents))

  #create common x and y labels
  y.grob <- grid::textGrob("UMAP_2",gp=grid::gpar(fontface="bold", col="black", fontsize=15), rot=90)
  x.grob <- grid::textGrob("UMAP_1",gp=grid::gpar(fontface="bold", col="black", fontsize=15))
  panel_figure <- gridExtra::grid.arrange(gridExtra::arrangeGrob(panel_figure, left = y.grob, bottom = x.grob))


  return(panel_figure)

}

