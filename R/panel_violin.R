#'
#' This functions generates violin plots for scGenePanel
#'
#'
#' @param object A Seurat object
#' @param cell_type_colname The metadata column name that contains the cell
#'   identity annotations
#' @param cell_type_name The cell type identity to highlight
#' @param meta_group The metadata column name of the variable to split with
#' @param gene Name of gene to explore gene expression
#' @param col_palette Color palettes to choose for violinplot panel. Options are
#'   "tableu" or RColorBrewer qualitative variables like "Dark2", "Paired",
#'   "Set1", "Set2", "Set3", "Accent" etc.
#'
#' @return ggplot2 violing plot object
#' @importFrom stats quantile
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @keywords internal
#' @noRd

violin_panel <- function(seurat_obj,
                         cell_type_colname,
                         cell_type_name,
                         meta_group,
                         gene,
                         col_palette,
                         levels_idents
) {

  select_col <- discrete_col_palette(num_colors = length(levels_idents), palette = col_palette)

  # Get the data from Seurat object
  plot_data <- Seurat::FetchData(seurat_obj, vars = c(gene, cell_type_colname, meta_group))

  # Rename columns for easier reference
  colnames(plot_data) <- c("gene_expression", "cell_type", "meta_group")

  # Create the plot from scratch
  #this is one-way ANOVA test, a significant p-value indicates that some of the group means are different, but we don’t know which pairs of groups are different.
  #ggpubr::stat_compare_means(method = "anova", label.y = max_value, label = "p.signif", size = 8) + # Add global p-value #currently this does not work on all plot
  # @JP yes, b/c it's not group-aware
                              
  plot2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cell_type,
                                                   y = gene_expression,
                                                   fill = meta_group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.5, scale = "width",
                         draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_fill_manual(values = select_col) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 20, face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 20, face = "bold"),
                   legend.text = ggplot2::element_text(size = 16, face = "bold"),
                   legend.title = ggplot2::element_text(size = 16, face = "bold")) +
    ggplot2::xlab(cell_type_colname) +
    ggplot2::ylab(paste("Expression of", gene)) +
    ggplot2::labs(fill = meta_group)

  return(plot2)
}


violin_sig_panel <- function(seurat_obj,
                          cell_type_colname,
                          cell_type_name,
                          meta_group,
                          gene,
                          col_palette,
                          levels_idents
) {


  # To perform stat tests within groups, we have to actually create the groups
  # and then perform the test within the group, as well as determine max value for
  # each
  loop_idents <- function(x) {

    #message(x)
    select_col <-  discrete_col_palette(num_colors = length(levels_idents), palette = col_palette)

    Seurat::Idents(seurat_obj) <- cell_type_colname
    obj_idents <- subset(seurat_obj, idents = x)

    #p <- suppressWarnings(Seurat::VlnPlot(object = obj_idents, features = gene, split.by = meta_group, cols = select_col, pt.size = -1))

    p <- scCustomize::VlnPlot_scCustom(seurat_object = obj_idents, features = gene, group.by = meta_group, pt.size = 0)

    #p$data

    # Get min/max values for plotting
    max_value <- max(p$data[,1])
    min_value <- min(p$data[,1])
    ylim_max <- max_value + (max_value * 0.5)
    ylim_min <- min_value - (min_value * 0.1)

    # Get unique values for current object subset
    levels_idents_subset <- as.character(unique(obj_idents@meta.data[[meta_group]]))

    # only perform pairwise tests if pairwise data exists
    if(length(levels_idents_subset) > 1){
      # Assemble all pairwise comparisons
      pairwise_combinations <- combn(sort(levels_idents_subset), 2, simplify = FALSE)
      #pairwise_combinations <- list( c("50y", "66y"), c("59y", "66y") )
      # pairwise_combinations HAVE to match x-axis values
      #
      if(TRUE){
        # Display ns, *, **, *** for significance
        # {p.format}{p.signif} also workds
        p <- p + ggpubr::geom_pwc(ggplot2::aes(group = ident), method = "t_test", label = "{p.signif}",   hide.ns = TRUE) + ggplot2::ylim(ylim_min, ylim_max)

      }else{
        # Display p-values, quite busy so disabled for now
        p <- p + ggpubr::geom_pwc(ggplot2::aes(group = ident), method = "t_test", label = "{p.format}",   hide.ns = TRUE) + ggplot2::ylim(ylim_min, ylim_max)
      }
    }


    p[[1]]$labels$x <- ''
    p[[1]]$labels$y <- ''
    p[[1]]$labels$title <- x


    suppressMessages(p)
  }

  cell_types <- unique(seurat_obj[[cell_type_colname]])
  cell_types <- cell_types$CellTypes
  panels <- lapply(cell_types, loop_idents)
  #panels[[2]]

  # Remove legends
  for(i in 1:length(panels)){
    panels[[i]] <- panels[[i]] + ggplot2::theme(legend.position = "none")
  }

  # transitioned to patchwork since cowplot wasn't aligning x-axis of each plot
  # correctly, due to different heighs of the rotated labels.
  panel_figure <- patchwork::wrap_plots(panels, ncol = length(panels)) + patchwork::plot_layout(guides = 'collect', ncol = length(panels))
  panel_figure

  if(FALSE){
    # Assemble multiplot panel
    panel_figure <- cowplot::plot_grid(plotlist = panels, ncol = length(cell_types))

    #create common x and y labels
    y.grob <- grid::textGrob("UMAP_2",gp=grid::gpar(fontface="bold", col="black", fontsize=15), rot=90)
    x.grob <- grid::textGrob("UMAP_1",gp=grid::gpar(fontface="bold", col="black", fontsize=15))
    panel_figure <- gridExtra::grid.arrange(gridExtra::arrangeGrob(panel_figure, left = y.grob, bottom = x.grob))
    panel_figure
  }

  return(panel_figure)
}
