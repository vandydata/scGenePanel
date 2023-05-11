

#' scRNAseq multipanel gene expression visual by sample groups
#'
#' `create_gene_panel` exports a multipanel scRNAseq gene expression plots of UMAP clusters highlighted with user defined cell types and split view on user defined sample groups, violin plot with user defined sample groups
#' along with tabular cell counts, ratios and expression values per sample group, in one visual
#'
#'
#' @param object A Seurat or SingleCellExperiment
#' @param gene Name of gene to explore gene expression in UMAP, violin plot, and cell frequency table
#' @param meta_group The metadata column name of the variable to split the UMAP, violinplot and cell frequency table by. for example, to split by disease condition
#' @param cell_type_name The cell type identity to highlight in UMAP
#' @param cell_type_colname The metadata column name that contains the cell identity annotations
#' @param col.palette Color palettes to choose for violinplot panel. Options are "tableu" or RColorBrewer qualitative variables like "Dark2","Paired","Set1","Set2", "Set3", "Accent" etc.
#' @param group_order User defined order of the meta_groups to be displayed
#' @param output_dir Output directory where the image will be saved
#' @return Multipanel plots
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
                              col.palette,
                              group_order = NULL,
                              output_dir="."){

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

  # Retrieve group idents to visualize umaps/violinplots by
  Seurat::Idents(seurat_obj) <- meta_group
  levels_idents <- unique(seurat_obj[[meta_group]][,1])
  levels_idents <- as.character(levels_idents)

  # Loop through group idents for respective umaps
  loop_idents <- function(x) {
    obj_idents <- subset(seurat_obj, idents = x)
    p1.d <- suppressMessages(Seurat::FeaturePlot(obj_idents, features = gene,  pt.size=1) +
                               viridis::scale_color_viridis(direction = -1) +
                               ggplot2::theme(plot.title = ggplot2::element_text(size=20),
                                              axis.title.x = ggplot2::element_text(size=20, face="bold"),
                                              axis.title.y = ggplot2::element_text(size=20, face="bold"),
                                              axis.text = ggplot2::element_text(size = 20,face = 'bold'),
                                              axis.text.x = ggplot2::element_text(size = 20,face = 'bold'),
                                              axis.text.y = ggplot2::element_text(size = 20,face = 'bold'),
                                              legend.text = ggplot2::element_text(size = 16,face = 'bold'),
                                              legend.title = ggplot2::element_text(size = 16,face = 'bold')
                               )+
                               ggplot2::labs(title = paste(x,"_subset", sep="")))
    p1.d[[1]]$layers[[1]]$aes_params$alpha = .5
    p1.d.umap <- data.frame(obj_idents@reductions$umap@cell.embeddings)
    Seurat::Idents(obj_idents) <- cell_type_colname
    plot1 <- p1.d + ggforce::geom_mark_ellipse(ggplot2::aes(x = p1.d.umap$UMAP_1, y = p1.d.umap$UMAP_2,
                                                   fill = Seurat::Idents(obj_idents),
                                                   filter = Seurat::Idents(obj_idents) == cell_type_name),
                                               alpha = 0.1) + ggplot2::scale_fill_discrete(name = "Cell Type")
    suppressMessages(plot1)
  }


  if (is.null(group_order)){
    panels <- lapply(levels_idents, loop_idents)
    panel_figure <- cowplot::plot_grid(plotlist = panels, ncol = length(levels_idents))
  } else {
    levels_idents <- levels_idents[order(match(levels_idents, group_order))]
    panels <- lapply(levels_idents, loop_idents)
    panel_figure <- cowplot::plot_grid(plotlist = panels, ncol = length(levels_idents))
    }


  # Violin plot by meta_group

  select_col <-  discrete_col_palette(num_colors = length(levels_idents), palette = col.palette)
  plot2 <- suppressWarnings(Seurat::VlnPlot(object = seurat_obj, features = gene, group.by = cell_type_colname, split.by = meta_group, cols = select_col, pt.size=-1)+
    # ggplot2::stat_summary(fun = median, fun.min = median, fun.max = median,
    #              geom = "crossbar",
    #              width = 0.2,
    #              position = ggplot2::position_dodge(width = 1))+
    ggpubr::stat_compare_means(method = "anova", label.y = 5, label = "p.signif")+
    #ggplot2::geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent")+
    ggplot2::geom_violin(trim = F, alph=0.5, scale = "width",draw_quantiles = c(0.25, 0.5, 0.75))+
    #ggplot2::geom_boxplot(width=0.1, color="black", alpha=0.2, position=ggplot2::position_dodge(1))+
    ggplot2::theme(plot.title = ggplot2::element_text(size=20),
                   axis.title.x = ggplot2::element_text(size=20, face="bold"),
                   axis.title.y = ggplot2::element_text(size=20, face="bold"),
                   axis.text = ggplot2::element_text(size = 20,face = 'bold'),
                   axis.text.x = ggplot2::element_text(size = 20,face = 'bold'),
                   axis.text.y = ggplot2::element_text(size = 20,face = 'bold'),
                   legend.text = ggplot2::element_text(size = 16,face = 'bold'),
                   legend.title = ggplot2::element_text(size = 16,face = 'bold')
                   )+
    ggplot2::xlab(cell_type_colname))



  # Table of cell counts/expression ratios
  Seurat::Idents(seurat_obj) <- cell_type_colname
  selected_cluster_cells <- subset(seurat_obj, idents=cell_type_name)
  #cc_tally <- selected_cluster_cells@meta.data %>% group_by(!!! syms(meta_group)) %>% tally()
  cell_counts_tally <- suppressMessages(selected_cluster_cells@meta.data %>% tidylog::group_by_at(meta_group) %>% dplyr::tally())

  # if gene has no expression over 0.25 ,print message of "No detection"
  selected_cells_expressed <- tryCatch({

    selected_cells_expressed <- selected_cluster_cells[, Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ] > 0.25]

  }, error = function(e)
  {
    selected_cells_expressed <- NULL
  }

  )

  if (is.null(selected_cells_expressed)){
    message(paste0("No expression detected for: ", gene))

  } else {

    selected_cells_expressed_tally <- suppressMessages(selected_cells_expressed@meta.data %>%
      tidylog::group_by_at(meta_group) %>%
      dplyr::tally())     # tally cells expressing the gene

   # Merge cell frequency and expression tally
    cc_table <- merge(cell_counts_tally, selected_cells_expressed_tally,
                      by = meta_group,
                      suffixes = c("_cells", "_expressing"),
                      all = TRUE)

    cc_table$ratio = cc_table$n_expressing / cc_table$n_cells * 100 # ratio
    cc_table[is.na(cc_table)] <- 0
    names(cc_table)[names(cc_table) == 'ratio'] <- 'pct.expressed'

   # Add quantiles for gene expression per meta_group
    celltype <- NULL # need this to remove “no visible binding” note, peculiarity of using dplyr
    obj_meta <- seurat_obj@meta.data
    meta_subset <- obj_meta[,c( cell_type_colname, meta_group)]
    names(meta_subset) <- c("celltype","group")
    meta_subset <- meta_subset %>% filter(celltype == cell_type_name)

    exp_orig <- Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ]
    exp_orig <- as.data.frame(exp_orig)

    reorder_exp <- exp_orig[match(rownames(meta_subset), rownames(exp_orig)),]

    meta_subset_final <- meta_subset %>% dplyr::mutate(gene_exp = reorder_exp)

    result <- do.call("cbind", tapply(meta_subset_final$gene_exp, meta_subset_final$group,quantile,probs=seq(0,1,0.25)))
    result <- t(result)[,2:5]


    combined_cc_table <- cbind(cc_table, result)
    names(combined_cc_table) <- c(meta_group,"n_cells", "n_expressing(>0.25)", "%expressed", "25%(quantile)", "50%(quantile)", "75%(quantile)", "100%(quantile)" )

    table_theme <- ggpubr::ttheme(
      base_style = "default",
      base_size = 20,
      base_colour = "black",
      padding = ggplot2::unit(c(6, 6), "mm"),
      colnames.style = ggpubr::colnames_style(size = 20),
      rownames.style = ggpubr::rownames_style(size = 20),
      tbody.style = ggpubr::tbody_style(size = 20)
    )

    t1 <- ggpubr::ggtexttable(combined_cc_table, rows = NULL, theme = table_theme)


    # Combine figures
    figure  <- suppressWarnings(ggpubr::ggarrange(panel_figure, plot2, t1, heights = c(1.8, 3, 1), ncol = 1, nrow = 3))
    figure <- ggpubr::annotate_figure(figure,
                                      top = ggpubr::text_grob(paste0("Gene: ", gene), face = "bold", size = 18),
                                      bottom = ggpubr::text_grob("Generated by genemap.R",
                                                         hjust = 1, x = 1, face = "italic", size = 10)
    )
    figure

    # Export figure
    filename <- paste0(output_dir, "genePanel__", gene, '_', cell_type_name, '.png')
    ggplot2::ggsave(filename, width = 12, height = 7, scale = 2, dpi=1200, units="in")


  }

}




