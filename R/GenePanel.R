






#' scRNAseq multipanel gene expression visual
#'
#' This functions exports a multipanel scRNAseq gene expression plots of UMAP, violin plot with usen defined cell type and condition/groups
#' along with tabular cell counts and ratios per chosen group in one visual
#'
#'
#' @param object A Seurat S4 class object
#' @param gene Name of gene to explore gene expression in UMAP and violin plot by condition of interest
#' @param meta_group The metadata column name that contains the groups to be compared
#' @param cell_type_name The cell type identity that will be highlighted all plots
#' @param col.palette color palettes to choose from. Options are "tableu","varibow" or RColorBrewer qualitative variables like "Dark2","Paired","Set1" etc
#' @param group_order user defined order of the groups to be displayed
#' @param output_dir Output directory where the image will be saved
#' @return Multipanel plots
#' @importFrom magrittr %>%
#' @importFrom Seurat Idemts FeaturePlot
#' @importFrom ggplot2
#' @importFrom dplyr tally
#' @importFrom ggpubr textgrob annotate

create_gene_panel <- function(seurat_object,
                              gene,
                              meta_group,
                              cell_type_name,
                              cell_type_colname,
                              col.palette,
                              group_order = NULL,
                              output_dir="."){

  # Check Seurat object
  Is_Seurat(seurat_object = seurat_object)

  # Check if gene included in object
  Is_gene(seurat_object = seurat_object)

  # Check if 'meta_group' exists
  Is_meta_group(seurat_object = seurat_object)

  # Check if 'cell_type_name' exists
  Is_cell_type_name(seurat_object = seurat_object)

  # Check if 'cell_type_colname' exists
  Is_celltype_colname(seurat_object = seurat_object)

  # retrieve group idents to visualize umaps/violinplots by
  Seurat::Idents(seurat_object) <- meta_group
  levels_idents<-unique(seurat_object[[meta_group]][,1])
  levels_idents<-as.character(levels_idents)

  # loop through group idents for respective umaps
  loop_idents <- function(x) {
    obj_idents <- subset(seurat_object, idents = x)
    p1.d <- suppressMessages(Seurat::FeaturePlot(obj_idents, features = gene,  pt.size=1) +
                               viridis::scale_color_viridis(direction = -1) +
                               ggplot2::labs(title = paste(x,"_subset", sep="")))
    p1.d[[1]]$layers[[1]]$aes_params$alpha = .5
    p1.d.umap <- data.frame(obj_idents@reductions$umap@cell.embeddings)
    Seurat::Idents(obj_idents) <- cell_type_colname
    plot1 <- p1.d + ggforce::geom_mark_ellipse(ggplot2::aes(x = p1.d.umap$UMAP_1, y = p1.d.umap$UMAP_2,
                                                   fill = Seurat::Idents(obj_idents),
                                                   filter = Seurat::Idents(obj_idents) == cell_type_name),
                                               alpha = 0.1)
    suppressMessages(plot1)
  }

  panels <- lapply(levels_idents, loop_idents)
  panel_figure <- cowplot::plot_grid(plotlist = panels, ncol = length(levels_idents))

  # Violin plot by meta_group

  select_col <-  discrete_col_palette(num_colors = length(levels_idents), palette = col.palette)
  plot2 <- Seurat::VlnPlot(object = seurat_object, features = gene,group.by = cell_type_colname, split.by = meta_group, cols = select_col, pt.size=-1)


  # table of cell counts/expression ratios
  Seurat::Idents(seurat_object) <- cell_type_colname
  selected_cluster_cells <- subset(seurat_object, idents=cell_type_name)
  #cc_tally <- selected_cluster_cells@meta.data %>% group_by(!!! syms(meta_group)) %>% tally()
  cell_counts_tally <- selected_cluster_cells@meta.data %>% tidylog::group_by_at(meta_group) %>% dplyr::tally()

  # if gene has no expression over 0.25 (no cells), it will error out. So tryCatch it
  selected_cells_expressed <- tryCatch({

    #p <- GetAssayData(object = selected_cluster_cells, assay = "RNA", slot = "data")[my.gene,]
    selected_cells_expressed <- selected_cluster_cells[, Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ] > 0.25]

  }, error = function(e)
  {
    selected_cells_expressed <- NULL
  }
  #cat('Yet another error replaced by NA \\n')

  )

  if(is.null(selected_cells_expressed)){
    message(paste0("No expression detected for: ", gene))
    #not_expressed_features <- append(not_expressed_features, gene)

  }else{

    selected_cells_expressed_tally <- selected_cells_expressed@meta.data %>% tidylog::group_by_at(meta_group) %>% dplyr::tally()     # tally cells expressing my.gene
    # merge DFs
    cc_table <- merge(cell_counts_tally, selected_cells_expressed_tally,
                      by = meta_group,
                      suffixes = c("_cells", "_expressing"),
                      all = TRUE)

    cc_table$ratio = cc_table$n_expressing / cc_table$n_cells # ratio
    cc_table[is.na(cc_table)] <- 0
    t1 <- ggpubr::ggtexttable(cc_table, rows = NULL)

    #filename <- paste0(output_dir, "genePanel__", my.gene, '_', cell_type_clean, '.csv')
    #write.csv(cc_table, file=filename)
    #t1

    figure  <- suppressWarnings(ggpubr::ggarrange(panel_figure, plot2, t1, heights = c(1.8, 3, 1), ncol = 1, nrow = 3))
    figure <- ggpubr::annotate_figure(figure,
                                      top = ggpubr::text_grob(paste0("Gene: ", gene), face = "bold", size = 18),
                                      bottom = ggpubr::text_grob("Generated by genemap.R",
                                                         hjust = 1, x = 1, face = "italic", size = 10)
    )
    figure

    filename <- paste0(output_dir, "genePanel__", gene, '_', cell_type_name, '.png')
    ggplot2::ggsave(filename, width = 12, height = 5, scale = 2, dpi=300, units="in")


  }

  #return(not_expressed_features)


}

#test with any seurat object
# create_gene_panel(seurat_object = data,
#                   gene = "INS",
#                   cell_type_name = "Beta",
#                   meta_group = "Age",
#                   cell_type_colname = "CellTypes",
#                   col.palette = "Dark2",
#                   output_dir = "/Users/shristishrestha/Documents/scViz/results/")



