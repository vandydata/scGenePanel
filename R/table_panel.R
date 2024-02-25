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

table_panel <- function(object,
                         cell_type_colname,
                         cell_type_name,
                         meta_group,
                         gene,
                         col_palette,
                         levels_idents
                         ) {


  # Table of cell counts/expression ratios
  Seurat::Idents(seurat_obj) <- cell_type_colname
  selected_cluster_cells <- subset(seurat_obj, idents = cell_type_name)
  #cc_tally <- selected_cluster_cells@meta.data %>% group_by(!!! syms(meta_group)) %>% tally()
  cell_counts_tally <- suppressMessages(selected_cluster_cells@meta.data %>% tidylog::group_by_at(meta_group) %>% dplyr::tally())

  # if gene has no expression over 0.25 ,print message of "No detection"
  selected_cells_expressed <- tryCatch({

    selected_cells_expressed <- selected_cluster_cells[, Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ] > 0.25]

  }, error = function(e) {
    selected_cells_expressed <- NULL
  }

  )

  # TODO shouldn't this be tested earlier?
  if (is.null(selected_cells_expressed)) {
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

    cc_table$ratio <- cc_table$n_expressing / cc_table$n_cells * 100 # ratio
    cc_table[is.na(cc_table)] <- 0
    names(cc_table)[names(cc_table) == 'ratio'] <- 'pct.expressed'

    # Add quantiles for gene expression per meta_group
    celltype <- NULL # need this to remove “no visible binding” note, peculiarity of using dplyr
    obj_meta <- seurat_obj@meta.data
    meta_subset <- obj_meta[, c(cell_type_colname, meta_group)]
    names(meta_subset) <- c("celltype", "group")
    meta_subset <- meta_subset %>% filter(celltype == cell_type_name)

    exp_orig <- Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ]
    exp_orig <- as.data.frame(exp_orig)

    reorder_exp <- exp_orig[match(rownames(meta_subset), rownames(exp_orig)), ]

    meta_subset_final <- meta_subset %>% dplyr::mutate(gene_exp = reorder_exp)

    result <- do.call("cbind", tapply(meta_subset_final$gene_exp, meta_subset_final$group, quantile, probs=seq(0, 1, 0.25)))
    result <- t(result)[,2:5]


    #combine the table and qunatile results
    combined_cc_table <- cbind(cc_table, result)
    names(combined_cc_table) <- c(meta_group, "n_cells", "n_expressing(>0.25)", "%expressed", "25%(quantile)", "50%(quantile)", "75%(quantile)", "100%(quantile)")

    #increase the size of table
    table_theme <- ggpubr::ttheme(
      base_style = "default",
      base_size = 20,
      base_colour = "black",
      padding = ggplot2::unit(c(6, 6), "mm"),
      colnames.style = ggpubr::colnames_style(size = 20),
      rownames.style = ggpubr::rownames_style(size = 20),
      tbody.style = ggpubr::tbody_style(size = 20)
    )

    # create table
    panel_table <- ggpubr::ggtexttable(combined_cc_table, rows = NULL, theme = table_theme)

    return(panel_table)
  }
}
