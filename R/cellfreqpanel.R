

#' scRNAseq multipanel gene expression visual
#'
#' This functions exports a multipanel scRNAseq gene expression plots of UMAP, violin plot with usen defined cell type and condition/groups
#' along with tabular cell counts and ratios per chosen group in one visual
#'
#'
#' @param object A Seurat or SingleCellExperiment
#' @param gene Name of gene to explore gene expression in UMAP, violin plot, and cell frequency table
#' @param meta_group The metadata column name of the variable to split the UMAP, violinplot and cell frequency table by. for example, to split by disease condition
#' @param cell_type_name The cell type identity to highlight
#' @param cell_type_colname The metadata column name that contains the cell identity annotations
#' @param col.palette Color palettes to choose for violinplot panel. Options are "tableu","varibow" or RColorBrewer qualitative variables like "Dark2","Paired","Set1" etc
#' @return Multipanel plots
#' @importFrom stats quantile
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @export

cellfreq_panel <- function(object,
                              cell_type_colname,
                              cell_type_name,
                              meta_group,
                              gene,
                              col.palette
                           ){

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

    selected_cells_expressed_tally <- selected_cells_expressed@meta.data %>%
      tidylog::group_by_at(meta_group) %>%
      dplyr::tally()     # tally cells expressing the gene

    # Merge cell frequency and expression tally
    cc_table <- merge(cell_counts_tally, selected_cells_expressed_tally,
                      by = meta_group,
                      suffixes = c("_cells", "_expressing"),
                      all = TRUE)

    cc_table$ratio = cc_table$n_expressing / cc_table$n_cells # ratio
    cc_table[is.na(cc_table)] <- 0
    #t1 <- ggpubr::ggtexttable(cc_table, rows = NULL)

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
    names(combined_cc_table) <- c(meta_group,"n_cells", "n_expressing(>0.25)", "ratio", "25%(quantile)", "50%(quantile)", "75%(quantile)", "100%(quantile)" )

    t1 <- ggpubr::ggtexttable(combined_cc_table, rows = NULL)
    print(t1)





  }

  #return(not_expressed_features)


}




