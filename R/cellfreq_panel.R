#'
#' This functions generates cell frequency table for scGenePanel
#'
#'
#' @param seurat_obj A Seurat object
#' @param cell_type_colname The metadata column name that contains the cell
#' identity annotations
#' @param cell_type_name The cell type identity to highlight
#' @param meta_group The metadata column name of the variable to split the UMAP,
#' violinplot and cell frequency table by. for example, to split by disease
#' condition
#' @param gene Name of gene to explore gene expression in UMAP, violin plot and
#' cell frequency table
#' @param col_palette Color palettes to choose for violinplot panel. Options
#' are "tableu","varibow" or RColorBrewer qualitative variables like "Dark2",
#' "Paired", "Set1" etc
#' @param levels_idents The levels of the metadata column to split with
#'
#' @return ggpubr table
#' @importFrom stats quantile
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @export

cellfreq_panel <- function(seurat_obj,
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
  cell_counts_tally <- suppressMessages(
    selected_cluster_cells@meta.data %>%
      dplyr::group_by_at(meta_group) %>%
      dplyr::tally()
  )

  # if gene has no expression over 0.25, print message of "No detection"
  selected_cells_expressed <- tryCatch({
    check <- Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ] > 0.25
    selected_cells_expressed <- selected_cluster_cells[, check]

  }, error = function(e) {
    selected_cells_expressed <- NULL
  }

  )

  if (is.null(selected_cells_expressed)) {
    message(paste0("No expression detected for: ", gene))
  } else {

    selected_cells_expressed_tally <- selected_cells_expressed@meta.data %>%
      dplyr::group_by_at(meta_group) %>%
      dplyr::tally()     # tally cells expressing the gene

    # Merge cell frequency and expression tally
    cc_table <- merge(cell_counts_tally, selected_cells_expressed_tally,
                      by = meta_group,
                      suffixes = c("_cells", "_expressing"),
                      all = TRUE)

    cc_table$ratio <- cc_table$n_expressing / cc_table$n_cells # ratio
    cc_table[is.na(cc_table)] <- 0
    names(cc_table)[names(cc_table) == 'ratio'] <- 'pct.expressed'
    #t1 <- ggpubr::ggtexttable(cc_table, rows = NULL)

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

    result <- do.call("cbind", tapply(
      meta_subset_final$gene_exp,
      meta_subset_final$group,
      quantile,
      probs = seq(0, 1, 0.25)
    ))
    result <- t(result)[, 2:5]

    combined_cc_table <- cbind(cc_table, result)
    names(combined_cc_table) <- c(
      meta_group,
      "n_cells",
      "n_expressing(>0.25)",
      "% expressed",
      "25% (quantile)",
      "50% (quantile)",
      "75% (quantile)",
      "100% (quantile)"
    )

    # Add titles andd footnote
    # :::::::::::::::::::::::::::::::::::::::::::::::::::
    # Wrap subtitle into multiple lines using strwrap()
    main_title <- paste0("Metrics of ", gene, " expression in ", cell_type_name, " cells per ", meta_group)
    subtitle <- paste(
      "n_cells = cell frequency per group",
      "n_expressing = cells expression above 0.25",
      "%expressed = n_expressing/n_cells*100", sep = "\n"
    ) %>%
      strwrap(width = 80) %>%
      paste(collapse = "\n")

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

    t1 <- ggpubr::ggtexttable(combined_cc_table, rows = NULL, theme = table_theme)
    t1 %>%
      ggpubr::tab_add_title(text = subtitle, face = "plain", size = 10) %>%
      ggpubr::tab_add_title(text = main_title, face = "bold", padding = unit(0.1, "line"))

  }

  return(t1)

}
