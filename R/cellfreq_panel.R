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
                           col_palette
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
    check <- Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ] > 0.25 # TODO parametize this, and a default
    selected_cells_expressed <- selected_cluster_cells[, check]
  }, error = function(e) {
    selected_cells_expressed <- NULL
  })

  if (is.null(selected_cells_expressed)) {
    message(paste0("No expression detected for: ", gene))
  } else {

    selected_cells_expressed_tally <- selected_cells_expressed@meta.data %>%
      dplyr::group_by_at(meta_group) %>%
      dplyr::tally()     # tally cells expressing the gene

    # Merge cell frequency and expression tally
    cc_table <- merge(cell_counts_tally,
                      selected_cells_expressed_tally,
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
    names(meta_subset) <- c("celltype", "group") # TODO - I don't understand why this is necessary after setting the above
    meta_subset <- meta_subset %>% dplyr::filter(celltype == cell_type_name)

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
    result <- round(t(result)[, 2:5], 2)


    combined_cc_table <- cbind(cc_table, result)
    names(combined_cc_table) <- c(
      meta_group,
      "n_cells",
      "n_expressing (>0.25)",
      "% expressed",
      "25%",
      "50%",
      "75%",
      "100%"
    )


    # round to one decimal point
    combined_cc_table$`% expressed` <- round((combined_cc_table$`% expressed` * 100), 1)

    if(TRUE){

      subset <- combined_cc_table[, c( "25%", "50%", "75%", "100%")]
      subset_min <- min(subset)
      subset_max <- max(subset)


      set_flextable_defaults(
        #font.family = "Consolas",
        font.color = "#000000",
        border.color = "#cccccc"
      )

      ft <- flextable::flextable(combined_cc_table)
      ft <- add_header_row(ft,
                           colwidths = c(1, 3, 4),
                           values = c("", "Gene Expression", "Quantiles")
      )
      ft <- theme_vanilla(ft)
      ft <- font(ft, part = "all", fontname = "Consolas")

      border <- fp_border_default()

      ft <- vline(ft, j = c(meta_group, '% expressed'), border = border, part = "all")

      ft <- footnote(ft,
               i = 2, j = 2:4,
               value = as_paragraph(
                 c(
                   "n_cells - Number of cells in group",
                   "n_expressing (>0.25%) - Number of cells with expression value at least 0.25",
                   "% expressed = (n_expressing / n_cells) * 100"
                 )
               ),
               ref_symbols = c("a", "b", "c"),
               part = "header",
               inline = TRUE
      )

      ft <- footnote(ft,
                     i = 1, j = 5,
                     value = as_paragraph(
                       c(
                         "Gene expression values per quantile for assessing distribution within each group."
                       )
                     ),
                     ref_symbols = c("d"),
                     part = "header",
                     inline = TRUE
      )

      colourer <- scales::col_numeric(
        palette = c("#fcee98", "#da4362"),
        domain = c(subset_min, subset_max)
      )
      ft <- bg(ft,
                 j = c(
                   "25%", "50%", "75%", "100%"
                 ),
                 bg = colourer, part = "body"
      )


      ft <- fontsize(ft, size = 20, part = "header")
      ft <- fontsize(ft, size = 18, part = "body")
      ft <- fontsize(ft, size = 18, part = "footer")


      #ft <- set_caption(ft, caption = "Title goes here")
      t1 <- gen_grob(ft)

    }


  }

  return(t1)

}
