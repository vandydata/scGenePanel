#' Calculate cell-frequency and expression summaries for a selected cell type
#'
#' This separates the deterministic summary calculation from its flextable
#' rendering, allowing the values displayed in the table to be unit tested.
#'
#' @return A data frame used to construct the cell-frequency table.
#' @keywords internal
#' @noRd
.cellfreq_summary <- function(seurat_obj,
                              cell_type_colname,
                              cell_type_name,
                              meta_group,
                              gene) {
  Seurat::Idents(seurat_obj) <- cell_type_colname
  selected_cluster_cells <- subset(seurat_obj, idents = cell_type_name)

  group_values <- selected_cluster_cells@meta.data[[meta_group]]
  group_levels <- if (is.factor(group_values)) {
    levels(droplevels(group_values))
  } else {
    unique(as.character(group_values))
  }
  expression <- as.numeric(
    Seurat::GetAssayData(selected_cluster_cells[["RNA"]])[gene, ]
  )

  summary_table <- do.call(rbind, lapply(group_levels, function(group) {
    group_expression <- expression[as.character(group_values) == group]
    data.frame(
      group = group,
      n_cells = length(group_expression),
      n_expressing = sum(group_expression > 0.25),
      pct_expressed = round(mean(group_expression > 0.25) * 100, 1),
      `25%` = round(stats::quantile(group_expression, 0.25), 2),
      `50%` = round(stats::quantile(group_expression, 0.50), 2),
      `75%` = round(stats::quantile(group_expression, 0.75), 2),
      `100%` = round(stats::quantile(group_expression, 1), 2),
      check.names = FALSE
    )
  }))
  names(summary_table)[1] <- meta_group
  names(summary_table)[3] <- "n_expressing (>0.25)"
  names(summary_table)[4] <- "% expressed"
  summary_table
}

#' Generate a cell-frequency table for scGenePanel
#'
#' @param seurat_obj A Seurat object.
#' @param cell_type_colname Metadata column containing cell-type annotations.
#' @param cell_type_name Cell type to summarize.
#' @param meta_group Metadata column defining groups.
#' @param gene Gene to summarize.
#' @param col_palette Unused; retained for compatibility with the panel API.
#' @return A grob containing cell counts and expression summaries.
#' @keywords internal
#' @noRd
cellfreq_panel <- function(seurat_obj,
                           cell_type_colname,
                           cell_type_name,
                           meta_group,
                           gene,
                           col_palette) {
  combined_cc_table <- .cellfreq_summary(
    seurat_obj = seurat_obj,
    cell_type_colname = cell_type_colname,
    cell_type_name = cell_type_name,
    meta_group = meta_group,
    gene = gene
  )

  if (all(combined_cc_table$`n_expressing (>0.25)` == 0)) {
    message(sprintf("No expression detected for: %s", gene))
  }

  expression_quantiles <- combined_cc_table[, c("25%", "50%", "75%", "100%")]
  colourer <- scales::col_numeric(
    palette = c("#fcee98", "#da4362"),
    domain = range(expression_quantiles)
  )

  flextable::set_flextable_defaults(
    font.color = "#000000",
    border.color = "#cccccc"
  )
  ft <- flextable::flextable(combined_cc_table)
  ft <- flextable::add_header_row(
    ft,
    colwidths = c(1, 3, 4),
    values = c("", "Gene Expression", "Quantiles")
  )
  ft <- flextable::theme_vanilla(ft)
  ft <- flextable::font(ft, part = "all", fontname = "mono")

  border <- flextable::fp_border_default()
  ft <- flextable::vline(
    ft,
    j = c(meta_group, "% expressed"),
    border = border,
    part = "all"
  )
  ft <- flextable::footnote(
    ft,
    i = 2,
    j = 2:4,
    value = flextable::as_paragraph(c(
      "n_cells - Number of cells in group",
      "n_expressing (>0.25%) - Number of cells with expression value at least 0.25",
      "% expressed = (n_expressing / n_cells) * 100"
    )),
    ref_symbols = c("a", "b", "c"),
    part = "header",
    inline = TRUE
  )
  ft <- flextable::footnote(
    ft,
    i = 1,
    j = 5,
    value = flextable::as_paragraph(
      "Gene expression values per quantile for assessing distribution within each group."
    ),
    ref_symbols = "d",
    part = "header",
    inline = TRUE
  )
  ft <- flextable::bg(
    ft,
    j = c("25%", "50%", "75%", "100%"),
    bg = colourer,
    part = "body"
  )
  ft <- flextable::fontsize(ft, size = 20, part = "header")
  ft <- flextable::fontsize(ft, size = 18, part = "body")
  ft <- flextable::fontsize(ft, size = 18, part = "footer")

  flextable::gen_grob(ft)
}
