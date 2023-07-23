
#' Discrete color paletters from multiple packages (varibow, Rcolorbrewer,
#' ggtheme)
#' @param num_colors Number of colors to be generated.
#' @param hex_cols hex code of the colors to be selected for violin plot
#' @param palette Options are "tableu" or RColorBrewer qualitative variables
#' like "Dark2","Paired","Set1","Set2", "Set3", "Accent" etc.
#' if empty, the function will use based on the entered n.
#' @return A color palette for plotting violin plot

#select_col <-  discrete_col_palette(num_colors = 5, palette = "Dark2", hex_cols = FALSE)

discrete_col_palette <- function(num_colors, hex_cols = NULL,  palette = NULL) {

  pal <- ggthemes::tableau_color_pal("Tableau 20")

  col_list <- list(
    Tableau = pal(attr(pal, "max_n")),
    Dark2 = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(num_colors),
    Paired = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(num_colors),
    Set1 = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_colors),
    Set2 = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(num_colors),
    Set3 = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(num_colors),
    Pastel1 = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))(num_colors),
    Pastel2 = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel2"))(num_colors),
    Accent = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(num_colors)

  )

  selected_palette <- col_list[[palette]]

  if (num_colors > length(selected_palette)) {
    warning("Not enough colors in specified palette. Choose another or enter NULL which will default to varibow")
  }



  if (is.null(x = palette) && !is.null(x = hex_cols)) {

    selected_palette <- hex_cols
    if (num_colors > length(selected_palette)) {
      warning("Not enough colors in specified hex_cols.")
    }
    selected_palette # TODO - do you really need this?
  }
  selected_palette

}
