
#' Discrete color palettes from multiple packages (varibow, Rcolorbrewer,
#' ggtheme)
#' @param num_colors Number of colors to be generated.
#' @param hex_cols hex code of the colors to be selected for violin plot
#' @param palette Options are "Tableau" or RColorBrewer qualitative variables
#' like "Dark2","Paired","Set1","Set2", "Set3", "Accent" etc.
#' if empty, the function will use based on the entered n.
#' @return A color palette for plotting violin plot

#select_col <-  discrete_col_palette(num_colors = 5, palette = "Dark2", hex_cols = FALSE)

discrete_col_palette <- function(num_colors, hex_cols = NULL,  palette = NULL) {

  pal <- ggthemes::tableau_color_pal("Tableau 20")

  col_list <- list(
    default = ggsci::pal_ucscgb()(num_colors),
    Tableau = pal(attr(pal, "max_n")),
    Dark2 = RColorBrewer::brewer.pal(8, "Dark2"),
    Paired = RColorBrewer::brewer.pal(12, "Paired"),
    Set1 = RColorBrewer::brewer.pal(9, "Set1"),
    Set2 = RColorBrewer::brewer.pal(8, "Set2"),
    Set3 = RColorBrewer::brewer.pal(12, "Set3"),
    Pastel1 = RColorBrewer::brewer.pal(9, "Pastel1"),
    Pastel2 = RColorBrewer::brewer.pal(8, "Pastel2"),
    Accent = RColorBrewer::brewer.pal(8, "Accent")

  )

  # Empty palette input? Set to default
  if(is.null(palette)){
    palette <- "default"
  }

  # Incorrect (unknown) color palette input? Set to default, warn user
  if(palette %in% names(col_list)){
    palette <- palette
  } else {
    palette <- "default"
    palette_str <- paste0(names(col_list), collapse = ", ")
    message <- paste0("You have provided an unknown palette so defaulting to pal_ucscgb. Please choose from: ", palette_str, ".")
    warning(message)
  }

  selected_palette <- col_list[[palette]]

  if (num_colors > length(selected_palette)) {
    warning("Not enough colors in specified palette. Choose another or enter NULL (default to ggsci UCSCB Genome Browser palette).")
  }

  if (is.null(x = palette) && !is.null(x = hex_cols)) {

    selected_palette <- hex_cols
    if (num_colors > length(selected_palette)) {
      warning("Not enough colors in specified hex_cols.")
    }
  }
  selected_palette

}
