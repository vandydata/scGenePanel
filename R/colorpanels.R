
#' Discrete color paletters from multiple packages (varibow,Rcolorbrewer,ggtheme)
#' @param num_colors Number of colors to be generated.
#' @param hex_cols hex code of the colors to be selected for violin plot
#' @param palette Options are "tableu","varibow" or RColorBrewer qualitative variables like "Dark2","Paired","Set1" etc
#' if empty, the function will use based on the entered n.
#' @return A color palette for plotting violin plot

#select_col <-  discrete_col_palette(num_colors = 5, palette = "Dark2", hex_cols = FALSE)

discrete_col_palette <- function(
  num_colors,
  hex_cols = NULL,
  palette = NULL

  ) {

  pal <- ggthemes::tableau_color_pal("Tableau 20")

  col_list <- list(
    tableau = pal(attr(pal, "max_n")),
    varibow = scrattch::varibow(n_colors = num_colors),
    # only qualitative palettes suitable for categorical data from Rcolorbrewer
    Dark2 = RColorBrewer::brewer.pal(num_colors, 'Dark2'),
    Paired = RColorBrewer::brewer.pal(num_colors, 'Paired'),
    Set1 = RColorBrewer::brewer.pal(num_colors, 'Set1'),
    Set2 = RColorBrewer::brewer.pal(num_colors, 'Set2'),
    Set3 = RColorBrewer::brewer.pal(num_colors, 'Set3'),
    Pastel1 = RColorBrewer::brewer.pal(num_colors, 'Pastel1'),
    Pastel2 = RColorBrewer::brewer.pal(num_colors, 'Pastel2'),
    Paired = RColorBrewer::brewer.pal(num_colors, 'Paired'),
    Accent = RColorBrewer::brewer.pal(num_colors, 'Accent')
  )

  if (is.null(x = palette) && is.null(x = hex_cols)) {
    if (num_colors <= 20) {
      palette <- "tableu"
    } else {
      palette <- "varibow"
    }
  }

  selected_palette <- col_list[[palette]]

  if (num_colors > length(selected_palette)) {
    warning("Not enough colors in specified palette. Choose another or enter NULL which will default to varibow")
  }



  if(is.null(x = palette) && !is.null(x = hex_cols)){

    selected_palette <- hex_cols
  if (num_colors > length(selected_palette)) {
    warning("Not enough colors in specified hex_cols.")
    }
    selected_palette
  }
  selected_palette

}
