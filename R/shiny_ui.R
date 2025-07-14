#' Create the Shiny UI for scGenePanel
#'
#' @param data A Seurat object used to populate gene choices
#' @return A Shiny UI object
#' @keywords internal
#' @importFrom shiny dashboardPage dashboardHeader dashboardSidebar dashboardBody
#' @importFrom shiny sidebarMenu menuItem selectizeInput textInput tabItems tabItem
#' @importFrom shiny fluidPage verticalLayout plotOutput tags HTML
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyWidgets addSpinner
#' @noRd
#'

create_genepanel_ui <- function(data) {

  # get filenmae of the object, to be used in metadata display block
  object_name <- attr(data, "object_name") %||% "Unknown Dataset"

  # Get metadata column names, filtering for categorical variables only
  meta_cols_all <- colnames(data@meta.data)
  meta_cols <- meta_cols_all[sapply(data@meta.data, function(x) {
    # Keep factors, characters, and logicals
    if(is.factor(x) || is.character(x) || is.logical(x)) {
      return(TRUE)
    }
    # For integers, check if they have 10 or fewer unique values (likely categorical)
    if(is.integer(x)) {
      return(length(unique(x[!is.na(x)])) <= 10)
    }
    # Exclude all other types (numeric, etc.)
    return(FALSE)
  })]

  # If no categorical columns found, fall back to all columns
  if(length(meta_cols) == 0) {
    meta_cols <- meta_cols_all
  }
  meta_cols <- sort(meta_cols)

  # Identify potential cell type columns (containing "cluster" or "cell")
  potential_celltype_cols <- meta_cols[grepl("cluster|cell", meta_cols, ignore.case = TRUE)]

  # Create formatted display for metadata columns
  create_metadata_display <- function(cols, potential_cols) {
    col_items <- lapply(cols, function(col) {
      if (col %in% potential_cols) {
        tags$li(
          tags$code(col, style = "background-color: #ffffcc; padding: 2px 4px; border-radius: 3px; font-weight: bold;")
        )
      } else {
        tags$li(tags$code(col, style = "background-color: #f8f9fa; padding: 2px 4px; border-radius: 3px;"))
      }
    })
    return(col_items)
  }

  #' Create color palette preview section
  #' @keywords internal
  #' @importFrom shiny tags
  #' @noRd
  create_color_palette_section <- function() {

    # Color palette data
    palettes <- list(
      "Accent" = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F"),
      "Dark2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"),
      "default" = c("#FF0000FF", "#FF9900FF", "#FFCC00FF", "#00FF00FF", "#6699FFFF", "#CC33FFFF"),
      "Paired" = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C"),
      "Pastel1" = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC"),
      "Pastel2" = c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE"),
      "Set1" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"),
      "Set2" = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"),
      "Set3" = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462"),
      "Tableau" = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D")
    )

    # Helper function to create a single palette row
    create_palette_row <- function(name, colors) {
      color_swatches <- lapply(colors, function(color) {
        tags$div(class = "color-swatch", style = paste0("background-color: ", color, ";"))
      })

      tags$div(
        class = "palette-row",
        tags$span(paste0(name, ":"), class = "palette-name"),
        tags$div(class = "palette-swatches", color_swatches)
      )
    }

    # Create all palette rows
    palette_rows <- mapply(create_palette_row, names(palettes), palettes, SIMPLIFY = FALSE)

    return(palette_rows)
  }

  # Header
  header <- shinydashboard::dashboardHeader(
    titleWidth = "100%",
    tags$li(class = "dropdown",
            tags$style(".main-header {max-height: 200px}"),
            tags$style(".main-header .logo {height: 0}")
    )
  )

  anchor <- tags$header(
    tags$div(
      tags$span("",
                style = "color: #FFFFFF; ")
    ),
    style = "height: 20px;"
  )

  header$children[[2]]$children <- tags$div(
    tags$head(tags$style(HTML(".name Gene-label { font-size:80%} "))),
    anchor, class = "name"
  )

  # UI
  ui <- shinydashboard::dashboardPage(
    header,
    title = "scGenePanel - multipanel plots for single cell expression data",
    # Sidebar
    shinydashboard::dashboardSidebar(width = 325,
                                        sidebarMenu(
                                        tags$div(
                                        style = "text-align: center; padding-bottom: 10px;",
                                        tags$a(
                                        href = "#",
                                        onclick = "document.querySelector('a[data-value=\"home\"]').click();",
                                        tags$img(src = "www/logo.svg", height = "150px")
                                        )
                                        ),
                                       id = "tabs",
                                       menuItem("🏡 Home", tabName = "home", selected = TRUE),
                                       tags$hr(),
                                       tags$h3("Step 1 - Configure:", style = "color: white; margin-left: 15px;"),
                                       # Use fluidRow and column for better alignment of inputs
                                       fluidRow(
                                         column(12,
                                                selectizeInput(inputId = "Gene",
                                                               label = HTML("Gene of interest"),
                                                               choices = NULL,
                                                               selected = NULL,
                                                               options = list(
                                                                 placeholder = "Type to search genes...",
                                                                 maxOptions = 50
                                                               ))
                                         )
                                       ),
                                       fluidRow(
                                         column(12,
                                                selectizeInput(inputId = "cell_type_colname",
                                                               label = "Metadata name for cell type annotation",
                                                               choices = meta_cols,
                                                               selected = if("CellTypes" %in% meta_cols) "CellTypes" else meta_cols[1],
                                                               options = list(
                                                                 placeholder = "Select metadata column...",
                                                                 maxOptions = 50
                                                               ))
                                         )
                                       ),

                                       fluidRow(
                                         column(12,
                                                selectizeInput(inputId = "cell_type_name",
                                                               label = "Cell type (# of cells)",
                                                               choices = NULL,  # Will be populated by server
                                                               selected = NULL,
                                                               options = list(
                                                                 placeholder = "Select cell type...",
                                                                 maxOptions = 50
                                                               ))
                                         )
                                       ),
                                       fluidRow(
                                         column(12,
                                                selectizeInput(inputId = "meta_group",
                                                               label = "Group - metadata column name for trait to explore",
                                                               choices = meta_cols,
                                                               selected = if("Source" %in% meta_cols) "Source" else meta_cols[1],
                                                               options = list(
                                                                 placeholder = "Select metadata column...",
                                                                 maxOptions = 50
                                                               ))
                                         )
                                       ),
                                       fluidRow(
                                         column(12,
                                                selectizeInput(inputId = "color",
                                                               label = "Color palette",
                                                                choices = c("Accent","Dark2","default","Paired","Pastel1","Pastel2","Set1","Set2","Set3","Tableau"),
                                                                selected = "Set1")
                                         )
                                       ),
                                       tags$h3("Step 2 - Choose a plot type:", style = "color: white; margin-left: 15px;"),
                                       menuItem("🎯 Multi-Panel", tabName = "fullpanel", badgeLabel = "UMAP + Violin + Table", badgeColor = "fuchsia"),
                                       menuItem("🗺️ UMAP plot", tabName = "umap"),
                                       menuItem("🎻 Violin plot", tabName = "vlnplot"),
                                       menuItem("📊 Table plot", tabName = "table")
                                     )
    ),

    # Body
    shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),

      # Mobile menu button
      tags$button(
        class = "mobile-menu-toggle",
        onclick = "toggleMobileMenu();",
        "☰"
      ),

      # Overlay for mobile
      tags$div(
        class = "sidebar-overlay",
        onclick = "toggleMobileMenu();"
      ),

      # JavaScript for mobile functionality
      tags$script(HTML("
        function toggleMobileMenu() {
          var sidebar = document.querySelector('.main-sidebar');
          var overlay = document.querySelector('.sidebar-overlay');

          if (sidebar.classList.contains('sidebar-open')) {
            sidebar.classList.remove('sidebar-open');
            overlay.classList.remove('show');
          } else {
            sidebar.classList.add('sidebar-open');
            overlay.classList.add('show');
          }
        }

        // Close menu when clicking on menu items
        document.addEventListener('DOMContentLoaded', function() {
          var menuItems = document.querySelectorAll('.sidebar-menu a');
          menuItems.forEach(function(item) {
            item.addEventListener('click', function() {
              if (window.innerWidth <= 768) {
                toggleMobileMenu();
              }
            });
          });
        });
      ")),

      tags$head(
        tags$title("scGenePanel"),
        tags$style(HTML("
          .skin-blue .main-header .logo { background-color: #196797; }
          .content-wrapper, .right-side { background-color: #FFFFFF; }
          .skin-blue .main-sidebar { font-size: 18px; background-color: #2f2f2f; }
          .main-sidebar { font-size: 18px; }
          .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{ background-color: #196797; }
          .left-side, .main-sidebar { padding-top: 10px; }
          .loading-spinner { left:45% !important; }
          header { padding-top:20px 0 0 0 !important; }
          .full-panel-container { padding: 20px; background-color: white; }
          .panel-title { text-align: center; font-size: 24px; font-weight: bold; margin-bottom: 20px; color: #333; }
          code { color: #dd1c77; background-color: #ffffcc; font-family: monospace;}
          .metadata-section { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; border: 1px solid #ddd; }
          /* hide header, adjust sidebar */
          header { display: none;}
          .main-sidebar { padding-top: 10px; }

          /* Mobile menu toggle button */
          .mobile-menu-toggle { display: none; position: fixed; top: 15px; left: 15px; z-index: 1001; background: #2f2f2f; color: white; border: none; padding: 10px 12px; border-radius: 3px; font-size: 18px; cursor: pointer;
          }
          .mobile-menu-toggle:hover { background: #196797; }

          /* Mobile responsiveness */
          @media (max-width: 768px) {
            .mobile-menu-toggle {display: block; }
            .main-sidebar { transform: translateX(-100%); transition: transform 0.3s ease; position: fixed; height: 100vh; z-index: 1000; }
            .main-sidebar.sidebar-open { transform: translateX(0); }
            .content-wrapper { margin-left: 0 !important; }
            /* Overlay for mobile menu */
            .sidebar-overlay { display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.5); z-index: 999;}
            .sidebar-overlay.show { display: block; }
          }
          /* Color Palette Preview Styles */
          .palette-row { margin-bottom: 8px; display: flex; align-items: center; }
          .palette-name { font-weight: bold; margin-right: 10px; display: inline-block; width: 80px; font-size: 14px; }
          .palette-swatches { display: flex; gap: 2px; }
          .color-swatch { width: 30px; height: 20px; border-radius: 2px; transition: transform 0.2s ease;
          }
        "))
      ),

      tabItems(
        tabItem(tabName = "home",
                fluidPage(

                  verticalLayout(

                    tags$h2("Welcome to scGenePanel Interactive Viewer!", style = "float:left"),
                    hr(),
                    tags$p("Quickly explore gene expression data interactively, focusing on specific genes across different cell types and experimental conditions. Generate a publication-ready figure that includes gene expression UMAP projections, violin plots, and summary statistics."),

                    fluidRow(
                      column(6,
                             tags$div(
                               class = "panel-description",
                               tags$ul(
                                 tags$li("🎯 ", tags$strong("Multi-Panel (UMAP + Violin + Table"), " - Complete integrated multi-panel visualization"),
                                 tags$li("🗺️ ", tags$strong("UMAP View"), " - Only ", tags$code("gene"), " expression feature plot of UMAP embedding"),
                                 tags$li("🎻 ", tags$strong("Violin View"), " - Only violin plots to assess ", tags$code("gene"), " expression distribution across ", tags$code("groups"), ""),
                                 tags$li("📊 ", tags$strong("Table View"), " - Only cell frequency and expression statistics of selected ", tags$code("gene"), " in selected ", tags$code("groups"), ""),
                                 tags$li("See ", tags$a(href="https://github.com/vandydata/scGenePanel", "scGenePanel GitHub Page"), " for more details.")
                               ),
                             ),
                             tags$div(
                               tags$h3("How to cite"),
                               tags$p("If you use this app in your research, please cite:"),
                               tags$blockquote(
                                 tags$p("Shrestha, R., et al. (2025). scGenePanel..."),
                                 style = "font-style: italic; margin-left: 20px; border-left: 4px solid #ddd; padding-left: 20px;"
                               )
                             )
                      ),
                      column(6,
                             tags$div(
                               class = "panel-image",
                               tags$a(href = "www/scGenePanel__ATF4_Beta_Age.jpg", target = "_blank",
                                      tags$img(src = "www/scGenePanel__ATF4_Beta_Age.jpg", width = "80%", style = "cursor: pointer; border: 2px solid #ddd; border-radius: 5px; hover: border-color: #007bff;")),
                             )
                      )
                    ),


                    tags$div(
                      tags$h3("Quick start"),
                      tags$ol(
                        tags$li("Select a ", tags$code("gene"), " of interest. ", tags$strong("Type"), " your gene of interet to find it (e.g. INS or GCG)"),
                        tags$li("Select the object's metadata column name for the  ", tags$code("cell type annotation"), "to use. For example, it could be 'celltypes' or 'celltype' - there is no standard but you must choose one that exists in the object."),
                        tags$li("Select the ", tags$code("cell type"), " (e.g., Beta or Alpha) - this must be present in item #2 below"),
                        tags$li("Select a ", tags$code("trait or variable"), " of interest in the object's metadata to split the data by. For example, 'age' or 'source'"),
                        tags$li("Select a color palette (optional)"),
                        tags$li("Click on the plot type desired"),
                      )
                    ),

                    fluidRow(
                      class = "metadata-section",
                      column(6,
                             tags$div(

                               tags$h3("Available metadata columns"),
                             tags$p("We extracted all metadata columns from the loaded object: ", tags$code(object_name), "."),
                             # Show highlighted potential cell type columns if any exist
                             if (length(potential_celltype_cols) > 0) {
                               tags$div(
                                 tags$p(
                                   tags$strong("💡 Candidate cell type columns: "),
                                   "Columns highlighted in yellow are likely cell type annotations. You may want to use one of these for the 'Metadata column name of cell type annotation' field above."
                                 ),
                                 tags$ul(
                                   style = "margin-bottom: 15px;",
                                   lapply(potential_celltype_cols, function(col) {
                                     tags$li(
                                       tags$code(col, style = "background-color: #ffffcc; padding: 4px 8px; border-radius: 3px; font-weight: bold;"),
                                     )
                                   })
                                 )
                               )
                             } else {
                               tags$div()
                             },
                             # Show all metadata columns
                             tags$strong("All Metadata Columns:"),
                             tags$div(
                               class = "metadata-columns",
                               tags$ul(
                                 create_metadata_display(meta_cols, potential_celltype_cols)
                               )
                             )
                           )
                      ),
                      column(6,
                             tags$h3("Available color palettes"),
                             tags$p("Preview of all available color palettes for plots:"),


                             create_color_palette_section()
                             # Accent palette

                      )


                    )
                  )
                )
        ),

        # Full Panel tab
        tabItem(tabName = "fullpanel",
                fluidPage(
                  div(class = "full-panel-container",
                      div(class = "panel-title", textOutput("panelTitle")),
                      shinyWidgets::addSpinner(
                        plotOutput("fullPanel", height = "800px", width = "100%"),
                        spin = "dots", color = "#2b6cb3"
                      ),
                      br(),
                      tags$div(
                        style = "text-align: center; color: #666; font-style: italic;",
                        "Generated by scGenePanel Interactive Viewer"
                      )
                  )
                )
        ),

        # UMAP tab
        tabItem(tabName = "umap",
                fluidPage(
                  verticalLayout(
                    tags$h3("UMAP Visualization", style = "text-align: center;"),
                    br(),
                    shinyWidgets::addSpinner(plotOutput("plot1"), spin = "dots", color = "#2b6cb3")
                  )
                )
        ),

        # Violin tab
        tabItem(tabName = "vlnplot",
                fluidPage(
                  verticalLayout(
                    tags$h3("Violin Plot Visualization", style = "text-align: center;"),
                    br(),
                    shinyWidgets::addSpinner(plotOutput("plot2"), spin = "dots", color = "#2b6cb3")
                  )
                )
        ),

        # Table tab
        tabItem(tabName = "table",
                fluidPage(
                  verticalLayout(
                    tags$h3("Cell Frequency Table", style = "text-align: center;"),
                    br(),
                    shinyWidgets::addSpinner(plotOutput("plot3", height = "600px", width = "100%"), spin = "dots", color = "#2b6cb3")
                  )
                )
        )
      )
    )
  )

  return(ui)
}
