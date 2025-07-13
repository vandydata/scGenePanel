
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

create_genepanel_ui <- function(data) {

  # Header
  header <- shinydashboard::dashboardHeader(
    titleWidth = "100%",
    tags$li(class = "dropdown",
            tags$style(".main-header {max-height: 200px}"),
            tags$style(".main-header .logo {height: 100px}")
    )
  )

  anchor <- tags$header(
    "scGenePanel : Visuals for Single Cell RNA-seq Data",
    style = "color: #FFFFFF; float:center; font-size: 35px; padding:20px; font-weight: bold"
  )

  header$children[[2]]$children <- tags$div(
    tags$head(tags$style(HTML(".name Gene-label { font-size:80%} "))),
    anchor, class = "name"
  )

  # UI
  ui <- shinydashboard::dashboardPage(
    header,
    title = "scGenePanel : Visuals for Single Cell RNA-seq Data",

    # Sidebar
    shinydashboard::dashboardSidebar(width = 400,
      sidebarMenu(
        id = "tabs",
        menuItem("Home", tabName = "home", selected = T),
        menuItem(startExpanded = TRUE,
          selectizeInput(inputId = "Gene",
                         label = "Enter Official Gene Symbol",
                         choices = rownames(data),
                         selected = "INS",
                         options = list(maxOptions = 100))),
        textInput(inputId = "cell_type_name", label = "Enter name of cell type", value = "Beta"),
        textInput(inputId = "cell_type_colname", label = "Enter column name of the cell type annotation", value = "CellTypes"),
        textInput(inputId = "meta_group", label = "Enter group name", value = "Source"),
        menuItem(startExpanded = TRUE,
                 selectizeInput(inputId = "color",
                                label = "Enter color palette",
                                choices = c("Set1","Set2","Set3","Paired","Dark2","Accent"),
                                selected = "Set3")),
        br(),
        tags$hr(),
        tags$h4("Views:", style = "color: white; margin-left: 15px;"),
        menuItem("🏠 Home", tabName = "home"),
        menuItem("🎯 Full Panel", tabName = "fullpanel", badgeLabel = "MAIN", badgeColor = "green"),
        menuItem("🗺️ UMAP Only", tabName = "umap"),
        menuItem("🎻 Violin Only", tabName = "vlnplot"),
        menuItem("📊 Table Only", tabName = "table")
      )
    ),

    # Body
    shinydashboard::dashboardBody(
      shinyjs::useShinyjs(),
      tags$head(
        tags$title("scGenePanel"),
        tags$style(HTML("
          .skin-blue .main-header .logo { background-color: #79A8B7; }
          .content-wrapper, .right-side { background-color: #FFFFFF; }
          .skin-blue .main-sidebar { font-size: 18px; background-color: #333333; }
          .main-sidebar { font-size: 18px; }
          .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{ background-color: #006683; }
          .left-side, .main-sidebar { padding-top: 110px; }
          .loading-spinner { left:45% !important; }
          header { padding-top:20px 0 0 0 !important; }
          .full-panel-container { padding: 20px; background-color: white; }
          .panel-title { text-align: center; font-size: 24px; font-weight: bold; margin-bottom: 20px; color: #333; }
        "))
      ),

      tabItems(
        # Home tab
        tabItem(tabName = "home",
          fluidPage(
            verticalLayout(
              tags$h2("Welcome to scGenePanel Interactive Viewer!"),
              hr(),
              tags$h4("This app provides interactive access to single cell RNA-Seq data visualization using the original scGenePanel functions."),
              tags$p("Features:"),
              tags$ul(
                tags$li("🎯 ", tags$strong("Full Panel View"), " - Complete integrated visualization (UMAP + Violin + Table)"),
                tags$li("🗺️ ", tags$strong("UMAP View"), " - Dimensionality reduction plots"),
                tags$li("🎻 ", tags$strong("Violin View"), " - Expression distribution plots"),
                tags$li("📊 ", tags$strong("Table View"), " - Cell frequency and expression statistics")
              ),
              tags$br(),
              tags$div(
                style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-top: 20px;",
                tags$h5("Quick Start:"),
                tags$p("1. Select a gene (e.g., INS, GCG, SST)"),
                tags$p("2. Choose cell type (e.g., Beta, Alpha, Delta)"),
                tags$p("3. Click on ", tags$strong("Full Panel"), " tab for complete visualization")
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
