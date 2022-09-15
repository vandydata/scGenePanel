
library(shiny)
library(shinydashboard)
library(shinythemes)
require(shinyjs)
library(shinyBS)
library(shinyWidgets)


# load data----
data <- readRDS("DATA/data.TK.Rds")

# required R functions
source("../../../R/checks.R")
source("../../../R/colorpanels.R")
source("../../../R/genepanel.R")


#1.Header----
header <- dashboardHeader(titleWidth = "100%",
                          # Set height of dashboardHeader
                          tags$li(class = "dropdown",
                                  tags$style(".main-header {max-height: 200px}"),
                                  tags$style(".main-header .logo {height: 100px}")
                          )
)

# webpage links to the images
anchor <- tags$header(
  tags$a(href='https://cds.vanderbilt.edu',
         tags$img(src='CDS-logo-600x85.png', width='200',style="float:left; margin:0 70px 10px 20px; margin-top:15px;" )),
  #style = "padding-top:100px; padding-bottom:100px;"),
  'Genemap Visual for Single Cell RNA-seq Data',
  style = "color: #FFFFFF;
  float:center;
  /*font-family: Avenir Light;*/
  font-size: 35px;
  padding:20px;
  font-weight: bold"
)

header$children[[2]]$children <- tags$div(
  tags$head(tags$style(HTML(".name Gene-label { font-size:80%} "))),
  anchor,
  class = 'name')

#2.User Interface----

#*  Dashboard header----
ui<-dashboardPage( header,
                   title = "Genemap visual for single cell RNAseq data",
                   #skin = "green",

                  #* Dashboard sidebar ----
                  dashboardSidebar(width = 400,

                                   sidebarMenu(
                                     id = "tabs",
                                     menuItem("Home", tabName = "home", selected = T),
                                     menuItem(startExpanded = TRUE,
                                       selectizeInput(inputId = "Gene",
                                                      label = "Enter Official Gene Symbol",
                                                      choices=NULL)),
                                     textInput(inputId = "cell_type_name",label = "Enter name of cell type"),
                                     textInput(inputId = "cell_type_colname",label = "Enter column name of the cell type annotation"),
                                     textInput(inputId = "meta_group",label = "Enter group name"),

                                     menuItem(startExpanded = TRUE,
                                              selectizeInput(inputId = "color",
                                                             label = "Enter color palette",
                                                             choices=NULL)),
                                     menuItem("UMAP", tabName = "umap"),
                                     menuItem("Violinplot", tabName = "vlnplot"),
                                     menuItem("Cell Frequency", tabName = "cellfreq")
                                   )
                  ),

                  #*  Dashboardbody----
                  dashboardBody(
                    useShinyjs(),

                    tags$head(
                      tags$title("Genepanel"),
                      tags$style(HTML('
                                            /* header */
                                            .skin-blue .main-header .logo {
                                            background-color: #79A8B7;
                                            }

                                            /* body */
                                            .content-wrapper, .right-side {
                                            background-color: #FFFFFF;
                                            }

                                            /* main sidebar */
                                            .skin-blue .main-sidebar { font-size: 20px;
                                                            background-color: #333333;
                                            }
                                            .main-sidebar { font-size: 20px; }

                                            /* active selected tab in the sidebarmenu */
                                            .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                            background-color: #006683;
                                            }

                                            .left-side, .main-sidebar {
                                            	padding-top: 110px;
                                            }

                                            /* Gene-label { font-size:70%;} */

                                            /* image {max-width: 60%; width: 60%; height: auto; } */

                                            /* fix for spinner showing up in right of plots in large monitors
                                            not elegant, but quick fix */
                                            .loading-spinner { left:25% !important;}

                                            header { padding-top:20px 0 0 0 !important;}

                                            '))),


                    #mainPanel(tableOutput("table1")), # gene description table

                    tabItems(
                      tabItem(tabName = "home",
                              fluidPage(
                                verticalLayout(tags$h2("Welcome!"),
                                               hr(),
                                               tags$h4("This app provides interactive access to our single cell RNA-Seq data that is reported in:")
                                               #                                      tags$div(
                                               #                                        HTML("<p style = 'font-size:20px;'><u><b>Combinatorial transcription factor profiles predict mature and functional human islet α and β cells.</u></b><br></p> Shristi Shrestha*, Diane C. Saunders*, John T. Walker*, Joan Camunas-Soler, Xiao-Qing Dai, Rachana Haliyur, Radhika Aramandla, Greg Poffenberger, Nripesh Prasad, Rita Bottino, Roland Stein, Jean-Philippe Cartailler, Stephen C. J. Parker, Patrick E. MacDonald, Shawn E. Levy, Alvin C. Powers, Marcela Brissova, <b>JCI Insight. 2021 Sep 22 doi: 10.1172/jci.insight.151621 </b><br> *first co-authors <blockquote style='font-size:15px'> Abstract <br> Islet-enriched transcription factors (TFs) exert broad control over cellular processes in pancreatic α and β cells and changes in their expression are associated with developmental state and diabetes. However, the implications of heterogeneity in TF expression across islet cell populations are not well understood. To define this TF heterogeneity and its consequences for cellular function, we profiled >40,000 cells from normal human islets by scRNA-seq and stratified α and β cells based on combinatorial TF expression. Subpopulations of islet cells co-expressing ARX/MAFB (α cells) and MAFA/MAFB (β cells) exhibited greater expression of key genes related to glucose sensing and hormone secretion relative to subpopulations expressing only one or neither TF. Moreover, all subpopulations were identified in native pancreatic tissue from multiple donors. By Patch-seq, MAFA/MAFB co-expressing β cells showed enhanced electrophysiological activity. Thus, these results indicate combinatorial TF expression in islet α and β cells predicts highly functional, mature subpopulations.</blockquote>
                                               # 	                                                      "))
                                )
                              )

                      ),

                    #For UMAP plot tab
                    tabItem(tabName = "umap",
                            fluidPage(
                              verticalLayout(tableOutput("umap"),
                                             br(),
                                             addSpinner(plotOutput("plot1"), spin = "dots", color = "#2b6cb3")
                              )
                            )
                    ),

                    #For violinplot tab
                    tabItem(tabName = "vlnplot",
                            fluidPage(
                              verticalLayout(tableOutput("vlnplot"),
                                             br(),
                                             addSpinner(plotOutput("plot2"), spin = "dots", color = "#2b6cb3"),
                                               # br(),br(),
                                               # sidebarPanel(
                                               #   sliderInput("Cellsize",
                                               #               "Increase Cell Size:",
                                               #               min=-1,
                                               #               max=1,
                                               #               value=0.1,
                                               #               step = 0.1,
                                               #               animate=TRUE),
                                               #   hr(),
                                               #   helpText("set to -1 to remove dots(cells)")
                                               #   )
                                               )
                                )
                              ),

                      # Cell frequency
                      tabItem(tabName = "cellfreq",
                              fluidPage(
                                verticalLayout(plotOutput("plot3"),
                                )
                              )
                      )

                    )))


#3.Server
server <- function(input, output, session) {

  observe({

    updateSelectizeInput(
      session,
      'Gene',
      choices = rownames(data),
      server = TRUE,
      selected=1)
    })

  observe({

    updateSelectizeInput(
      session,
      'color',
      choices = c("Set1","Set2","Set3","Paste1","Paste2","Paired","Dark2","Accent","tableu","varibow"),
      server = TRUE,
      selected=1)
  })



  # 1. Umap plot ----
  output$plot1<- renderPlot({

    req(input$Gene)
    req(input$cell_type_name)
    req(input$meta_group)
    req(input$cell_type_colname)
    req(input$color)

    umap_panel(seurat_object = data,
                      gene = input$Gene,
                      cell_type_name = input$cell_type_name,
                      meta_group = input$meta_group,
                      cell_type_colname = input$cell_type_colname,
                      col.palette = input$color,
                      output_dir = getwd())

  }, height = 700, width = 600)

  # 2. Vlnplot ----
  output$plot2<- renderPlot({

    req(input$Gene)
    req(input$cell_type_name)
    req(input$meta_group)
    req(input$cell_type_colname)
    req(input$color)

    violin_panel(seurat_object = data,
               gene = input$Gene,
               cell_type_name = input$cell_type_name,
               meta_group = input$meta_group,
               cell_type_colname = input$cell_type_colname,
               col.palette = input$color,
               output_dir = getwd())

  }, height = 600, width = 1200)

  # 3. Cell frequency table ----
  output$plot3<-renderPlot({
    req(input$Gene)
    req(input$cell_type_name)
    req(input$meta_group)
    req(input$cell_type_colname)
    req(input$color)

    cellfreq_panel(seurat_object = data,
                 gene = input$Gene,
                 cell_type_name = input$cell_type_name,
                 meta_group = input$meta_group,
                 cell_type_colname = input$cell_type_colname,
                 output_dir = getwd())

  }, height = 1000, width = 1200)

}

shinyApp(ui=ui,server=server)
