#' Create the Shiny server for scGenePanel
#'
#' @param data A Seurat object containing the single-cell data
#' @return A Shiny server function
#' @keywords internal
#' @importFrom shiny renderText renderPlot req updateSelectizeInput observe
#' @importFrom ggplot2 ggplot theme_void annotation_custom xlim ylim theme element_rect margin annotate
#' @importFrom ggpubr ggarrange annotate_figure text_grob rremove
#' @importFrom grid annotation_custom
#' @noRd

create_genepanel_server <- function(data) {

  server <- function(input, output, session) {

    addResourcePath("www", file.path(getwd(), "www"))

    # Fast gene search
    updateSelectizeInput(
      session = session,
      inputId = "Gene",
      choices = sort(rownames(data)),
      selected = "INS",
      server = TRUE
    )

    # Dynamic title for full panel
    output$panelTitle <- renderText({
      req(input$Gene, input$cell_type_name)
      paste0(input$Gene, " expression in ", input$cell_type_name, " cells")
    })


    # Reactive expression to get unique cell types from selected column
    get_cell_types <- reactive({
      req(input$cell_type_colname)

      # Get the selected metadata column
      if(input$cell_type_colname %in% colnames(data@meta.data)) {
        cell_types <- unique(data@meta.data[[input$cell_type_colname]])
        # Remove NA values and sort alphabetically
        cell_types <- sort(cell_types[!is.na(cell_types)])
        return(cell_types)
      }
      return(NULL)
    })

    # Reactive value to track if dropdowns are synchronized
    values <- reactiveValues(dropdowns_synced = TRUE)

    # Observer to update cell type choices when metadata column changes
    observeEvent(input$cell_type_colname, {
      # Mark dropdowns as not synced
      values$dropdowns_synced <- FALSE

      cell_types <- get_cell_types()

      if(!is.null(cell_types)) {
        # Set default to "Beta" if it exists, otherwise first option
        default_selection <- if("Beta" %in% cell_types) "Beta" else cell_types[1]

        updateSelectizeInput(
          session,
          "cell_type_name",
          choices = cell_types,
          selected = default_selection
        )
      }
    })

    # Initialize cell type dropdown on app start
    observe({
      # Only run once when app starts
      if(is.null(input$cell_type_name)) {
        cell_types <- get_cell_types()

        if(!is.null(cell_types)) {
          default_selection <- if("Beta" %in% cell_types) "Beta" else cell_types[1]

          updateSelectizeInput(
            session,
            "cell_type_name",
            choices = cell_types,
            selected = default_selection
          )
        }
      }
    })

    # Observer to mark dropdowns as synced when cell type is updated
    observeEvent(input$cell_type_name, {
      if(!is.null(input$cell_type_name) && input$cell_type_name != "") {
        values$dropdowns_synced <- TRUE
      }
    })

    # Validation function to check if inputs are ready for plotting
    inputs_valid <- reactive({
      req(input$cell_type_colname, input$cell_type_name)

      # Check if dropdowns are synced
      if(!values$dropdowns_synced) {
        return(FALSE)
      }

      # Check if selected cell type exists in the selected metadata column
      if(input$cell_type_colname %in% colnames(data@meta.data)) {
        available_types <- unique(data@meta.data[[input$cell_type_colname]])
        available_types <- available_types[!is.na(available_types)]
        return(input$cell_type_name %in% as.character(available_types))
      }

      return(FALSE)
    })

    # Initialize cell type dropdown on app start
    observe({
      # Only run once when app starts
      if(is.null(input$cell_type_name)) {
        cell_types <- get_cell_types()

        if(!is.null(cell_types)) {
          default_selection <- if("Beta" %in% cell_types) "Beta" else cell_types[1]

          updateSelectizeInput(
            session,
            "cell_type_name",
            choices = cell_types,
            selected = default_selection
          )
        }
      }
    })

    # Full Panel Integration
    output$fullPanel <- renderPlot({
      req(input$Gene, input$cell_type_name, input$meta_group, input$cell_type_colname, input$color)
      req(inputs_valid())

      tryCatch({
        # Get the levels for the metadata group
        levels_idents <- unique(data[[input$meta_group]][, 1])
        levels_idents <- as.character(levels_idents)

        message("Generating UMAP panel...")
        umap <- umap_panel(seurat_obj = data,
                           cell_type_colname = input$cell_type_colname,
                           cell_type_name = input$cell_type_name,
                           meta_group = input$meta_group,
                           gene = input$Gene,
                           levels_idents = levels_idents)

        message("Generating violin panel...")
        violin <- violin_sig_panel(seurat_obj = data,
                                   cell_type_colname = input$cell_type_colname,
                                   cell_type_name = input$cell_type_name,
                                   meta_group = input$meta_group,
                                   gene = input$Gene,
                                   col_palette = input$color,
                                   levels_idents = levels_idents)

        message("Generating cell frequency table...")
        table <- cellfreq_panel(seurat_obj = data,
                                cell_type_colname = input$cell_type_colname,
                                cell_type_name = input$cell_type_name,
                                meta_group = input$meta_group,
                                gene = input$Gene,
                                col_palette = input$color)

        message("Combining panels...")

        figure <- suppressWarnings(
          ggpubr::ggarrange(
            umap,
            NULL,
            violin,
            NULL,
            ggpubr::ggarrange(NULL, table, NULL, ncol=3, widths = c(1, 18, 1)),
            heights = c(1.8, 0.1, 3, 0.1, 1.5),
            ncol = 1,
            nrow = 5)
        )

        figure <- ggpubr::annotate_figure(figure,
                                          top = ggpubr::text_grob(
                                            label = paste0(input$Gene, " expression in ", input$cell_type_name, " cells"),
                                            face = "bold", size = 20
                                          ),
                                          bottom = ggpubr::text_grob(
                                            label = "Generated by scGenePanel Interactive Viewer",
                                            hjust = 1,
                                            x = 1,
                                            face = "italic",
                                            size = 12)
        )

        figure <- figure + ggpubr::rremove("grid")

        message("Full panel generated successfully!")
        return(figure)

      }, error = function(e) {
        message("Error in full panel: ", e$message)
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                   label = paste("Error generating full panel:", e$message),
                   size = 5, color = "red")
      })
    }, height = 800, width = 1000)

    # Individual UMAP plot
    output$plot1 <- renderPlot({
      req(input$Gene, input$cell_type_name, input$meta_group, input$cell_type_colname, input$color)
      req(inputs_valid())

      tryCatch({
        levels_idents <- unique(data[[input$meta_group]][, 1])
        levels_idents <- as.character(levels_idents)

        umap_panel(seurat_obj = data,
                   cell_type_colname = input$cell_type_colname,
                   cell_type_name = input$cell_type_name,
                   meta_group = input$meta_group,
                   gene = input$Gene,
                   levels_idents = levels_idents)
      }, error = function(e) {
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                   label = paste("Error:", e$message),
                   size = 5, color = "red")
      })
    }, height = 600, width = 800)

    # Individual Violin plot
    output$plot2 <- renderPlot({
      req(input$Gene, input$cell_type_name, input$meta_group, input$cell_type_colname, input$color)
      req(inputs_valid())

      tryCatch({
        levels_idents <- unique(data[[input$meta_group]][, 1])
        levels_idents <- as.character(levels_idents)

        violin_sig_panel(seurat_obj = data,
                         cell_type_colname = input$cell_type_colname,
                         cell_type_name = input$cell_type_name,
                         meta_group = input$meta_group,
                         gene = input$Gene,
                         col_palette = input$color,
                         levels_idents = levels_idents)
      }, error = function(e) {
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                   label = paste("Error:", e$message),
                   size = 5, color = "red")
      })
    }, height = 500, width = 1000)

    # Individual Table
    output$plot3 <- renderPlot({
      req(input$Gene, input$cell_type_name, input$meta_group, input$cell_type_colname, input$color)
      req(inputs_valid())

      tryCatch({
        table_result <- cellfreq_panel(seurat_obj = data,
                                      cell_type_colname = input$cell_type_colname,
                                      cell_type_name = input$cell_type_name,
                                      meta_group = input$meta_group,
                                      gene = input$Gene,
                                      col_palette = input$color)

        # Convert grob to a plot that Shiny can render
        if (!is.null(table_result) && inherits(table_result, "grob")) {
          p <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::theme(
              plot.margin = ggplot2::margin(20, 20, 20, 20),
              panel.background = ggplot2::element_rect(fill = "white", color = NA)
            ) +
            ggplot2::annotation_custom(table_result, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
            ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)  # explicit limits for full stretching
          return(p)
        } else {
          ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
                     label = "No expression data found for this combination",
                     size = 6, color = "gray50")
        }

      }, error = function(e) {
        # Error plot
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                   label = paste("Error:", e$message),
                   size = 5, color = "red")
      })
    }, height = 600)
  }

  return(server)
}
