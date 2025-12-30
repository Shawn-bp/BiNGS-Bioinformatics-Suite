# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DESeq2
BiocManager::install("DESeq2")

# Install SummarizedExperiment
BiocManager::install("SummarizedExperiment")

require(shiny)
require(shinythemes)
require(data.table)
require(ggplot2)
require(dplyr)
require(DT)
require(plotly)
require(fst)
require(ggpubr)
require(matrixStats)
require(heatmaply)
library(RColorBrewer)
library(DESeq2)
library(SummarizedExperiment)
source("helper.R")
source("global.R")

ui <- fluidPage(
  titlePanel("BiNGS Transcriptomics Suite"),
  
  tabsetPanel(
    
    # ------------------ UPLOAD TAB ------------------
    tabPanel("Upload",
             fluidRow(
               column(6,
                      br(),
                      fileInput("counts_csv", "Select gene counts file to import", accept = ".csv"),
                      fileInput("metadata_csv", "Select metadata file to import", accept = ".csv")
               ),
               column(6,
                      div(
                        style = "text-align: center; 
                             padding: 60px 40px; 
                             background-color: #337ab7; 
                             border-radius: 5px;
                             margin: 20px;",
                        img(src = "bings_logo_white.png", height = "80px")
                      )
               )
             )
    ),
    
    # ------------------ PCA TAB ------------------
    tabPanel("PCA",
             fluidRow(
               column(3,
                      br(),
                      selectizeInput(
                        'pca_samples_to_remove',
                        label = "Select samples to remove:",
                        choices = NULL,
                        selected = NULL,
                        options = list(
                          `actions-box` = TRUE,
                          placeholder = "Select samples to exclude"
                        ),
                        multiple = TRUE
                      ),
                      selectizeInput(
                        'x_pc',
                        label = "X axis PC",
                        choices = NULL,
                        selected = NULL, 
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      selectizeInput(
                        'y_pc',
                        label = "Y axis PC",
                        choices = NULL,
                        selected = NULL, 
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      selectizeInput(
                        "color_var",
                        label = "Color samples by:",
                        choices = NULL,
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      selectizeInput(
                        "shape_var",
                        label = "Change sample shapes by:",
                        choices = NULL,
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      radioButtons("pca_color_palette",
                                   "Select the color palette to use:",
                                   choices = as.list(pca_color_palette_list),
                                   selected = pca_color_palette_list[[1]]),
                      radioButtons("PCA_Plot_type",
                                   "Interactive or static:",
                                   choices = as.list(type_list),
                                   selected = type_list[[1]]),
                      actionButton("PCA_run_button", "Create PCA", 
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                      br(), 
                      br(), 
                      br(),
                      downloadButton("download_PCA", "Download_PCA_coords.csv")
               ),
               # ------------------ PCA outputs ------------------
               column(9,
                      br(),
                      br(),
                      conditionalPanel("input.PCA_Plot_type == 'plotly'",
                                       plotlyOutput("pca_plotly", height = "600px")),
                      conditionalPanel("input.PCA_Plot_type == 'ggplot'",
                                       plotOutput("pca_ggplot", height = "600px"))
               ),
             )
    ),
    # ------------------ BOXPLOT TAB ------------------
    tabPanel("Boxplot",
             fluidRow(
               column(3,
                      br(),
                      selectizeInput(
                        'boxplot_samples_to_remove',
                        label = "Select samples to remove:",
                        choices = NULL,
                        selected = NULL,
                        options = list(
                          `actions-box` = TRUE,
                          placeholder = "Select samples to exclude"
                        ),
                        multiple = TRUE
                      ),
                      selectizeInput(
                        'gene_var',
                        label = "Select Gene:",
                        choices = NULL,
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      selectizeInput(
                        'fact_var',
                        label = "Sample groups",
                        choices = NULL,
                        selected = NULL, 
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      radioButtons("boxplot_color_palette",
                                   "Select the color palette to use:",
                                   choices = as.list(color_palette_list),
                                   selected = color_palette_list[[1]]),
                      radioButtons("QC_check",
                                   "Remove low quality samples:",
                                   choices = as.list(QC_list),
                                   selected = QC_list[[1]]),
                      radioButtons("Log",
                                   "Log expression:",
                                   choices = as.list(Log_list),
                                   selected = Log_list[[1]]),
                      radioButtons("Boxplot_Plot_type",
                                   "Interactive or static:",
                                   choices = as.list(type_list),
                                   selected = type_list[[1]]),
                      radioButtons("Box_violin",
                                   "Plot type:",
                                   choices = as.list(box_plot_list),
                                   selected = box_plot_list[[1]]),
                      br(),
                      actionButton("run_button", "Create Boxplot", 
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                      br(), 
                      br(), 
                      br(),
                      downloadButton("download_exp", "Download_expression.csv"),
                      br(), 
                      br(),
                      downloadButton("download_pval", "Download_pvalue_table.csv"),
               ),
               # ------------------ BOXPLOT OUTPUTS ------------------
               column(9,
                      br(),
                      br(),
                      conditionalPanel(condition = "input.Boxplot_Plot_type == 'plotly'",
                                       plotlyOutput("boxplot_plotly")),
                      conditionalPanel(condition = "input.Boxplot_Plot_type == 'ggplot'",
                                       plotOutput("boxplot_ggplot")),
                      h3("P-value Table"),
                      DT::dataTableOutput("table_p"),
                      h3("ANOVA Table"),
                      DT::dataTableOutput("table_panova"),
                      h3("Gene Expression Table"),
                      DT::dataTableOutput("table_data")
               ),
             )
    ),
    # ------------------ SAMPLE DISTANCE HEATMAP TAB ------------------
    tabPanel("Sample Distance Heatmap",
             fluidRow(
               column(3,
                      br(),
                      selectizeInput(
                        'sample_distance_samples_to_remove',
                        label= "Select samples to remove:",
                        choices= NULL,
                        selected = NULL,
                        options = list(
                          `actions-box` = TRUE,
                          placeholder= "Select samples to exclude"
                        ),
                        multiple = TRUE
                      ),
                      selectizeInput(
                        'metadata_color_bars',
                        label = "Select metadata field to color by:",
                        choices = NULL,
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      radioButtons("sample_distance_heatmap_show_names",
                                   "Show sample names:",
                                   choices = list("None" = "none",
                                                  "X-axis only" = "x",
                                                  "Y-axis only" = "y",
                                                  "Both" = "both"),
                                   selected = "both"),
                      radioButtons("sample_distance_heatmap_color_palette",
                                   "Select the color palette to use:",
                                   choices = as.list(heatmap_color_scheme_list),
                                   selected = heatmap_color_scheme_list[[1]]),
                      radioButtons("sample_distance_heatmap_scaling_type",
                                   "Select the type of scaling to perform:",
                                   choices = as.list(scaling_list),
                                   selected = scaling_list[[2]]),
                      radioButtons("sample_distance_heatmap_dendrogram_list",
                                   "Select the type of clustering to perform:",
                                   choices = as.list(dendrogram_list),
                                   selected = dendrogram_list[[4]]),
                      radioButtons("Sample_Distance_Heatmap_Plot_type",
                                   "Interactive or static:",
                                   choices = as.list(heatmap_type_list),
                                   selected = heatmap_type_list[[1]]),
                      br(),
                      actionButton("sample_distance_run_button", "Create Heatmap", 
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                      br(),
                      br(),
                      downloadButton("download_distance_matrix", "Download_distance_matrix.csv")
               ),
               column(9,
                      br(),
                      br(),
                      conditionalPanel(
                        "input.Sample_Distance_Heatmap_Plot_type == 'ggplot'",
                        plotOutput("sample_distance_heatmap", height = "700px")
                      ),
                      conditionalPanel(
                        "input.Sample_Distance_Heatmap_Plot_type == 'heatmaply'",
                        plotlyOutput("sample_distance_heatmaply", height = "700px")
                      ))
             )
    ),
    
    # ------------------ GENE EXPRESSION HEATMAP TAB ------------------
    tabPanel("Gene Expression Heatmap",
             fluidRow(
               column(3,
                      br(),
                      selectizeInput(
                        'samples_to_remove_select',
                        label = "Select samples to remove:",
                        choices = NULL,
                        selected = NULL,
                        options = list(
                          `actions-box` = TRUE,
                          placeholder = "Select samples to exclude"
                        ),
                        multiple = TRUE
                      ),
                      selectizeInput(
                        'gene_list_select',
                        label = "Select genes:",
                        choices = NULL,
                        selected = NULL,
                        options = list(
                          `actions-box` = TRUE,
                          placeholder = "Select one or more genes"
                        ),
                        multiple = TRUE
                      ),
                      textInput(
                        'gene_heatmap_title',
                        label = "Heatmap title:",
                        value = "Gene Expression Heatmap"
                      ),
                      selectizeInput(
                        'gene_heatmap_color_by',
                        label = "Select metadata field to color by:",
                        choices = NULL,
                        options = list(`actions-box` = TRUE),
                        multiple = FALSE
                      ),
                      radioButtons("gene_heatmap_sidebar_color_palette",
                                   "Sidebar color palette:",
                                   choices = as.list(color_palette_list),
                                   selected = color_palette_list[[1]]),
                      radioButtons("gene_heatmap_color_scheme",
                                   "Select the color palette to use:",
                                   choices = as.list(heatmap_color_scheme_list),
                                   selected = heatmap_color_scheme_list[[1]]),
                      radioButtons("gene_heatmap_scaling_type",
                                   "Select the type of scaling to perform:",
                                   choices = as.list(scaling_list),
                                   selected = scaling_list[[2]]),
                      radioButtons("gene_heatmap_dendrogram_list",
                                   "Select the type of clustering to perform:",
                                   choices = as.list(dendrogram_list),
                                   selected = dendrogram_list[[2]]),
                      radioButtons("gene_heatmap_show_names",
                                   "Show axis labels:",
                                   choices = as.list(gene_heatmap_show_names_list),
                                   selected = gene_heatmap_show_names_list[[4]]),
                      radioButtons("gene_heatmap_plot_type",
                                   "Interactive or static",
                                   choices = as.list(heatmap_type_list),
                                   selected = heatmap_type_list[[1]]),
                      br(),
                      actionButton("gene_heatmap_run_button", "Create Heatmap", 
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                      br(),
                      br(),
                      downloadButton("download_gene_heatmap_data", "Download_heatmap_data.csv")
               ),
               column(9,
                      br(),
                      br(),
                      conditionalPanel(
                        "input.gene_heatmap_plot_type == 'ggplot'",
                        plotOutput("gene_expression_heatmap_ggplot", height = "700px")
                      ),
                      conditionalPanel(
                        "input.gene_heatmap_plot_type == 'heatmaply'",
                        plotlyOutput("gene_expression_heatmap_heatmaply", height = "700px")
                      )
               )
             )
    )
  ),
  tags$div(
    style = "background-color: rgba(128, 128, 128, 0.3); 
             color: #666; 
             text-align: center; 
             padding: 10px; 
             margin-top: 20px;
             font-size: 12px;
             border-radius: 5px;",
    "Version 1.0 | Developed by Shawn Fridlyand"
  )
)

#----- SERVER -----
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100 * 1024^2)
  
  # --- Upload reactives ---
  count_data <- reactive({
    req(input$counts_csv)
    formatted_counts <- read.csv(input$counts_csv$datapath, check.names = FALSE)
    colnames(formatted_counts) <- make.names(colnames(formatted_counts))
    return(formatted_counts)
  })
  
  sample_metadata <- reactive({
    req(input$metadata_csv)
    formatted_metadata <- read.csv(input$metadata_csv$datapath)
    formatted_metadata$sample_id <- make.names(formatted_metadata$sample_id)
    return(formatted_metadata)
  })
  
  # --- Sync sample removal across all tabs ---
  removed_samples <- reactiveVal(NULL)
  updating_samples <- reactiveVal(FALSE)  # Flag to prevent loops
  
  # Update all sample removal dropdowns when data loads
  observeEvent(count_data(), {
    sample_choices <- colnames(count_data())
    sample_choices <- sample_choices[!sample_choices %in% c("gene_id", "gene_name")]
    
    updating_samples(TRUE)  # Set flag before updating
    
    updateSelectizeInput(session, "pca_samples_to_remove", 
                         choices = sample_choices, 
                         selected = removed_samples(),
                         server = TRUE)
    updateSelectizeInput(session, "boxplot_samples_to_remove", 
                         choices = sample_choices, 
                         selected = removed_samples(),
                         server = TRUE)
    updateSelectizeInput(session, "sample_distance_samples_to_remove", 
                         choices = sample_choices, 
                         selected = removed_samples(),
                         server = TRUE)
    updateSelectizeInput(session, "samples_to_remove_select", 
                         choices = sample_choices, 
                         selected = removed_samples(),
                         server = TRUE)
    
    updating_samples(FALSE)  # Reset flag after updating
  })
  
  # Debounced version of each input to slow down rapid changes
  pca_samples_debounced <- reactive({
    input$pca_samples_to_remove
  }) %>% debounce(500)  # Wait 500ms after user stops selecting
  
  boxplot_samples_debounced <- reactive({
    input$boxplot_samples_to_remove
  }) %>% debounce(500)
  
  distance_samples_debounced <- reactive({
    input$sample_distance_samples_to_remove
  }) %>% debounce(500)
  
  heatmap_samples_debounced <- reactive({
    input$samples_to_remove_select
  }) %>% debounce(500)
  
  # Sync PCA removal to all other tabs
  observeEvent(pca_samples_debounced(), {
    if (updating_samples()) return()  # Skip if currently updating
    
    updating_samples(TRUE)
    removed_samples(pca_samples_debounced())
    updateSelectizeInput(session, "boxplot_samples_to_remove", selected = pca_samples_debounced())
    updateSelectizeInput(session, "sample_distance_samples_to_remove", selected = pca_samples_debounced())
    updateSelectizeInput(session, "samples_to_remove_select", selected = pca_samples_debounced())
    updating_samples(FALSE)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  # Sync Boxplot removal to all other tabs
  observeEvent(boxplot_samples_debounced(), {
    if (updating_samples()) return()
    
    updating_samples(TRUE)
    removed_samples(boxplot_samples_debounced())
    updateSelectizeInput(session, "pca_samples_to_remove", selected = boxplot_samples_debounced())
    updateSelectizeInput(session, "sample_distance_samples_to_remove", selected = boxplot_samples_debounced())
    updateSelectizeInput(session, "samples_to_remove_select", selected = boxplot_samples_debounced())
    updating_samples(FALSE)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  # Sync Sample Distance removal to all other tabs
  observeEvent(distance_samples_debounced(), {
    if (updating_samples()) return()
    
    updating_samples(TRUE)
    removed_samples(distance_samples_debounced())
    updateSelectizeInput(session, "pca_samples_to_remove", selected = distance_samples_debounced())
    updateSelectizeInput(session, "boxplot_samples_to_remove", selected = distance_samples_debounced())
    updateSelectizeInput(session, "samples_to_remove_select", selected = distance_samples_debounced())
    updating_samples(FALSE)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  # Sync Gene Expression Heatmap removal to all other tabs
  observeEvent(heatmap_samples_debounced(), {
    if (updating_samples()) return()
    
    updating_samples(TRUE)
    removed_samples(heatmap_samples_debounced())
    updateSelectizeInput(session, "pca_samples_to_remove", selected = heatmap_samples_debounced())
    updateSelectizeInput(session, "boxplot_samples_to_remove", selected = heatmap_samples_debounced())
    updateSelectizeInput(session, "sample_distance_samples_to_remove", selected = heatmap_samples_debounced())
    updating_samples(FALSE)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  # ---------------- PCA ----------------
  pca_res <- reactive({
    req(count_data(), sample_metadata())
    run_pca(count_data(), sample_metadata(), ntop = 500, remove_samples = input$pca_samples_to_remove)
  })
  
  observeEvent(pca_res(), {
    pcs <- colnames(pca_res()$coords)[grepl("^PC", colnames(pca_res()$coords))]
    if (length(pcs) >= 1) {
      updateSelectizeInput(session, "x_pc", choices = pcs, selected = pcs[1], server = TRUE)
    }
    if (length(pcs) >= 2) {
      updateSelectizeInput(session, "y_pc", choices = pcs, selected = pcs[2], server = TRUE)
    }
  })
  
  observeEvent(sample_metadata(), {
    updateSelectizeInput(session, "color_var", label = "Color samples by:",
                         choices = setdiff(colnames(sample_metadata()), pca_metadata_columns_to_remove),
                         selected = if ("condition" %in% colnames(sample_metadata())) "condition" else setdiff(colnames(sample_metadata()), pca_metadata_columns_to_remove)[1],
                         server = TRUE)
  })
  observeEvent(sample_metadata(), {
    updateSelectizeInput(session, "shape_var", label = "Change sample shapes by:",
                         choices = setdiff(colnames(sample_metadata()), pca_metadata_columns_to_remove),
                         selected = if ("condition" %in% colnames(sample_metadata())) "condition" else setdiff(colnames(sample_metadata()), pca_metadata_columns_to_remove)[1],
                         server = TRUE)
  })
  
  observeEvent(input$PCA_run_button, {
    output$pca_ggplot <- renderPlot({
      req(pca_res(), input$x_pc, input$y_pc, input$color_var)
      plot_pca(
        pca_coords = pca_res()$coords,
        variance   = pca_res()$variance,
        x_pc       = input$x_pc,
        y_pc       = input$y_pc,
        color_var  = input$color_var,
        palette    = input$pca_color_palette,
        plot_type  = "ggplot",
        shape_var  = input$shape_var
      )
    })
    
    output$pca_plotly <- renderPlotly({
      req(pca_res(), input$x_pc, input$y_pc, input$color_var)
      plot_pca(
        pca_coords = pca_res()$coords,
        variance   = pca_res()$variance,
        x_pc       = input$x_pc,
        y_pc       = input$y_pc,
        color_var  = input$color_var,
        palette    = input$pca_color_palette,
        plot_type  = "plotly",
        shape_var  = input$shape_var
      )
    })
  })
  
  output$download_PCA <- downloadHandler(
    filename = function() { paste0("PCA_coordinates_", input$color_var, ".csv") },
    content = function(file) {
      req(pca_res())
      coords <- pca_res()$coords
      meta_field <- input$color_var
      if (!is.null(meta_field) && meta_field %in% colnames(sample_metadata())) {
        meta <- sample_metadata()[, meta_field, drop = FALSE]
        # Align by sample_id
        rownames(meta) <- sample_metadata()$sample_id
        meta <- meta[match(coords$sample_id, rownames(meta)), , drop = FALSE]
        export_df <- cbind(coords, meta)
      } else {
        export_df <- coords
      }
      write.csv(export_df, file, row.names = FALSE)
    }
  )
  
  # ---------------- BOXPLOT ----------------
  observeEvent(count_data(), {
    updateSelectizeInput(session, "gene_var", label = "Select gene:", 
                         choices = if ("gene_name" %in% colnames(count_data())) unique(count_data()$gene_name) else colnames(count_data()),
                         selected = NULL, server = TRUE)
  })
  
  observeEvent(sample_metadata(), {
    updateSelectizeInput(session, "fact_var", label = "Sample groups", 
                         choices = setdiff(colnames(sample_metadata()), metadata_columns_to_remove), 
                         selected = if ("condition" %in% colnames(sample_metadata())) "condition" else setdiff(colnames(sample_metadata()), metadata_columns_to_remove)[1],
                         server = TRUE)
  })
  
  expression_df <- reactive({
    modify_df(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var, remove_samples = input$boxplot_samples_to_remove)
  })
  
  expression_table_df <- reactive({
    modify_table(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var, input$fact_var, remove_samples = input$boxplot_samples_to_remove)
  })
  
  pvalue_df <- reactive({
    table_pvalue(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var, input$fact_var, remove_samples = input$boxplot_samples_to_remove)
  })
  
  panova_df <- reactive({
    table_panova(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var, input$fact_var, remove_samples = input$boxplot_samples_to_remove)
  })
  
  observeEvent(input$run_button, {
    req(input$gene_var, input$fact_var)
    
    output$boxplot_plotly <- renderPlotly({
      req(expression_df())
      draw_boxplot(expression_df(), input$gene_var, input$fact_var, input$boxplot_color_palette, input$Log, input$Box_violin, input$Boxplot_Plot_type)
    })
    
    output$boxplot_ggplot <- renderPlot({
      req(expression_df())
      draw_boxplot(expression_df(), input$gene_var, input$fact_var, input$boxplot_color_palette, input$Log, input$Box_violin, input$Boxplot_Plot_type)
    })
    
    output$table_data <- renderDataTable(expression_table_df(), options = list(lengthMenu = c(5,30,50), pageLength = 5))
    output$table_p <- renderDataTable(pvalue_df(), options = list(lengthMenu = c(5,30,50), pageLength = 5))
    output$table_panova <- renderDataTable(panova_df(), options = list(lengthMenu = c(5,30,50), pageLength = 5))
  })
  
  output$download_exp <- downloadHandler(
    filename = function() { paste0(input$gene_var, "_expression.csv") },
    content = function(file) { write.csv(expression_table_df(), file, row.names = FALSE) }
  )
  
  output$download_pval <- downloadHandler(
    filename = function() { paste0(input$gene_var, "_pvalue.csv") },
    content = function(file) { write.csv(pvalue_df(), file, row.names = FALSE) }
  )
  
  # ---------------- SAMPLE DISTANCE HEATMAP ----------------
  observeEvent(count_data(), {
    # sample names (exclude gene_id / gene_name columns if present)
    sample_choices <- colnames(count_data())
    sample_choices <- sample_choices[!sample_choices %in% c("gene_id", "gene_name")]
    updateSelectizeInput(session, "sample_distance_samples_to_remove",
                         label = "Select samples to remove:",
                         choices = sample_choices,
                         selected = NULL,
                         server = TRUE)
  })
  
  observeEvent(sample_metadata(), {
    md_choices <- colnames(sample_metadata())
    md_choices_no_id <- setdiff(md_choices, "sample_id")
    updateSelectizeInput(session, "metadata_color_bars", choices = md_choices_no_id, selected = if ("condition" %in% md_choices_no_id) "condition" else md_choices_no_id[1], server = TRUE)
  })
  
  observeEvent(sample_metadata(), {
    updateSelectizeInput(session, "metadata_color_bars", label = "Select metadata field to color by:",
                         choices = setdiff(colnames(sample_metadata()), metadata_columns_to_remove),
                         selected = if ("condition" %in% colnames(sample_metadata())) "condition" else setdiff(colnames(sample_metadata()), metadata_columns_to_remove)[1],
                         server = TRUE)
  })
  
  # ---- Reactive: compute distance matrix on button press ----
  dist_matrix_reactive <- eventReactive(input$sample_distance_run_button, {
    req(count_data(), sample_metadata())
    
    samples_remove <- input$sample_distance_samples_to_remove
    
    # produce full dist matrix 
    run_sample_distance(
      counts = count_data(),
      metadata = sample_metadata(),
      remove_samples = samples_remove,
      scale_type = input$sample_distance_heatmap_scaling_type,
      ntop = 500
    )
  })
  
  # ---- Render static ggplot heatmap ----
  output$sample_distance_heatmap <- renderPlot({
    req(count_data(), sample_metadata())
    
    p <- plot_sample_distance_heatmap(
      dist_matrix = dist_matrix_reactive(),
      metadata = sample_metadata(),
      color_by = input$metadata_color_bars,
      sidebar_color_scheme = "Set1",
      heatmap_title = "Sample Distance Heatmap",
      heatmap_color_scheme = input$sample_distance_heatmap_color_palette,
      column_text_angle = 45,
      show_dendrogram = c(
        input$sample_distance_heatmap_dendrogram_list %in% c("row", "both"),
        input$sample_distance_heatmap_dendrogram_list %in% c("column", "both")
      ),
      plot_type = "ggplot",
      show_tick_labels = c(
        input$sample_distance_heatmap_show_names %in% c("row", "both"),
        input$sample_distance_heatmap_show_names %in% c("column", "both")
      )
    )
    print(p)
  })
  
  # ---- Render heatmaply ----
  output$sample_distance_heatmaply <- renderPlotly({
    req(count_data(), sample_metadata())
    
    hm <- plot_sample_distance_heatmap(
      dist_matrix = dist_matrix_reactive(),
      metadata = sample_metadata(),
      color_by = input$metadata_color_bars,
      sidebar_color_scheme = "Set1",
      heatmap_title = "Sample Distance Heatmap",
      heatmap_color_scheme = input$sample_distance_heatmap_color_palette,
      column_text_angle = 45,
      show_dendrogram = c(
        input$sample_distance_heatmap_dendrogram_list %in% c("row", "both"),
        input$sample_distance_heatmap_dendrogram_list %in% c("column", "both")
      ),
      plot_type = "heatmaply",
      show_tick_labels = c(
        input$sample_distance_heatmap_show_names %in% c("row", "both"),
        input$sample_distance_heatmap_show_names %in% c("column", "both")
      ),
      colorbar_xpos = 1.02,
      colorbar_ypos = 0.5,
      colorbar_len = 0.4
    )
    return(hm)
  })
  
  # ---- Download distance matrix ----
  output$download_distance_matrix <- downloadHandler(
    filename = function() { paste0("sample_distance_matrix_", Sys.Date(), ".csv") },
    content = function(file) {
      req(dist_matrix_reactive())
      write.csv(dist_matrix_reactive(), file, row.names = TRUE)
    }
  )
  
  # ------------------ GENE EXPRESSION HEATMAP ------------------
  
  # Update gene list choices when count data loads
  observeEvent(count_data(), {
    gene_choices <- if ("gene_name" %in% colnames(count_data())) {
      unique(count_data()$gene_name)
    } else {
      rownames(count_data())
    }
    updateSelectizeInput(session, "gene_list_select", 
                         label = "Select genes:",
                         choices = gene_choices,
                         selected = NULL,
                         server = TRUE)
    
    # Update sample removal choices (exclude gene_id/gene_name columns)
    sample_choices <- colnames(count_data())
    sample_choices <- sample_choices[!sample_choices %in% c("gene_id", "gene_name")]
    updateSelectizeInput(session, "samples_to_remove_select",
                         label = "Select samples to remove:",
                         choices = sample_choices,
                         selected = NULL,
                         server = TRUE)
  })
  
  # Update color_by choices when metadata loads
  observeEvent(sample_metadata(), {
    updateSelectizeInput(session, "gene_heatmap_color_by", 
                         label = "Color samples by:",
                         choices = setdiff(colnames(sample_metadata()), metadata_columns_to_remove),
                         selected = if ("condition" %in% colnames(sample_metadata())) "condition" else setdiff(colnames(sample_metadata()), metadata_columns_to_remove)[1],
                         server = TRUE)
  })
  
  # NEW: Add gene annotations reactive (if you have gene annotation data)
  gene_annotations_reactive <- reactive({
    # If you have a gene annotations file input
    if (!is.null(input$gene_annotations_file)) {
      read.csv(input$gene_annotations_file$datapath)
    } else if (exists("gene_annotations")) {
      # Or if you have it pre-loaded as a global variable
      gene_annotations
    } else {
      # Return NULL if no annotations available
      NULL
    }
  })
  
  # Reactive: prepare gene list and sample list
  gene_sample_lists_reactive <- eventReactive(input$gene_heatmap_run_button, {
    req(count_data(), input$gene_list_select)
    
    # Get selected genes
    gene_list <- input$gene_list_select
    
    # Get sample list (all samples minus removed ones)
    all_samples <- colnames(count_data())
    all_samples <- all_samples[!all_samples %in% c("gene_id", "gene_name")]
    
    samples_remove <- input$samples_to_remove_select
    sample_list <- if (!is.null(samples_remove) && length(samples_remove) > 0) {
      setdiff(all_samples, samples_remove)
    } else {
      all_samples
    }
    
    list(
      gene_list = gene_list,
      sample_list = sample_list
    )
  })
  
  # Render static ggplot heatmap
  output$gene_expression_heatmap_ggplot <- renderPlot({
    req(gene_sample_lists_reactive())
    
    lists <- gene_sample_lists_reactive()
    
    p <- plot_gene_expression_heatmap(
      counts = count_data(),
      gene_list = lists$gene_list,
      gene_annotations = gene_annotations_reactive(),
      sample_list = lists$sample_list,
      metadata = sample_metadata(),
      heatmap_title = input$gene_heatmap_title,
      color_by = input$gene_heatmap_color_by,
      sidebar_color_scheme = input$gene_heatmap_sidebar_color_palette,
      heatmap_color_scheme = input$gene_heatmap_color_scheme,
      xlab = "",
      ylab = "Genes",
      column_text_angle = 90,
      legend_title = "Expression",
      cex_row = 0.5,
      cex_col = 0.5,
      scaling = input$gene_heatmap_scaling_type,
      cluster = input$gene_heatmap_clustering_type,
      dendrogram = input$gene_heatmap_dendrogram_list,
      show_names = input$gene_heatmap_show_names,
      heatmap_type = "ggheatmap"
    )
    print(p)
  })
  
  # Render interactive heatmaply
  output$gene_expression_heatmap_heatmaply <- renderPlotly({
    req(gene_sample_lists_reactive())
    
    lists <- gene_sample_lists_reactive()
    
    plot_gene_expression_heatmap(
      counts = count_data(),
      gene_list = lists$gene_list,
      gene_annotations = gene_annotations_reactive(),
      sample_list = lists$sample_list,
      metadata = sample_metadata(),
      heatmap_title = input$gene_heatmap_title,
      color_by = input$gene_heatmap_color_by,
      sidebar_color_scheme = input$gene_heatmap_sidebar_color_palette,
      heatmap_color_scheme = input$gene_heatmap_color_scheme,
      xlab = "",
      ylab = "Genes",
      column_text_angle = 90,
      legend_title = "Expression",
      cex_row = 0.5,
      cex_col = 0.5,
      scaling = input$gene_heatmap_scaling_type,
      cluster = input$gene_heatmap_clustering_type,
      dendrogram = input$gene_heatmap_dendrogram_list,
      show_names = input$gene_heatmap_show_names,
      heatmap_type = "heatmaply"
    )
  })
  
  # Download gene expression matrix (updated to use new function)
  output$download_gene_heatmap_data <- downloadHandler(
    filename = function() { 
      paste0("gene_expression_heatmap_", input$gene_heatmap_scaling_type, "_", Sys.Date(), ".csv") 
    },
    content = function(file) {
      req(gene_sample_lists_reactive())
      
      lists <- gene_sample_lists_reactive()
      
      # Prepare the matrix
      expr_matrix <- prepare_gene_expression_matrix(
        counts = count_data(),
        metadata = sample_metadata(),
        gene_list = lists$gene_list,
        remove_samples = input$samples_to_remove_select
      )
      
      write.csv(expr_matrix, file, row.names = TRUE)
    }
  )
  
}
shinyApp(ui = ui, server = server)
