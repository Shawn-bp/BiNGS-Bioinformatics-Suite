require(shiny)
require(shinythemes)
require(data.table)
require(ggplot2)
require(dplyr)
require(DT)
require(plotly)
require(fst)
require(ggpubr)
require(heatmaply)
library(RColorBrewer)
source("helper.R")
source("global.R")

ui <- fluidPage(
  titlePanel("BiNGS Transcriptomics Suite"),
  
  tabsetPanel(
    
    # ------------------ UPLOAD TAB ------------------
    tabPanel("Upload",
             fileInput("counts_csv", "Select normalized counts file to import", accept = ".csv"),
             fileInput("metadata_csv", "Select metadata file to import", accept = ".csv")
    ),
    
    # ------------------ PCA TAB ------------------
    tabPanel("PCA",
             fluidRow(
               column(3,
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
             radioButtons("pca_color_palette",
                          "Select the color palette to use:",
                          choices = as.list(color_palette_list),
                          selected = color_palette_list[[1]]),
             radioButtons("PCA_Plot_type",
                          "Plotly or GGplot:",
                          choices = as.list(type_list),
                          selected = type_list[[1]]),
             actionButton("PCA_run_button", "Create PCA", 
                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
             br(), 
             br(), 
             br(),
             downloadButton("download_PCA", "Download_PCA_coords.csv"),
             br(), 
             br(),
             downloadButton("download_variance", "Download_variance_table.csv")
               ),
    # ------------------ PCA outputs ------------------
             column(9,
             conditionalPanel("input.PCA_Plot_type == 'plotly'",
                                plotlyOutput("pca_plotly", height = "600px")),
             conditionalPanel("input.PCA_Plot_type == 'ggplot'",
                                plotOutput("pca_ggplot", height = "600px")),
               
               h3("Variance Explained"),
               DT::dataTableOutput("pca_variance")
             ),
             )
    ),
    # ------------------ BOXPLOT TAB ------------------
    tabPanel("Boxplot",
             fluidRow(
               column(3,
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
                          "Plotly or GGplot:",
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
  # ------------------ HEATMAP TABS ------------------
  tabPanel("Sample Distance Heatmap",
           fluidRow(
             column(3,
                    selectizeInput(
                      'samples_x_var',
                      label = "Select samples on X-axis:",
                      choices = NULL,
                      options = list(`actions-box` = TRUE),
                      multiple = TRUE
                    ),
                    selectizeInput(
                      'samples_y_var',
                      label = "Select samples on Y-axis:",
                      choices = NULL,
                      options = list(`actions-box` = TRUE),
                      multiple = TRUE
                    ),
                    selectizeInput(
                      'metadata_color_bars',
                      label = "Select metadata field to color by:",
                      choices = NULL,
                      options = list(`actions-box` = TRUE),
                      multiple = FALSE
                    ),
                    radioButtons("sample_distance_heatmap_color_palette",
                                 "Select the color palette to use:",
                                 choices = as.list(color_palette_list),
                                 selected = color_palette_list[[1]]),
                    radioButtons("sample_distance_heatmap_scaling_type",
                                 "Select the type of scaling to perform:",
                                 choices = as.list(scaling_list),
                                 selected = scaling_list[[1]]),
                    radioButtons("sample_distance_heatmap_clustering_type",
                                 "Select the type of clustering to perform:",
                                 choices = as.list(clustering_list),
                                 selected = clustering_list[[1]]),
                    radioButtons("sample_distance_heatmap_dendrogram_list",
                                 "Select the dendrograms to show:",
                                 choices = as.list(dendrogram_list),
                                 selected = dendrogram_list[[1]]),
                    radioButtons("Sample_Distance_Heatmap_Plot_type",
                                 "Heatmaply or GGplot:",
                                 choices = as.list(heatmap_type_list),
                                 selected = heatmap_type_list[[1]]),
                    br(),
                    actionButton("sample_distance_run_button", "Create Heatmap", 
                                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
             ),
             # ------------------ HEATMAP OUTPUTS ------------------
             column(9,
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
    tabPanel("Gene Expression Heatmap",
             plotOutput("gene_expression_heatmap_plot")
    )
  )
)
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
  
  # ---------------- PCA ----------------
  pca_res <- reactive({
    req(count_data(), sample_metadata())
    run_pca(count_data(), sample_metadata(), scale_data = TRUE)
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
                         choices = setdiff(colnames(sample_metadata()), metadata_columns_to_remove),
                         selected = if ("condition" %in% colnames(sample_metadata())) "condition" else setdiff(colnames(sample_metadata()), metadata_columns_to_remove)[1],
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
        plot_type  = "ggplot"
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
        plot_type  = "plotly"
      )
    })
    
    output$pca_variance <- DT::renderDataTable({
      req(pca_res())
      data.frame(
        PC = paste0("PC", seq_along(pca_res()$variance)),
        Variance = round(pca_res()$variance, 2)
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
  
  output$download_variance <- downloadHandler(
    filename = function() { "PCA_variance_table.csv" },
    content = function(file) {
      req(pca_res())
      df <- data.frame(PC = paste0("PC", seq_along(pca_res()$variance)), Variance = pca_res()$variance)
      write.csv(df, file, row.names = FALSE)
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
    modify_df(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var)
  })
  
  expression_table_df <- reactive({
    modify_table(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var, input$fact_var)
  })
  
  pvalue_df <- reactive({
    table_pvalue(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var, input$fact_var)
  })
  
  panova_df <- reactive({
    table_panova(count_data(), sample_metadata(), input$Log, input$QC_check, input$gene_var, input$fact_var)
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
  
  # ---------------- HEATMAP ----------------
  observeEvent(count_data(), {
    # sample names (exclude gene_id / gene_name columns if present)
    sample_choices <- colnames(count_data())
    sample_choices <- sample_choices[!sample_choices %in% c("gene_id", "gene_name")]
    updateSelectizeInput(session, "samples_x_var", choices = sample_choices, selected = sample_choices, server = TRUE)
    updateSelectizeInput(session, "samples_y_var", choices = sample_choices, selected = sample_choices, server = TRUE)
  })
  
  observeEvent(sample_metadata(), {
    md_choices <- colnames(sample_metadata())
    # exclude sample_id from color choices (but you can include it if desired)
    md_choices_no_id <- setdiff(md_choices, "sample_id")
    updateSelectizeInput(session, "metadata_color_bars", choices = md_choices_no_id, selected = if ("condition" %in% md_choices_no_id) "condition" else md_choices_no_id[1], server = TRUE)
  })
  
  # ---- Reactive: compute distance matrix on button press ----
  dist_matrix_reactive <- eventReactive(input$sample_distance_run_button, {
    req(count_data(), sample_metadata())
    
    # produce full dist matrix (will be subset for plotting below if needed)
    run_sample_distance(
      counts = count_data(),
      metadata = sample_metadata(),
      remove_samples = NULL,
      scale_type = input$sample_distance_heatmap_scaling_type
    )
  })
  
  # ---- Render static ggplot heatmap ----
  output$sample_distance_heatmap <- renderPlot({
    req(dist_matrix_reactive())
    dist_mat <- dist_matrix_reactive()
    
    # Optionally subset samples_x_var / samples_y_var (if user selected specific)
    sel_x <- input$samples_x_var
    sel_y <- input$samples_y_var
    # If UI select inputs are single (non-multiple) you provided, allow NULL or single
    if (!is.null(sel_x) && length(sel_x) > 0) {
      keep_cols_x <- intersect(colnames(dist_mat), sel_x)
    } else {
      keep_cols_x <- colnames(dist_mat)
    }
    if (!is.null(sel_y) && length(sel_y) > 0) {
      keep_cols_y <- intersect(rownames(dist_mat), sel_y)
    } else {
      keep_cols_y <- rownames(dist_mat)
    }
    # subset and re-order to selected axes intersection
    sub_mat <- dist_mat[keep_cols_y, keep_cols_x, drop = FALSE]
    
    # Create static ggplot heatmap via your helper
    p <- plot_sample_distance_heatmap(
      dist_matrix = sub_mat,
      metadata = sample_metadata(),
      color_scheme = input$sample_distance_heatmap_color_palette,
      cluster = input$sample_distance_heatmap_clustering_type,
      dendrogram = input$sample_distance_heatmap_dendrogram_list,
      show_names = "both",
      color_by = input$metadata_color_bars,
      heatmap_type = "ggplot"
    )
    print(p)
  })
  
  # ---- Render interactive heatmaply ----
  output$sample_distance_heatmaply <- renderPlotly({
    req(dist_matrix_reactive())
    dist_mat <- dist_matrix_reactive()
    
    sel_x <- input$samples_x_var
    sel_y <- input$samples_y_var
    if (!is.null(sel_x) && length(sel_x) > 0) {
      keep_cols_x <- intersect(colnames(dist_mat), sel_x)
    } else {
      keep_cols_x <- colnames(dist_mat)
    }
    if (!is.null(sel_y) && length(sel_y) > 0) {
      keep_cols_y <- intersect(rownames(dist_mat), sel_y)
    } else {
      keep_cols_y <- rownames(dist_mat)
    }
    sub_mat <- dist_mat[keep_cols_y, keep_cols_x, drop = FALSE]
    
    hm <- plot_sample_distance_heatmap(
      dist_matrix = sub_mat,
      metadata = sample_metadata(),
      color_scheme = input$sample_distance_heatmap_color_palette,
      cluster = input$sample_distance_heatmap_clustering_type,
      dendrogram = input$sample_distance_heatmap_dendrogram_list,
      show_names = "both",
      color_by = input$metadata_color_bars,
      heatmap_type = "heatmaply"
    )
    # heatmaply returns an htmlwidget (plotly); return it directly
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
  
}


shinyApp(ui = ui, server = server)
