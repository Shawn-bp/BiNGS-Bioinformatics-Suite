require(shiny)
require(shinythemes)
require(data.table)
require(ggplot2)
require(dplyr)
require(DT)
require(plotly)
require(fst)
require(ggpubr)
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
             plotOutput("sample_distance_heatmap_plot")
    ),
    
    tabPanel("Gene Expression Heatmap",
             plotOutput("gene_expression_heatmap_plot")
    )
  )
)
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=100*1024^2) 
  
  count_data <- reactive({
    req(input$counts_csv)
    formatted_counts = read.csv(input$counts_csv$datapath, check.names = FALSE)
    formatted_colnames = make.names(colnames(formatted_counts))
    colnames(formatted_counts) = formatted_colnames
    return(formatted_counts)
  })
  
  sample_metadata <- reactive({
    req(input$metadata_csv)
    formatted_metadata = read.csv(input$metadata_csv$datapath)
    formatted_metadata$sample_id = make.names(formatted_metadata$sample_id)
    return(formatted_metadata)
  })
  
  # Reactive PCA
  pca_res <- reactive({
    req(count_data(), sample_metadata())
    run_pca(count_data(), sample_metadata(),
            scale_data = TRUE)
  })
  
  # Populate PC choices
  observeEvent(pca_res(), {
    pcs <- colnames(pca_res()$coords)[grepl("^PC", colnames(pca_res()$coords))]
    updateSelectizeInput(session, "x_pc", choices = pcs, selected = pcs[1])
    updateSelectizeInput(session, "y_pc", choices = pcs, selected = pcs[2])
  })
  
  # Populate color variable
  observeEvent(sample_metadata(), {
    updateSelectizeInput(session, "color_var", label = "Color samples by:",
                         choices = setdiff(colnames(sample_metadata()), metadata_columns_to_remove),
                         selected = "condition",
                         server = TRUE
    )
  })
  
  # Render plot
  observeEvent(input$PCA_run_button, {
    
    output$pca_ggplot = renderPlot({
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
    
    output$pca_plotly = renderPlotly({
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
    
    # Variance explained table
    output$pca_variance <- DT::renderDataTable({
      req(pca_res())
      data.frame(
        PC = paste0("PC", seq_along(pca_res()$variance)),
        Variance = round(pca_res()$variance, 2)
      )
    })
    
  })
  
  output$download_PCA <- downloadHandler(
    filename = function() {
      paste0("PCA_coordinates_", input$color_var, ".csv")
    },
    content = function(file) {
      # Extract PCA coordinates
      coords <- pca_res()$coords   # assumes samples as rows, PCs as columns
      
      # Get selected metadata column
      meta_field <- input$color_var
      meta <- sample_metadata()[, meta_field, drop = FALSE]  # returns a one-column data.frame
      
      # Make sure rows align by sample
      meta <- meta[rownames(coords), , drop = FALSE]
      
      # Combine coords + metadata
      export_df <- cbind(coords, meta)
      
      # Write to CSV
      write.csv(export_df, file, row.names = TRUE)
    }
  )
  
  output$download_variance <- downloadHandler(
    filename = function() {
      "PCA_variance_table.csv"
    },
    content = function(file) {
      # grab your variance table
      variance_df <- pca_res()$variance
      
      # write it to CSV
      write.csv(variance_df, file, row.names = TRUE)
    }
  )
  
  
  observeEvent(count_data(), {012
    updateSelectizeInput(session, "gene_var", label = "Select gene:", 
                         choices = unique(count_data()$gene_name), 
                         selected = NULL, 
                         server = TRUE)
  })
  
  observeEvent(sample_metadata(), {
    updateSelectizeInput(session, "fact_var", label = "Sample groups", 
                         choices = setdiff(colnames(sample_metadata()), metadata_columns_to_remove), 
                         selected = "condition", 
                         server = TRUE)
  })
  
  expression_df = reactive({
    modify_df(count_data(), 
              sample_metadata(), 
              input$Log, 
              input$QC_check, 
              input$gene_var)
  })
  
  expression_table_df = reactive({
    modify_table(count_data(), 
                 sample_metadata(), 
                 input$Log, 
                 input$QC_check, 
                 input$gene_var, 
                 input$fact_var)
  })
  
  pvalue_df = reactive({
    table_pvalue(count_data(), 
                 sample_metadata(), 
                 input$Log, 
                 input$QC_check, 
                 input$gene_var, 
                 input$fact_var)
  })
  
  panova_df = reactive({
    table_panova(count_data(), 
                 sample_metadata(), 
                 input$Log, 
                 input$QC_check, 
                 input$gene_var, 
                 input$fact_var)
  })
  
  observeEvent(input$run_button, {
    req(input$gene_var, input$fact_var)  # <--- ensures inputs exist
    
    output$boxplot_plotly = renderPlotly({
      req(expression_df())
      draw_boxplot(
        expression_df(), 
        input$gene_var, 
        input$fact_var, 
        input$boxplot_color_palette, 
        input$Log, 
        input$Box_violin, 
        input$Boxplot_Plot_type
      )
    })
    
    output$boxplot_ggplot = renderPlot({
      req(expression_df())
      draw_boxplot(
        expression_df(), 
        input$gene_var, 
        input$fact_var, 
        input$boxplot_color_palette, 
        input$Log, 
        input$Box_violin, 
        input$Boxplot_Plot_type
      )
    })
    
    output$table_data = renderDataTable(
      expression_table_df(), 
      options = list(lengthMenu = c(5, 30, 50), pageLength = 5)
    )
    output$table_p = renderDataTable(
      pvalue_df(), 
      options = list(lengthMenu = c(5, 30, 50), pageLength = 5)
    )
    output$table_panova = renderDataTable(
      panova_df(), 
      options = list(lengthMenu = c(5, 30, 50), pageLength = 5)
    )
  })
  
  
  output$download_exp <- downloadHandler(filename = function() {paste0(input$gene_var, "_expression.csv")}, 
                                         content = function(file) {write.csv(expression_table_df(), file, row.names = FALSE)})
  
  output$download_pval <- downloadHandler(filename = function() {paste0(input$gene_var, "_pvalue.csv")}, 
                                          content = function(file) { write.csv(pvalue_df(), file, row.names = FALSE)})
  
}

shinyApp(ui = ui, server = server)
