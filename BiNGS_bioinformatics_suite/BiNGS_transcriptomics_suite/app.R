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
             plotOutput("pca_plot")  # placeholder for PCA plot
    ),
    
    # ------------------ BOXPLOT TAB ------------------
    tabPanel("Boxplot",
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
             radioButtons("color_palette",
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
             radioButtons("Plot_type",
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
             br(), br(), br(),
             downloadButton("download_exp", "Download_expression.csv"),
             br(), br(),
             downloadButton("download_pval", "Download_pvalue_table.csv"),
             
             # ------------------ BOXPLOT OUTPUTS ------------------
             conditionalPanel(condition = "input.Plot_type == 'plotly'",
                              plotlyOutput("boxplot_plotly")),
             conditionalPanel(condition = "input.Plot_type == 'ggplot'",
                              plotOutput("boxplot_ggplot")),
             h3("P-value Table"),
             DT::dataTableOutput("table_p"),
             h3("ANOVA Table"),
             DT::dataTableOutput("table_panova"),
             h3("Gene Expression Table"),
             DT::dataTableOutput("table_data")
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

# Define server logic required to draw a histogram
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
    
    output$boxplot_plotly = renderPlotly(
      draw_boxplot(expression_df(), 
                   input$gene_var, 
                   input$fact_var, 
                   input$color_palette,
                   input$Log, 
                   input$Box_violin, 
                   input$Plot_type
      )
    )
    output$boxplot_ggplot = renderPlot(
      draw_boxplot(expression_df(), 
                   input$gene_var, 
                   input$fact_var, 
                   input$color_palette,
                   input$Log, 
                   input$Box_violin, 
                   input$Plot_type
      )
    )
    output$table_data = renderDataTable(
      expression_table_df(), 
      options = list(lengthMenu = c(5, 30, 50), 
                     pageLength = 5)
    )
    output$table_p = renderDataTable(
      pvalue_df(), 
      options = list(lengthMenu = c(5, 30, 50), 
                     pageLength = 5)
    )
    output$table_panova = renderDataTable(
      panova_df(), 
      options = list(lengthMenu = c(5, 30, 50), 
                     pageLength = 5)
    )
  })
  
  output$download_exp <- downloadHandler(filename = function() {paste0(input$gene_var, "_expression.csv")}, 
                                         content = function(file) {write.csv(expression_table_df(), file, row.names = FALSE)})
  
  output$download_pval <- downloadHandler(filename = function() {paste0(input$gene_var, "_pvalue.csv")}, 
                                          content = function(file) { write.csv(pvalue_df(), file, row.names = FALSE)})
  
}

shinyApp(ui = ui, server = server)
