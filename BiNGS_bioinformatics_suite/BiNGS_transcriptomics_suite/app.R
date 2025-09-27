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
  titlePanel("Transcriptomics Suite"),
  
  sidebarLayout(
    sidebarPanel(
      title = "Upload",
      fileInput("counts_csv", "Select normalized counts file to import", accept = ".csv"),
      fileInput("metadata_csv", "Select metadata file to import", accept = ".csv"),
      selectizeInput(
        'gene_var',
        label = "Select Gene:",
        choices = NULL,
        options = list(`actions-box` = TRUE),
        multiple = FALSE
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("PCA", plotlyOutput("pca_plot")),
        tabPanel("Boxplot", selectizeInput(
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
                     choices = as.list(plot_list),
                     selected = plot_list[[1]]),
        br(),
        actionButton("run_button", "Create Boxplot", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
        br(),
        br(),
        br(),
        downloadButton("download_exp", "Download_expression.csv"),
        br(),
        br(),
        downloadButton("download_pval", "Download_pvalue_table.csv")
        ),
        # Show a plot of the generated distribution
        mainPanel(
          
          #conditional statement for which type of plot based on input
          #https://shiny.posit.co/r/reference/shiny/1.4.0/conditionalpanel
          conditionalPanel(condition = "input.Plot_type == 'plotly'",
                           plotlyOutput("boxplot_plotly")),
          conditionalPanel(condition = "input.Plot_type == 'ggplot'",
                           plotOutput("boxplot_ggplot")),
          h3("P-value Table"),
          DT::dataTableOutput("table_p"),
          h3("ANOVA Table"),
          DT::dataTableOutput("table_panova"),
          h3("Gene Expression Table"),
          DT::dataTableOutput("table_data"),
          
        ),
        tabPanel("Sample Distance Heatmap", plotOutput("sample_distance_heatmap_plot")),
        tabPanel("Gene Expression Heatmap", plotOutput("gene_expression_heatmap_plot"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive: load counts
  counts_data <- reactive({
    req(input$counts_csv)
    read.csv(input$counts_csv$datapath, row.names = 1)
  })
  
  # Reactive: load metadata
  metadata <- reactive({
    req(input$metadata_csv)
    read.csv(input$metadata_csv$datapath)
  })
  
  # Update gene dropdown when counts are uploaded
  observeEvent(counts_data(), {
    updateSelectizeInput(
      session,
      "gene_var",
      choices = counts_data()$gene_name,   # uses the gene_name column
      server = TRUE
    )
  })
}


shinyApp(ui = ui, server = server)