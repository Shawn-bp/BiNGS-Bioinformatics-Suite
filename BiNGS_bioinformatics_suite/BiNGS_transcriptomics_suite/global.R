# ------------------ PCA ------------------
run_pca <- function(counts, metadata, ntop = 500, remove_samples = NULL) {
  
  gene_cols <- intersect(colnames(counts), c("gene_id", "gene_name"))
  
  sample_cols <- setdiff(colnames(counts), gene_cols)
  
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    sample_cols <- setdiff(sample_cols, remove_samples)
    metadata <- metadata[metadata$sample_id %in% sample_cols, ]
  }
  
  expr_matrix <- as.matrix(counts[, sample_cols, drop = FALSE])
  
  if ("gene_id" %in% colnames(counts)) {
    rownames(expr_matrix) <- counts$gene_id
  }
  
  rv <- matrixStats::rowVars(expr_matrix)
  
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca <- prcomp(t(expr_matrix[select, ]))
  
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  
  pca_coords <- as.data.frame(pca$x)
  pca_coords$sample_id <- rownames(pca_coords)
  
  if (!is.null(metadata)) {
    pca_coords <- dplyr::left_join(pca_coords, metadata, by = "sample_id")
  }
  
  return(list(
    pca = pca,
    coords = pca_coords,
    variance = percentVar * 100
  ))
}

plot_pca <- function(pca_coords, variance, x_pc, y_pc, color_var = NULL, shape_var = NULL,
                     palette = "Set1", plot_type = "ggplot") {
  
  x_label <- paste0(x_pc, " (", round(variance[as.numeric(gsub("PC", "", x_pc))], 2), "%)")
  y_label <- paste0(y_pc, " (", round(variance[as.numeric(gsub("PC", "", y_pc))], 2), "%)")
  
  if (plot_type == "ggplot") {
    p <- ggplot(pca_coords, aes_string(x = x_pc, y = y_pc, color = color_var, shape = shape_var)) +
      geom_point(size = 6, alpha = 0.8) +
      scale_color_brewer(palette = palette, name = color_var) +
      labs(
        title = "Principal Component Analysis",
        x = x_label, 
        y = y_label
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.key = element_rect(fill = "white"),
      )
    
    if (!is.null(shape_var) && shape_var != "") {
      p <- p + scale_shape_discrete(name = shape_var)
    }
    
    return(p)
    
  } else if (plot_type == "plotly") {
    
    symbol_var <- NULL
    if (!is.null(shape_var) && shape_var != "") {
      symbol_var <- pca_coords[[shape_var]]
    }
    
    p <- plotly::plot_ly(
      pca_coords,
      x = ~get(x_pc),
      y = ~get(y_pc),
      color = if (!is.null(color_var)) pca_coords[[color_var]] else NULL,
      symbol = symbol_var,
      symbols = c('circle', 'square', 'diamond', 'cross', 'triangle-up', 'triangle-down'),
      colors = RColorBrewer::brewer.pal(max(3, length(unique(pca_coords[[color_var]]))), palette),
      type = "scatter",
      mode = "markers",
      marker = list(size = 15),
      text = ~sample_id,
      hovertemplate = paste0(
        "<b>%{text}</b><br>",
        x_pc, ": %{x:.2f}<br>",
        y_pc, ": %{y:.2f}<br>",
        if (!is.null(shape_var) && shape_var != "") paste0(shape_var, ": ", pca_coords[[shape_var]], "<br>") else "",
        "<extra></extra>"
      )
    ) %>% plotly::layout(
      title = list(text = "Principal Component Analysis", x = 0.5, xanchor = "center"),
      xaxis = list(
        title = x_label,
        zerolinecolor = '#969696',
        zerolinewidth = 2,
        linecolor = '#636363',
        linewidth = 2,
        gridcolor = '#e0e0e0'
      ),
      yaxis = list(
        title = y_label,
        zerolinecolor = '#969696',
        zerolinewidth = 2,
        linecolor = '#636363',
        linewidth = 2,
        gridcolor = '#e0e0e0'
      ),
      legend = list(
        title = list(text = paste0('<b>', if (!is.null(color_var)) color_var else "", '</b>')),
        x = 1.02,
        y = 1,
        xanchor = "left",
        bgcolor = 'rgba(255,255,255,0.9)',
        bordercolor = '#636363',
        borderwidth = 1
      )
    ) %>% configure_plotly_panel(exportFormat = "svg")
    return(p)
  }
}
# ------------------ BOXPLOT ------------------
get_factor_comparisons = function(condition = c(), 
                                  reference_group = NULL) {
  condition = sort(unique(as.character(condition)))
  if (length(condition) <= 1) {
    comparisons = data.frame(comparison_name = character(), 
                             comparison_group = character(), 
                             reference_group = character())
    return(comparisons)    
  } 
  
  if (is.null(reference_group) == TRUE) {
    comparisons = t(utils::combn(x = condition, m = 2, simplify = TRUE))
  } else {
    comparisons = t(sapply(condition, function(cond) {c(cond, reference_group)}))
  }
  comparisons = as.data.frame(comparisons, stringsAsFactors = FALSE)
  colnames(comparisons) = c("comparison_group", "reference_group")
  
  comparisons = comparisons[which(comparisons$comparison_group != comparisons$reference_group), ]
  comparisons = comparisons[order(comparisons$reference_group), , drop = FALSE]
  
  comparisons$comparison_name = paste0(comparisons$comparison_group, "_vs_", comparisons$reference_group)
  comparisons = comparisons[, c("comparison_name", "comparison_group", "reference_group"), drop = FALSE]
  rownames(comparisons) = NULL
  
  return(comparisons)
}

modify_df <- function(counts_df, metadata, log, QC_check, gene, remove_samples = NULL){
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    keep_samples <- setdiff(metadata$sample_id, remove_samples)
    metadata <- metadata[metadata$sample_id %in% keep_samples, ]
    keep_cols <- c(intersect(colnames(counts_df), c("gene_id", "gene_name")), keep_samples)
    counts_df <- counts_df[, colnames(counts_df) %in% keep_cols, drop = FALSE]
  }
  
  if(QC_check == "no"){
    if(log == "yes"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
      combination_df[gene] = log2(combination_df[gene]+1)
    }
    if(log == "no"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
    }}
  if(QC_check == "yes"){
    if(log == "yes"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
      combination_df[gene] = log2(combination_df[gene]+1)
      combination_df = filter(combination_df, quality_check == "pass")
    } 
    if(log == "no"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
      combination_df = filter(combination_df, quality_check == "pass")
    }}
  return(combination_df)
}


modify_table <- function(counts_df, metadata, log, QC_check, gene, fact_var, remove_samples = NULL){
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    keep_samples <- setdiff(metadata$sample_id, remove_samples)
    metadata <- metadata[metadata$sample_id %in% keep_samples, ]
    keep_cols <- c(intersect(colnames(counts_df), c("gene_id", "gene_name")), keep_samples)
    counts_df <- counts_df[, colnames(counts_df) %in% keep_cols, drop = FALSE]
  }
  if(QC_check == "no"){
    if(log == "yes"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
      combination_df = combination_df[,unique(c("sample_id", "condition", "replicate","quality_check", fact_var, gene))]
      combination_df
      combination_df[gene] = log2(combination_df[gene]+1)
      rownames(combination_df) = NULL
    }
    if(log == "no"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
      combination_df = combination_df[,unique(c("sample_id", "condition", "replicate","quality_check", fact_var, gene))]
      rownames(combination_df) = NULL
    }}
  if(QC_check == "yes"){
    if(log == "yes"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
      combination_df = combination_df[,unique(c("sample_id", "condition", "replicate","quality_check", fact_var, gene))]
      combination_df[gene] = log2(combination_df[gene]+1)
      combination_df = filter(combination_df, quality_check == "pass")
      rownames(combination_df) = NULL
    } 
    if(log == "no"){
      counts_table_logged = counts_df[,metadata$sample_id]
      transposed_counts = as.data.frame(t(counts_table_logged))
      colnames(transposed_counts) = counts_df$gene_name
      combination_df = cbind(metadata, transposed_counts)
      combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
      combination_df = combination_df[,unique(c("sample_id", "condition", "replicate","quality_check", fact_var, gene))]
      combination_df = filter(combination_df, quality_check == "pass")
      rownames(combination_df) = NULL
    }}
  return(combination_df)
}

table_pvalue = function(counts_df, metadata, log, QC_check, gene, fact_var, remove_samples = NULL){
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    keep_samples <- setdiff(metadata$sample_id, remove_samples)
    metadata <- metadata[metadata$sample_id %in% keep_samples, ]
    keep_cols <- c(intersect(colnames(counts_df), c("gene_id", "gene_name")), keep_samples)
    counts_df <- counts_df[, colnames(counts_df) %in% keep_cols, drop = FALSE]
  }
  dups = table(metadata[,fact_var])
  dups_df = as.data.frame(dups)
  save_s = c()
  if(sum(dups_df[,2]) == length(unique(metadata[,fact_var]))){
    ptab = data.frame(Error = c("P-value requires more than 1 replicate per group!"))
  } else if(length(unique(metadata[,fact_var]))<=1){
    ptab = data.frame(Error = c("P-value requires more than 1 grouping factor!"))
  }else{
    for(i in 1:nrow(dups_df)){
      if(dups_df[i,2] == 1){
        g = dups_df[i,1]
        g = as.character(g)
        save_s = append(save_s, g)}
      sid = metadata %>% filter(.data[[fact_var]] %in% save_s)
      snames = sid$sample_id
      metadata = filter(metadata, !sample_id %in% snames)
      counts_df = counts_df[,!(names(counts_df) %in% snames)]}  
    if(QC_check == "no"){
      if(log == "yes"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df[gene] = log2(combination_df[gene]+1)
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        ptab = compare_means(x ~ y, combination_df,method = "t.test")
        rownames(ptab) = NULL
        ptab = ptab[,-1]
      }
      if(log == "no"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        ptab = compare_means(x ~ y, combination_df,method = "t.test")
        model <- aov(x ~ y, data = combination_df)
        panova = summary(model)
        rownames(ptab) = NULL
        ptab = ptab[,-1]
      }}
    if(QC_check == "yes"){
      if(log == "yes"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df[gene] = log2(combination_df[gene]+1)
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        combination_df = filter(combination_df, quality_check == "pass")
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        ptab = compare_means(x ~ y, combination_df,method = "t.test")
        model <- aov(x ~ y, data = combination_df)
        panova = summary(model)
        rownames(ptab) = NULL
        ptab = ptab[,-1]
      } 
      if(log == "no"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        combination_df = filter(combination_df, quality_check == "pass")
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        ptab = compare_means(x ~ y, combination_df,method = "t.test")
        model <- aov(x ~ y, data = combination_df)
        panova = summary(model)
        rownames(ptab) = NULL
        ptab = ptab[,-1]
      }}
    if(is.null(save_s) != TRUE){
      for(i in 1:length(save_s)){
        vec = c(save_s[i], "*", "NA", "NA", "NA", "NA","NA")
        ptab = rbind(ptab, vec)}
      return(ptab)
    }else{
      return(ptab)
    }
  }}

table_panova = function(counts_df, metadata, log, QC_check, gene, fact_var, remove_samples = NULL){
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    keep_samples <- setdiff(metadata$sample_id, remove_samples)
    metadata <- metadata[metadata$sample_id %in% keep_samples, ]
    keep_cols <- c(intersect(colnames(counts_df), c("gene_id", "gene_name")), keep_samples)
    counts_df <- counts_df[, colnames(counts_df) %in% keep_cols, drop = FALSE]
  }
  if(length(unique(metadata[,fact_var]))<=1){
    ptab = data.frame(Error = c("P-value requires more than 1 grouping factor!"))
  } else {
    if(QC_check == "no"){
      if(log == "yes"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df[gene] = log2(combination_df[gene]+1)
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        model <- aov(x ~ y, data = combination_df)
        panova = summary(model)
        panova_df = as.data.frame(panova[[1]])
      }
      if(log == "no"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        model <- aov(x ~ y, data = combination_df)
        panova = summary(model)
        panova_df = as.data.frame(panova[[1]])
      }}
    if(QC_check == "yes"){
      if(log == "yes"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df[gene] = log2(combination_df[gene]+1)
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        combination_df = filter(combination_df, quality_check == "pass")
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        model <- aov(x ~ y, data = combination_df)
        panova = summary(model)
        panova_df = as.data.frame(panova[[1]])
      } 
      if(log == "no"){
        counts_table_logged = counts_df[,metadata$sample_id]
        transposed_counts = as.data.frame(t(counts_table_logged))
        colnames(transposed_counts) = counts_df$gene_name
        combination_df = cbind(metadata, transposed_counts)
        combination_df = combination_df[, !colnames(combination_df) %in% metadata_columns_to_remove, drop = FALSE]
        combination_df = combination_df[,unique(c(gene, fact_var, "sample_id", "condition", "replicate","quality_check"))]
        combination_df = filter(combination_df, quality_check == "pass")
        colnames(combination_df)[1] = "x"
        colnames(combination_df)[2] ="y"
        model <- aov(x ~ y, data = combination_df)
        panova = summary(model)
        panova_df = as.data.frame(panova[[1]])
      }}
    return(panova_df)
  }}

configure_plotly_panel = function(
    plotly_panel, 
    exportFormat = "svg", 
    annotationPosition = TRUE,
    annotationTail = TRUE, 
    annotationText = TRUE,
    axisTitleText = FALSE,
    colorbarPosition = TRUE,
    colorbarTitleText = TRUE,
    legendPosition = TRUE,
    legendText = TRUE,
    shapePosition = FALSE,
    titleText = TRUE) {
  
  plotly_panel = plotly_panel %>%
    plotly::config(
      toImageButtonOptions = list(
        format = exportFormat
      ),
      edits = list(
        annotationPosition = annotationPosition,
        annotationTail = annotationTail, 
        annotationText = annotationText,
        axisTitleText = axisTitleText,
        colorbarPosition = colorbarPosition,
        colorbarTitleText = colorbarTitleText,
        legendPosition = legendPosition,
        legendText = legendText,
        shapePosition = shapePosition,
        titleText = titleText
      )
    ) 
  return(plotly_panel)
}

draw_boxplot <- function(table, gene, factor, color_palette, log, box_type, plot_type){
  nb.cols <- length(unique(table[,factor]))
  my_colors <- colorRampPalette(brewer.pal(8, color_palette))(nb.cols)
  
  if(plot_type == "plotly"){
    if(grepl("-", gene) == TRUE){
      gene2 <- sub('-','_',gene)
      generow <- which(colnames(table) == gene)
      colnames(table)[generow] <- gene2
      gene2 <- as.symbol(gene2)
    } else {
      gene2 <- gene
    }
    if(factor == "replicate"){
      table$replicate <- as.character(table$replicate)
    }
    
    if(box_type == "Boxplot"){
      plotly::plot_ly(
        table,
        type = "box",
        x = as.formula(paste0("~", factor)),
        y = as.formula(paste0("~", gene2)),
        color = as.formula(paste0("~", factor)),
        colors = my_colors
      ) %>%
        layout(
          title = list(text = sprintf("<b>%s Gene Expression</b>", gene), x = 0.5, xanchor = "center"),
          xaxis = list(
            title = list(text = "<b>Sample Groups</b>"),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          yaxis = list(
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Gene Expression</b>")),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          legend = list(
            title = list(text = paste0('<b>', factor, '</b>')),
            x = 1.02,
            y = 1,
            xanchor = "left",
            bgcolor = 'rgba(255,255,255,0.9)',
            bordercolor = '#636363',
            borderwidth = 1
          )
        ) %>% configure_plotly_panel(exportFormat = "svg")
      
    } else if(box_type == "box_points"){
      plotly::plot_ly(
        table,
        type = "box",
        boxpoints = 'all',
        pointpos = 0,
        x = as.formula(paste0("~", factor)),
        y = as.formula(paste0("~", gene2)),
        color = as.formula(paste0("~", factor)),
        colors = my_colors
      ) %>%
        layout(
          title = list(text = sprintf("<b>%s Gene Expression</b>", gene), x = 0.5, xanchor = "center"),
          xaxis = list(
            title = list(text = "<b>Sample Groups</b>"),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          yaxis = list(
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Gene Expression</b>")),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          legend = list(
            title = list(text = paste0('<b>', factor, '</b>')),
            x = 1.02,
            y = 1,
            xanchor = "left",
            bgcolor = 'rgba(255,255,255,0.9)',
            bordercolor = '#636363',
            borderwidth = 1
          )
        ) %>% configure_plotly_panel(exportFormat = "svg")
      
    } else if(box_type == "Violin plot"){
      plotly::plot_ly(
        table,
        type = "violin",
        box = list(visible = TRUE),
        meanline = list(visible = TRUE),
        x = as.formula(paste0("~", factor)),
        y = as.formula(paste0("~", gene2)),
        color = as.formula(paste0("~", factor)),
        colors = my_colors
      ) %>%
        layout(
          title = list(text = sprintf("<b>%s Gene Expression</b>", gene), x = 0.5, xanchor = "center"),
          xaxis = list(
            title = list(text = "<b>Sample Groups</b>"),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          yaxis = list(
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Gene Expression</b>")),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          legend = list(
            title = list(text = paste0('<b>', factor, '</b>')),
            x = 1.02,
            y = 1,
            xanchor = "left",
            bgcolor = 'rgba(255,255,255,0.9)',
            bordercolor = '#636363',
            borderwidth = 1
          )
        ) %>% configure_plotly_panel(exportFormat = "svg")
      
    } else if(box_type == "violin_points"){
      plotly::plot_ly(
        table,
        type = "violin",
        points = 'all',
        pointpos = 0,
        box = list(visible = TRUE),
        meanline = list(visible = TRUE),
        x = as.formula(paste0("~", factor)),
        y = as.formula(paste0("~", gene2)),
        color = as.formula(paste0("~", factor)),
        colors = my_colors
      ) %>%
        layout(
          title = list(text = sprintf("<b>%s Gene Expression</b>", gene), x = 0.5, xanchor = "center"),
          xaxis = list(
            title = list(text = "<b>Sample Groups</b>"),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          yaxis = list(
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Gene Expression</b>")),
            zerolinecolor = '#969696',
            zerolinewidth = 2,
            linecolor = '#636363',
            linewidth = 2,
            gridcolor = '#e0e0e0'
          ),
          legend = list(
            title = list(text = paste0('<b>', factor, '</b>')),
            x = 1.02,
            y = 1,
            xanchor = "left",
            bgcolor = 'rgba(255,255,255,0.9)',
            bordercolor = '#636363',
            borderwidth = 1
          )
        ) %>% configure_plotly_panel(exportFormat = "svg")
    }
    
  } else if(plot_type == "ggplot"){
    if(factor == "replicate"){
      table$replicate <- as.character(table$replicate)
    }
    
    if(length(unique(table[,factor])) <= 4){
      comparisons <- get_factor_comparisons(condition = table[,factor])
      my_comparisons <- list()
      for(i in 1:nrow(comparisons)){
        v <- c(comparisons[i,2], comparisons[i,3])
        my_comparisons[[length(my_comparisons)+1]] <- v
      }
      
      if (box_type == "Boxplot"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_boxplot(notch = FALSE, alpha = 0.2, outlier.shape = NA) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE) +
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
        
      } else if(box_type == "box_points"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_boxplot(notch = FALSE, alpha = 0.2, outlier.shape = NA) +
          geom_jitter(size = 2, alpha = 0.6, width = 0.2) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE) +
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
        
      } else if(box_type == "Violin plot"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_violin(alpha = 0.2) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE) +
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
        
      } else if(box_type == "violin_points"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_violin(alpha = 0.2) +
          geom_jitter(size = 2, alpha = 0.6, width = 0.2) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE) +
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
      } 
      
    } else {
      if (box_type == "Boxplot"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_boxplot(notch = FALSE, alpha = 0.2, outlier.shape = NA) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
        
      } else if(box_type == "box_points"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_boxplot(notch = FALSE, alpha = 0.2, outlier.shape = NA) +
          geom_jitter(size = 2, alpha = 0.6, width = 0.2) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
        
      } else if(box_type == "Violin plot"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_violin(alpha = 0.2) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
        
      } else if(box_type == "violin_points"){
        ggplot(table, aes_string(x = factor, y = as.symbol(gene), color = factor, fill = factor)) +
          geom_violin(alpha = 0.2) +
          geom_jitter(size = 2, alpha = 0.6, width = 0.2) +
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Gene Expression")) + 
          xlab("Sample Groups") +
          ggtitle(paste(gene, "Gene Expression")) +
          scale_color_manual(values = my_colors, name = factor) + 
          scale_fill_manual(values = my_colors, name = factor) + 
          theme_classic() + 
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            legend.key = element_rect(fill = "white")
          )
      }
    }
  }
}

## ------------------ SAMPLE DISTANCE HEATMAP ------------------

run_sample_distance <- function(counts, metadata = NULL, remove_samples = NULL, scale_type = "none", ntop = 500) {
  
  gene_cols <- intersect(colnames(counts), c("gene_id", "gene_name"))
  
  sample_cols <- setdiff(colnames(counts), gene_cols)
  
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    sample_cols <- setdiff(sample_cols, remove_samples)
    if (!is.null(metadata)) {
      metadata <- metadata[metadata$sample_id %in% sample_cols, ]
    }
  }
  
  expr_matrix <- as.matrix(counts[, sample_cols, drop = FALSE])
  
  if ("gene_id" %in% colnames(counts)) {
    rownames(expr_matrix) <- counts$gene_id
  }
  
  if (scale_type == "row") {
    expr_matrix <- t(scale(t(expr_matrix)))
  } else if (scale_type == "column") {
    expr_matrix <- scale(expr_matrix)
  } else if (scale_type == "log") {
    expr_matrix <- log10(expr_matrix + 1)
  }
  
  rv <- rowVars(as.matrix(expr_matrix))
  
  select <- order(rv, decreasing = TRUE)
  
  sample_set <- colnames(expr_matrix)
  sampleDists <- dist(t(as.matrix(expr_matrix[select, ])))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- sample_set
  colnames(sampleDistMatrix) <- NULL
  
  sampleDistMatrix2 <- as.data.frame(sampleDistMatrix)
  rownames(sampleDistMatrix2) <- stringr::str_to_title(gsub("_", " ", rownames(sampleDistMatrix)))
  colnames(sampleDistMatrix2) <- stringr::str_to_title(gsub("_", " ", rownames(sampleDistMatrix)))
  
  return(sampleDistMatrix2)
}

plot_sample_distance_heatmap <- function(dist_matrix, 
                                         metadata = NULL,
                                         color_by = NULL,
                                         sidebar_color_scheme = "Set1",
                                         heatmap_title = "Sample Distance Heatmap",
                                         heatmap_color_scheme = "viridis",
                                         xlab = "",
                                         ylab = "",
                                         column_text_angle = 45,
                                         legend_title = "Euclidean\nsample\ndistance",
                                         cex_row = 1,
                                         cex_col = 1,
                                         show_dendrogram = c(TRUE, TRUE),
                                         plot_type = "heatmaply",
                                         colorbar_xpos = 1.02,
                                         colorbar_ypos = 0.5,
                                         show_tick_labels = c(TRUE, TRUE),
                                         colorbar_len = 0.4) {
  
  sampleDistMatrix <- dist_matrix
  
  row_side_colors <- NULL
  row_side_palette <- NULL
  
  if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
    formatted_names <- rownames(sampleDistMatrix)
    original_names <- tolower(gsub(" ", "_", formatted_names))
    
    ann_colors <- metadata[, c("sample_id", color_by), drop = FALSE]
    rownames(ann_colors) <- ann_colors$sample_id
    
    common_samples <- intersect(original_names, ann_colors$sample_id)
    if (length(common_samples) > 0) {
      sample_order <- match(original_names, common_samples)
      sample_order <- sample_order[!is.na(sample_order)]
      
      ann_colors_ordered <- ann_colors[common_samples[sample_order], , drop = FALSE]
      row_side_colors <- data.frame(ann_colors_ordered[, color_by, drop = FALSE])
      colnames(row_side_colors) <- stringr::str_to_title(gsub("_", " ", color_by))
      
      rownames(row_side_colors) <- formatted_names[original_names %in% common_samples]
      
      unique_vals <- unique(row_side_colors[[stringr::str_to_title(gsub("_", " ", color_by))]])
      n_colors <- length(unique_vals)
      
      if (n_colors <= 8) {
        base_colors <- RColorBrewer::brewer.pal(max(3, n_colors), sidebar_color_scheme)
      } else {
        base_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, sidebar_color_scheme))(n_colors)
      }
      
      row_side_palette <- setNames(base_colors[1:n_colors], unique_vals)
    }
  }
  
  selected_colors <- c()
  if (length(heatmap_color_scheme) == 1) {
    if (heatmap_color_scheme == "viridis") {
      selected_colors <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "magma") {
      selected_colors <- viridis::magma(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "inferno") {
      selected_colors <- viridis::inferno(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "plasma") {
      selected_colors <- viridis::plasma(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "cividis") {
      selected_colors <- viridis::cividis(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "rocket") {
      selected_colors <- viridis::rocket(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "mako") {
      selected_colors <- viridis::mako(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "turbo") {
      selected_colors <- viridis::turbo(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "RdYlBu") {
      selected_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(256)
    } else {
      tryCatch({
        selected_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, heatmap_color_scheme))(256)
      }, error = function(e) {
        selected_colors <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1)
      })
    }
  } else {
    selected_colors <- heatmap_color_scheme
  }
  
  if (tolower(plot_type) == "heatmaply") {
    p <- heatmaply::heatmaply(
      sampleDistMatrix,
      colors = selected_colors,
      xlab = xlab,
      ylab = ylab,
      main = heatmap_title,
      column_text_angle = column_text_angle,
      key.title = legend_title,
      cexRow = cex_row,
      cexCol = cex_col,
      Rowv = show_dendrogram[1],
      Colv = show_dendrogram[2],
      row_side_colors = row_side_colors,
      row_side_palette = row_side_palette,
      plot_method = "plotly",
      colorbar_xpos = colorbar_xpos,
      colorbar_ypos = colorbar_ypos,
      showticklabels = show_tick_labels,
      colorbar_len = colorbar_len
    )
    
    p <- p %>% plotly::layout(
      xaxis = list(title = xlab),
      yaxis = list(title = ylab),
      legend = list(
        orientation = "v",
        x = 1.15,
        y = 0.5,
        xanchor = "left",
        yanchor = "middle",
        bgcolor = 'rgba(255,255,255,0.9)',
        bordercolor = '#636363',
        borderwidth = 1
      )
    )
    
    if (exists("configure_plotly_panel")) {
      p <- configure_plotly_panel(p, exportFormat = "svg")
    }
    return(p)
    
  } else if (tolower(plot_type) == "ggplot") {
    
    hc_row <- NULL
    hc_col <- NULL
    
    if (show_dendrogram[1] || show_dendrogram[2]) {
      dist_mat <- as.matrix(sampleDistMatrix)
      
      if (show_dendrogram[1]) {
        hc_row <- hclust(as.dist(dist_mat))
      }
      if (show_dendrogram[2]) {
        hc_col <- hclust(as.dist(t(dist_mat)))
      }
    }
    
    dist_melt <- reshape2::melt(as.matrix(sampleDistMatrix))
    colnames(dist_melt) <- c("Sample1", "Sample2", "Distance")
    
    if (!is.null(hc_row)) {
      dist_melt$Sample1 <- factor(dist_melt$Sample1, 
                                  levels = rownames(sampleDistMatrix)[hc_row$order])
    } else {
      dist_melt$Sample1 <- factor(dist_melt$Sample1, 
                                  levels = rownames(sampleDistMatrix))
    }
    
    if (!is.null(hc_col)) {
      dist_melt$Sample2 <- factor(dist_melt$Sample2, 
                                  levels = colnames(sampleDistMatrix)[hc_col$order])
    } else {
      dist_melt$Sample2 <- factor(dist_melt$Sample2, 
                                  levels = colnames(sampleDistMatrix))
    }
    
    p <- ggplot2::ggplot(dist_melt, ggplot2::aes(x = Sample2, y = Sample1, fill = Distance)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradientn(colors = selected_colors, name = legend_title) +
      ggplot2::labs(title = heatmap_title, x = xlab, y = ylab) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = ggplot2::element_text(
          angle = column_text_angle, 
          vjust = 0.5, 
          size = 10
        ),
        axis.text.y = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 9),
        legend.position = "right",
        axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 12, face = "bold")
      )
    
    if (!show_tick_labels[2]) {
      p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                              axis.ticks.y = ggplot2::element_blank())
    }
    if (!show_tick_labels[1]) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                              axis.ticks.x = ggplot2::element_blank())
    }
    
    return(p)
  }
}

# ------------------ GENE EXPRESSION HEATMAP ------------------

prepare_gene_expression_matrix <- function(counts, 
                                           metadata, 
                                           gene_list, 
                                           remove_samples = NULL) {
  
  if ("gene_name" %in% colnames(counts)) {
    gene_data <- counts[counts$gene_name %in% gene_list, ]
    rownames(gene_data) <- gene_data$gene_name
    gene_data <- gene_data[, !colnames(gene_data) %in% c("gene_id", "gene_name"), drop = FALSE]
  } else {
    gene_data <- counts[rownames(counts) %in% gene_list, ]
  }
  
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    keep_cols <- setdiff(colnames(gene_data), remove_samples)
    gene_data <- gene_data[, keep_cols, drop = FALSE]
  }
  
  expr_matrix <- as.matrix(gene_data)
  
  expr_matrix[is.na(expr_matrix)] <- 0
  
  return(expr_matrix)
}

plot_gene_expression_heatmap <- function(counts,
                                         gene_list,
                                         gene_annotations,
                                         sample_list,
                                         metadata = NULL,
                                         color_by = NULL,
                                         sidebar_color_scheme = "Set1",
                                         heatmap_title = "Gene Expression Heatmap",
                                         heatmap_color_scheme = "RdYlBu",
                                         xlab = "",
                                         ylab = "Genes",
                                         column_text_angle = 90,
                                         legend_title = "Expression",
                                         cex_row = 0.5,
                                         cex_col = 0.5,
                                         cluster = "both",
                                         dendrogram = "both",
                                         scaling = "none",
                                         show_names = "both",
                                         heatmap_type = "heatmaply",
                                         vst_data = NULL) { 
  
  if ("gene_id" %in% colnames(counts) || "gene_name" %in% colnames(counts)) {
    gene_cols <- intersect(colnames(counts), c("gene_id", "gene_name"))
    
    data_to_use <- if (!is.null(vst_data)) vst_data else counts
    
    dat <- data_to_use[data_to_use$gene_name %in% gene_list | data_to_use$gene_id %in% gene_list, 
                       sample_list, drop = FALSE]
    rownames(dat) <- if("gene_name" %in% colnames(data_to_use)) {
      data_to_use$gene_name[data_to_use$gene_name %in% gene_list | data_to_use$gene_id %in% gene_list]
    } else {
      data_to_use$gene_id[data_to_use$gene_name %in% gene_list | data_to_use$gene_id %in% gene_list]
    }
  } else {
    data_to_use <- if (!is.null(vst_data)) vst_data else counts
    dat <- data_to_use[rownames(data_to_use) %in% gene_list, sample_list, drop = FALSE]
  }
  
  gene_sd <- apply(dat, 1, sd, na.rm = TRUE)
  dat <- dat[gene_sd != 0, , drop = FALSE]
  gene_list_filtered <- rownames(dat)
  
  col_side_colors <- NULL
  col_side_palette <- NULL
  
  if (!is.null(metadata) && !is.null(color_by) && color_by != "" && color_by %in% colnames(metadata)) {
    col_side_colors <- metadata[match(sample_list, metadata$sample_id), color_by, drop = FALSE]
    rownames(col_side_colors) <- sample_list
    colnames(col_side_colors) <- color_by
    
    unique_vals <- unique(col_side_colors[[color_by]])
    n_colors <- length(unique_vals)
    col_side_palette <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(min(8, max(3, n_colors)), sidebar_color_scheme))(n_colors),
      unique_vals
    )
  }
  
  gene_annotation_columns <- c("gene_name", "gene_id", "gene_biotype")
  available_cols <- intersect(gene_annotation_columns, colnames(gene_annotations))
  
  if (!is.null(gene_annotations) && length(available_cols) > 0) {
    dat_annotation <- gene_annotations[match(gene_list_filtered, gene_annotations$gene_name), 
                                       available_cols, drop = FALSE]
    if (nrow(dat_annotation) == 0 || all(is.na(dat_annotation$gene_name))) {
      dat_annotation <- gene_annotations[match(gene_list_filtered, gene_annotations$gene_id), 
                                         available_cols, drop = FALSE]
    }
    rownames(dat_annotation) <- gene_list_filtered
    
    hover_text <- sapply(sample_list, function(sid) {
      paste0("<b><i>", dat_annotation$gene_name, "</i></b><br>", 
             "<b>Sample: </b>", sid, "<br>",
             if("gene_id" %in% colnames(dat_annotation)) paste0("<b>ID: </b>", dat_annotation$gene_id, "<br>") else "",
             if("gene_biotype" %in% colnames(dat_annotation)) paste0("<b>Biotype: </b>", dat_annotation$gene_biotype, "<br>") else "")
    })
    if (length(gene_list_filtered) == 1) {
      hover_text <- matrix(hover_text, ncol = length(sample_list))
    }
  } else {
    hover_text <- NULL
  }
  
  show_rowv <- cluster %in% c("row", "both")
  show_colv <- cluster %in% c("column", "both")
  
  if (dendrogram == "none") {
    show_rowv <- FALSE
    show_colv <- FALSE
  } else if (dendrogram == "row") {
    show_rowv <- TRUE
    show_colv <- FALSE
  } else if (dendrogram == "column") {
    show_rowv <- FALSE
    show_colv <- TRUE
  } else if (dendrogram == "both") {
    show_rowv <- TRUE
    show_colv <- TRUE
  }
  
  if (show_names == "none") {
    show_tick_labels <- c(FALSE, FALSE)
  } else if (show_names == "column" || show_names == "x") {
    show_tick_labels <- c(TRUE, FALSE)
  } else if (show_names == "row" || show_names == "y") {
    show_tick_labels <- c(FALSE, TRUE)
  } else {
    show_tick_labels <- c(TRUE, TRUE)
  }
  
  # Convert scaling parameter
  scale_param <- if (tolower(scaling) %in% c("row", "z-score")) "row" 
  else if (tolower(scaling) == "log") "none"
  else "none"
  
  if (tolower(scaling) == "log") {
    dat <- log2(dat + 1)
    scale_param <- "none"
  }
  
  selected_colors <- c()
  if (length(heatmap_color_scheme) == 1) {
    if (heatmap_color_scheme == "viridis") {
      heatmap_colors <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "magma") {
      heatmap_colors <- viridis::magma(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "inferno") {
      heatmap_colors <- viridis::inferno(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "plasma") {
      heatmap_colors <- viridis::plasma(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "cividis") {
      heatmap_colors <- viridis::cividis(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "rocket") {
      heatmap_colors <- viridis::rocket(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "mako") {
      heatmap_colors <- viridis::mako(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "turbo") {
      heatmap_colors <- viridis::turbo(n = 256, alpha = 1, begin = 0, end = 1)
    } else if (heatmap_color_scheme == "RdYlBu") {
      heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(256)
    } else {
      tryCatch({
        heatmap_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, heatmap_color_scheme))(256)
      }, error = function(e) {
        heatmap_colors <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1)
      })
    }
  } else {
    heatmap_colors <- heatmap_color_scheme
  }
  
  if (tolower(heatmap_type) == "heatmaply") {
    hm <- heatmaply::heatmaply(
      dat,
      colors = heatmap_colors,
      xlab = xlab,
      ylab = ylab,
      main = heatmap_title,
      custom_hovertext = hover_text,
      column_text_angle = column_text_angle,
      key.title = legend_title,
      cexRow = cex_row,
      cexCol = cex_col,
      show_dendrogram = c(show_rowv, show_colv),
      col_side_colors = col_side_colors,
      col_side_palette = col_side_palette,
      plot_method = "plotly",
      showticklabels = show_tick_labels,
      colorbar_len = 0.3,
      side_color_colorbar_len = 0.25,
      Rowv = ifelse(nrow(dat) == 1, FALSE, show_rowv),
      Colv = ifelse(ncol(dat) == 1, FALSE, show_colv),
      scale = scale_param,
      width = 1500,
      height = 1000
    )
    
    hm <- hm %>% plotly::layout(
      legend = list(
        orientation = "v",
        x = 1.15,
        y = 0.5,
        xanchor = "left",
        yanchor = "middle",
        bgcolor = 'rgba(255,255,255,0.9)',
        bordercolor = '#636363',
        borderwidth = 1,
        tracegroupgap = 10
      )
    )
    
    hm <- configure_plotly_panel(hm, exportFormat = "svg")
    return(hm)
    
  } else if (tolower(heatmap_type) == "ggheatmap" || tolower(heatmap_type) == "ggplot") {
    p_heatmap <- heatmaply::ggheatmap(
      dat,
      xlab = xlab,
      ylab = ylab,
      main = heatmap_title,
      column_text_angle = column_text_angle,
      key.title = legend_title,
      cexRow = cex_row,
      cexCol = cex_col,
      show_dendrogram = c(show_rowv, show_colv),
      col_side_colors = col_side_colors,
      col_side_palette = col_side_palette,
      showticklabels = show_tick_labels,
      colorbar_len = 0.3,
      side_color_colorbar_len = 0.25,
      Rowv = ifelse(nrow(dat) == 1, FALSE, show_rowv),
      Colv = ifelse(ncol(dat) == 1, FALSE, show_colv),
      scale = scale_param
    )
    return(p_heatmap)
  } else {
    warning(paste0("Invalid heatmap type: ", heatmap_type))
    return(NULL)
  }
}
