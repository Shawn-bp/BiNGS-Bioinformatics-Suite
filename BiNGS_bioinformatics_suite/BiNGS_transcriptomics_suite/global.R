
# ------------------ PCA ------------------
# Run PCA
run_pca <- function(counts, metadata, scale_data = TRUE, remove_samples = NULL) {
  
  # Remove specified samples
  if (!is.null(remove_samples) && length(remove_samples) > 0) {
    keep_cols <- setdiff(colnames(counts), c("gene_id", "gene_name", remove_samples))
    gene_cols <- intersect(colnames(counts), c("gene_id", "gene_name"))
    counts <- counts[, c(gene_cols, keep_cols), drop = FALSE]
    metadata <- metadata[metadata$sample_id %in% keep_cols, ]
  }
  
  expr <- counts[, !(colnames(counts) %in% c("gene_id", "gene_name"))]
  expr_t <- t(expr)
  expr_t <- expr_t[, apply(expr_t, 2, var) > 0, drop = FALSE]
  
  pca_res <- prcomp(expr_t, scale. = scale_data)
  
  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  
  pca_coords <- as.data.frame(pca_res$x)
  pca_coords$sample_id <- rownames(pca_coords)
  
  if (!is.null(metadata)) {
    pca_coords <- dplyr::left_join(pca_coords, metadata, by = "sample_id")
  }
  
  return(list(
    pca = pca_res,
    coords = pca_coords,
    variance = var_explained
  ))
}

# Plot PCA
plot_pca <- function(pca_coords, variance, x_pc, y_pc, color_var = NULL, shape_var = NULL,
                     palette = "Set1", plot_type = "ggplot") {
  
  x_label <- paste0(x_pc, " (", round(variance[as.numeric(gsub("PC", "", x_pc))], 2), "%)")
  y_label <- paste0(y_pc, " (", round(variance[as.numeric(gsub("PC", "", y_pc))], 2), "%)")
  
  if (plot_type == "ggplot") {
    p <- ggplot(pca_coords, aes_string(x = x_pc, y = y_pc, color = color_var, shape = shape_var)) +
      geom_point(size = 4, alpha = 0.8) +
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
        axis.title.x = element_text(size = 40, face = "bold"),
        axis.title.y = element_text(size = 40, face = "bold"),
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.key = element_rect(fill = "white"),
      )
    
    # Add shape scale if shape_var is provided
    if (!is.null(shape_var) && shape_var != "") {
      p <- p + scale_shape_discrete(name = shape_var)
    }
    
    return(p)
    
  } else if (plot_type == "plotly") {
    
    # Prepare symbol mapping for plotly if shape_var is provided
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
      marker = list(size = 10),
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
  
  # remove identical condition comparisons
  comparisons = comparisons[which(comparisons$comparison_group != comparisons$reference_group), ]
  comparisons = comparisons[order(comparisons$reference_group), , drop = FALSE]
  
  # remove reverse comparisons
  comparisons$comparison_name = paste0(comparisons$comparison_group, "_vs_", comparisons$reference_group)
  comparisons = comparisons[, c("comparison_name", "comparison_group", "reference_group"), drop = FALSE]
  rownames(comparisons) = NULL
  
  return(comparisons)
}

modify_df <- function(counts_df, metadata, log, QC_check, gene, remove_samples = NULL){
  #Remove specified samples
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
  # Remove specified samples
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
  # Remove specified samples
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
  # Remove specified samples
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
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression</b>")),
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
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression</b>")),
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
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression</b>")),
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
            title = list(text = paste0("<b>", ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression</b>")),
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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
          ylab(paste0(ifelse(log == "yes", "Log₂ ", ""), "Normalized Gene Expression")) + 
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

# Calculate distance matrix (used by reactive)
run_sample_distance <- function(counts, metadata = NULL, remove_samples = NULL, scale_type = "none", ntop = NULL) {
  
  expr <- counts[, !(colnames(counts) %in% c("gene_id", "gene_name"))]
  
  if (!is.null(remove_samples)) {
    expr <- expr[, !(colnames(expr) %in% remove_samples), drop = FALSE]
  }
  
  expr <- expr[apply(expr, 1, var) > 0, , drop = FALSE]
  
  # Select top variable genes if specified
  if (!is.null(ntop) && ntop < nrow(expr)) {
    gene_vars <- apply(expr, 1, var)
    select <- order(gene_vars, decreasing = TRUE)[1:min(ntop, nrow(expr))]
    expr <- expr[select, , drop = FALSE]
  }
  
  # Apply scaling
  if (scale_type == "row") {
    expr <- t(scale(t(expr)))
  } else if (scale_type == "column") {
    expr <- scale(expr)
  } else if (scale_type == "log") {
    expr <- log10(expr + 1)
  }
  
  dist_matrix <- dist(t(expr), method = "euclidean")
  dist_matrix <- as.matrix(dist_matrix)
  
  rownames(dist_matrix) <- colnames(expr)
  colnames(dist_matrix) <- colnames(expr)
  
  return(dist_matrix)
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
                                         show_dendrogram = c(TRUE, TRUE),
                                         plot_type = "heatmaply",
                                         colorbar_xpos = 1.02,
                                         colorbar_ypos = 0.5,
                                         show_tick_labels = c(TRUE, TRUE),
                                         colorbar_len = 0.4) {
  
  # Use the pre-calculated distance matrix
  sampleDistMatrix <- dist_matrix
  
  # Get sample names from the matrix
  sample_names <- rownames(sampleDistMatrix)
  
  # Prepare side colors if metadata and color_by are provided
  row_side_colors <- NULL
  row_side_palette <- NULL
  
  if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
    ann_colors <- metadata[, c("sample_id", color_by), drop = FALSE]
    rownames(ann_colors) <- ann_colors$sample_id
    
    common_samples <- intersect(sample_names, ann_colors$sample_id)
    if (length(common_samples) > 0) {
      ann_colors <- ann_colors[common_samples, , drop = FALSE]
      row_side_colors <- data.frame(ann_colors[, color_by, drop = FALSE])
      colnames(row_side_colors) <- color_by
      rownames(row_side_colors) <- ann_colors$sample_id
      
      # Create named color palette for side colors
      unique_vals <- unique(row_side_colors[[color_by]])
      n_colors <- length(unique_vals)
      row_side_palette <- setNames(
        colorRampPalette(RColorBrewer::brewer.pal(min(8, max(3, n_colors)), sidebar_color_scheme))(n_colors),
        unique_vals
      )
    }
  }
  
  # Select color scheme
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
      # Try as RColorBrewer palette
      selected_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, heatmap_color_scheme))(256)
    }
  } else {
    selected_colors <- heatmap_color_scheme
  }
  
  if (tolower(plot_type) == "heatmaply") {
    # Create heatmaply plot
    p <- heatmaply::heatmaply(
      sampleDistMatrix,
      colors = selected_colors,
      xlab = xlab,
      ylab = ylab,
      main = heatmap_title,
      column_text_angle = column_text_angle,
      key.title = "Distance",
      Rowv = show_dendrogram[1],
      Colv = show_dendrogram[2],
      row_side_colors = row_side_colors,
      row_side_palette = row_side_palette,
      plot_method = "plotly",
      colorbar_xpos = colorbar_xpos,
      colorbar_ypos = colorbar_ypos,
      showticklabels = show_tick_labels,
      colorbar_len = colorbar_len,
      width = 900,
      height = 700
    )
    
    # Enhance layout
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
    
    p <- configure_plotly_panel(p, exportFormat = "svg")
    return(p)
    
  } else if (tolower(plot_type) == "ggplot") {
    # Create ggplot version
    dist_melt <- reshape2::melt(sampleDistMatrix)
    colnames(dist_melt) <- c("Sample1", "Sample2", "Distance")
    
    p <- ggplot(dist_melt, aes(x = Sample2, y = Sample1, fill = Distance)) +
      geom_tile() +
      scale_fill_gradientn(colors = selected_colors, name = "Distance") +
      labs(title = heatmap_title, x = xlab, y = ylab) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = column_text_angle, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "right"
      )
    
    if (!show_tick_labels[1]) {
      p <- p + theme(axis.text.y = element_blank())
    }
    if (!show_tick_labels[2]) {
      p <- p + theme(axis.text.x = element_blank())
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

plot_gene_expression_heatmap <- function(expr_matrix,
                                         metadata = NULL,
                                         color_by = NULL,
                                         sidebar_color_scheme = "Set1",
                                         heatmap_title = "Gene Expression Heatmap",
                                         heatmap_color_scheme = "RdYlBu",
                                         cluster = "both",
                                         dendrogram = "both",
                                         scaling = "none",
                                         show_names = "both",
                                         heatmap_type = "ggplot") {
  
  expr_matrix <- as.matrix(expr_matrix)
  
  # Prepare annotation colors
  annotation_df <- NULL
  row_side_palette <- NULL
  
  if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
    ann_colors <- metadata[, c("sample_id", color_by), drop = FALSE]
    rownames(ann_colors) <- ann_colors$sample_id
    
    common_samples <- intersect(colnames(expr_matrix), ann_colors$sample_id)
    if (length(common_samples) > 0) {
      annotation_df <- ann_colors[common_samples, , drop = FALSE]
      colnames(annotation_df) <- c("sample_id", "group")
      
      # Create color palette
      unique_vals <- unique(annotation_df$group)
      n_colors <- length(unique_vals)
      row_side_palette <- setNames(
        colorRampPalette(RColorBrewer::brewer.pal(min(8, max(3, n_colors)), sidebar_color_scheme))(n_colors),
        unique_vals
      )
    }
  }
  
  # ==== HEATMAPLY VERSION ====
  if (tolower(heatmap_type) == "heatmaply") {
    
    # Prepare side colors for heatmaply
    row_side_colors <- NULL
    
    if (!is.null(annotation_df)) {
      row_side_colors <- data.frame(annotation_df[, "group", drop = FALSE])
      colnames(row_side_colors) <- color_by
      rownames(row_side_colors) <- annotation_df$sample_id
    }
    
    # Select color scheme
    selected_colors <- c()
    if (length(heatmap_color_scheme) == 1) {
      if (heatmap_color_scheme == "viridis") {
        selected_colors <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1)
      } else if (heatmap_color_scheme == "magma") {
        selected_colors <- viridis::magma(n = 256, alpha = 1, begin = 0, end = 1)
      } else if (heatmap_color_scheme == "RdYlBu") {
        selected_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(256)
      } else {
        selected_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, heatmap_color_scheme))(256)
      }
    } else {
      selected_colors <- heatmap_color_scheme
    }
    
    # Apply scaling
    if (tolower(scaling) == "row") {
      expr_matrix <- t(scale(t(expr_matrix)))
      expr_matrix[is.na(expr_matrix)] <- 0
      expr_matrix[is.infinite(expr_matrix)] <- 0
    } else if (tolower(scaling) == "column") {
      expr_matrix <- scale(expr_matrix)
      expr_matrix[is.na(expr_matrix)] <- 0
      expr_matrix[is.infinite(expr_matrix)] <- 0
    } else if (tolower(scaling) == "log") {
      expr_matrix <- log2(expr_matrix + 1)
    }
    
    # Determine dendrograms
    show_row_dend <- dendrogram %in% c("row", "both")
    show_col_dend <- dendrogram %in% c("column", "both")
    
    # Determine tick labels
    show_row_labels <- show_names %in% c("row", "y", "both")
    show_col_labels <- show_names %in% c("column", "x", "both")
    
    # Create heatmaply
    p <- heatmaply::heatmaply(
      expr_matrix,
      colors = selected_colors,
      main = heatmap_title,
      xlab = "",
      ylab = "",
      Rowv = if (cluster %in% c("row", "both")) TRUE else FALSE,
      Colv = if (cluster %in% c("column", "both")) TRUE else FALSE,
      dendrogram = dendrogram,
      col_side_colors = row_side_colors,
      col_side_palette = row_side_palette,
      showticklabels = c(show_row_labels, show_col_labels),
      key.title = "Expression",
      plot_method = "plotly",
      width = 900,
      height = 700
    )
    
    # Enhance layout
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
    
    p <- configure_plotly_panel(p, exportFormat = "svg")
    return(p)
    
  }
  
  # ==== GGPLOT VERSION ====
  else if (tolower(heatmap_type) == "ggplot") {
    
    # Perform clustering on UNSCALED data
    row_order <- rownames(expr_matrix)
    col_order <- colnames(expr_matrix)
    row_hclust <- NULL
    col_hclust <- NULL
    
    if (cluster %in% c("row", "both") || dendrogram %in% c("row", "both")) {
      tryCatch({
        # Remove any rows with NA, Inf, or zero variance
        valid_rows <- apply(expr_matrix, 1, function(x) {
          !any(is.na(x)) && !any(is.infinite(x)) && var(x, na.rm = TRUE) > 0
        })
        
        if (sum(valid_rows) > 1) {
          row_dist <- dist(expr_matrix[valid_rows, ])
          row_hclust <- hclust(row_dist)
          row_order <- rownames(expr_matrix)[valid_rows][row_hclust$order]
          # Add back any invalid rows at the end
          if (sum(!valid_rows) > 0) {
            row_order <- c(row_order, rownames(expr_matrix)[!valid_rows])
          }
        }
      }, error = function(e) {
        warning("Row clustering failed: ", e$message)
      })
    }
    
    if (cluster %in% c("column", "both") || dendrogram %in% c("column", "both")) {
      tryCatch({
        # Remove any columns with NA, Inf, or zero variance
        valid_cols <- apply(expr_matrix, 2, function(x) {
          !any(is.na(x)) && !any(is.infinite(x)) && var(x, na.rm = TRUE) > 0
        })
        
        if (sum(valid_cols) > 1) {
          col_dist <- dist(t(expr_matrix[, valid_cols]))
          col_hclust <- hclust(col_dist)
          col_order <- colnames(expr_matrix)[valid_cols][col_hclust$order]
          # Add back any invalid columns at the end
          if (sum(!valid_cols) > 0) {
            col_order <- c(col_order, colnames(expr_matrix)[!valid_cols])
          }
        }
      }, error = function(e) {
        warning("Column clustering failed: ", e$message)
      })
    }
    
    # Apply scaling
    if (tolower(scaling) == "row" || tolower(scaling) == "z-score") {
      expr_matrix <- t(scale(t(expr_matrix)))
      expr_matrix[is.na(expr_matrix)] <- 0
      expr_matrix[is.infinite(expr_matrix)] <- 0
    } else if (tolower(scaling) == "column") {
      expr_matrix <- scale(expr_matrix)
      expr_matrix[is.na(expr_matrix)] <- 0
      expr_matrix[is.infinite(expr_matrix)] <- 0
    } else if (tolower(scaling) == "log") {
      expr_matrix <- log2(expr_matrix + 1)
    }
    
    # Reorder matrix based on clustering
    expr_matrix <- expr_matrix[row_order, col_order]
    
    # Melt for ggplot
    expr_melt <- reshape2::melt(expr_matrix)
    colnames(expr_melt) <- c("Gene", "Sample", "Expression")
    
    # Set factor levels to maintain order
    expr_melt$Gene <- factor(expr_melt$Gene, levels = row_order)
    expr_melt$Sample <- factor(expr_melt$Sample, levels = col_order)
    
    # Create base heatmap
    p_heat <- ggplot(expr_melt, aes(x = Sample, y = Gene, fill = Expression)) +
      geom_tile() +
      scale_fill_distiller(palette = heatmap_color_scheme, direction = -1, name = "Expression") +
      labs(title = heatmap_title, x = "", y = "") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "right"
      )
    
    # Add color bar annotation 
    if (!is.null(annotation_df) && !is.null(row_side_palette)) {
      # Match annotation to ordered samples
      col_annotation <- annotation_df[match(col_order, annotation_df$sample_id), ]
      col_annotation <- col_annotation[!is.na(col_annotation$sample_id), ]
      col_annotation$Sample <- factor(col_annotation$sample_id, levels = col_order)
      col_annotation$y_pos <- length(row_order) + 0.5
      
      # Add annotation bar at top
      p_heat <- p_heat +
        geom_tile(data = col_annotation,
                  aes(x = Sample, y = y_pos, fill = NULL, color = group),
                  height = 0.3, size = 3, inherit.aes = FALSE, show.legend = TRUE) +
        scale_color_manual(values = row_side_palette, name = color_by) +
        guides(
          fill = guide_colorbar(order = 1, title = "Expression"),
          color = guide_legend(order = 2, title = color_by)
        )
    }
    
    # Handle show_names
    if (show_names == "none") {
      p_heat <- p_heat + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    } else if (show_names == "column" || show_names == "x") {
      p_heat <- p_heat + theme(axis.text.y = element_blank())
    } else if (show_names == "row" || show_names == "y") {
      p_heat <- p_heat + theme(axis.text.x = element_blank())
    }
    
    # Create dendrograms if needed
    p_row_dend <- NULL
    p_col_dend <- NULL
    
    if (dendrogram %in% c("row", "both") && !is.null(row_hclust)) {
      library(ggdendro)
      row_dend_data <- dendro_data(as.dendrogram(row_hclust), type = "rectangle")
      
      p_row_dend <- ggplot(segment(row_dend_data)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_x_continuous(expand = c(0, 0.5)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_void()
    }
    
    if (dendrogram %in% c("column", "both") && !is.null(col_hclust)) {
      library(ggdendro)
      col_dend_data <- dendro_data(as.dendrogram(col_hclust), type = "rectangle")
      
      p_col_dend <- ggplot(segment(col_dend_data)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        scale_x_continuous(expand = c(0, 0.5)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_void()
    }
    
    # Store the legend separately before combining
    if (!is.null(annotation_df) && !is.null(row_side_palette)) {
      library(cowplot)
      
      # Extract legend
      legend <- get_legend(p_heat)
      
      # Remove legend from main plot for combining
      p_heat_no_legend <- p_heat + theme(legend.position = "none")
      
      # Combine plots with dendrograms
      if (!is.null(p_row_dend) && !is.null(p_col_dend)) {
        library(gridExtra)
        empty <- ggplot() + theme_void()
        p_combined <- grid.arrange(
          p_col_dend, empty, empty,
          p_heat_no_legend, p_row_dend, legend,
          ncol = 3, nrow = 2,
          widths = c(4, 0.5, 0.5),
          heights = c(0.5, 4)
        )
      } else if (!is.null(p_row_dend)) {
        library(gridExtra)
        p_combined <- grid.arrange(
          p_heat_no_legend, p_row_dend, legend,
          ncol = 3,
          widths = c(4, 0.5, 0.5)
        )
      } else if (!is.null(p_col_dend)) {
        library(gridExtra)
        empty <- ggplot() + theme_void()
        p_combined <- grid.arrange(
          p_col_dend, empty,
          p_heat_no_legend, legend,
          ncol = 2, nrow = 2,
          widths = c(4, 0.5),
          heights = c(0.5, 4)
        )
      } else {
        p_combined <- p_heat
      }
      return(p_combined)
      
    } else {

      if (!is.null(p_row_dend) && !is.null(p_col_dend)) {
        library(gridExtra)
        empty <- ggplot() + theme_void()
        p_combined <- grid.arrange(
          p_col_dend, empty,
          p_heat, p_row_dend,
          ncol = 2, nrow = 2,
          widths = c(4, 0.5),
          heights = c(0.5, 4)
        )
        return(p_combined)
      } else if (!is.null(p_row_dend)) {
        library(gridExtra)
        p_combined <- grid.arrange(
          p_heat, p_row_dend,
          ncol = 2,
          widths = c(4, 0.5)
        )
        return(p_combined)
      } else if (!is.null(p_col_dend)) {
        library(gridExtra)
        p_combined <- grid.arrange(
          p_col_dend,
          p_heat,
          nrow = 2,
          heights = c(0.5, 4)
        )
        return(p_combined)
      } else {
        return(p_heat)
      }
    }
  }
}
