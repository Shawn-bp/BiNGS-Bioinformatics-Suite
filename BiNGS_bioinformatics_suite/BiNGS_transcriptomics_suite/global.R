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
plot_pca <- function(pca_coords, variance, x_pc, y_pc, color_var = NULL,
                     palette = "Set1", plot_type = "ggplot") {
  
  x_label <- paste0(x_pc, " (", round(variance[as.numeric(gsub("PC", "", x_pc))], 2), "%)")
  y_label <- paste0(y_pc, " (", round(variance[as.numeric(gsub("PC", "", y_pc))], 2), "%)")
  
  if (plot_type == "ggplot") {
    p <- ggplot(pca_coords, aes_string(x = x_pc, y = y_pc, color = color_var)) +
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
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.key = element_rect(fill = "white")
      )
    return(p)
    
  } else if (plot_type == "plotly") {
    p <- plotly::plot_ly(
      pca_coords,
      x = ~get(x_pc),
      y = ~get(y_pc),
      color = if (!is.null(color_var)) pca_coords[[color_var]] else NULL,
      colors = RColorBrewer::brewer.pal(max(3, length(unique(pca_coords[[color_var]]))), palette),
      type = "scatter",
      mode = "markers",
      marker = list(size = 10),
      text = ~sample_id,
      hovertemplate = paste0(
        "<b>%{text}</b><br>",
        x_pc, ": %{x:.2f}<br>",
        y_pc, ": %{y:.2f}<br>",
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
        title = list(text = paste0('<b>', color_var, '</b>')),
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

# ------------------ SAMPLE DISTANCE HEATMAP ------------------
run_sample_distance <- function(counts, metadata = NULL, remove_samples = NULL, scale_type = "none") {
  
  expr <- counts[, !(colnames(counts) %in% c("gene_id", "gene_name"))]
  
  if (!is.null(remove_samples)) {
    expr <- expr[, !(colnames(expr) %in% remove_samples), drop = FALSE]
  }
  
  expr <- expr[apply(expr, 1, var) > 0, , drop = FALSE]
  
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
                                         color_scheme = "Dark2",
                                         cluster = "none",
                                         dendrogram = "none",
                                         show_names = "both",
                                         scaling = "none",
                                         color_by = NULL,
                                         heatmap_type = "ggplot") {
  
  dist_matrix <- as.matrix(dist_matrix)
  
  if (tolower(scaling) == "row") {
    dist_matrix <- t(scale(t(dist_matrix)))
  } else if (tolower(scaling) == "column") {
    dist_matrix <- scale(dist_matrix)
  } else if (tolower(scaling) == "log") {
    dist_matrix <- log10(dist_matrix + 1)
  }
  
  dist_matrix[is.na(dist_matrix)] <- 0
  
  ann_colors <- NULL
  if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
    ann_colors <- metadata[, c("sample_id", color_by), drop = FALSE]
    rownames(ann_colors) <- ann_colors$sample_id
    ann_colors <- ann_colors[rownames(dist_matrix), , drop = FALSE]
  }
  
  if (heatmap_type == "ggplot") {
    dist_melt <- reshape2::melt(dist_matrix)
    colnames(dist_melt) <- c("Sample1", "Sample2", "Distance")
    
    p <- ggplot(dist_melt, aes(x = Sample1, y = Sample2, fill = Distance)) +
      geom_tile() +
      scale_fill_distiller(palette = color_scheme, direction = 1, name = "Distance") +
      labs(
        title = "Sample Distance Heatmap", 
        x = "", 
        y = ""
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "right"
      )
    
    if (show_names == "none") {
      p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    } else if (show_names == "x") {
      p <- p + theme(axis.text.y = element_blank())
    } else if (show_names == "y") {
      p <- p + theme(axis.text.x = element_blank())
    }
    
    return(p)
  }
  
  else if (tolower(heatmap_type) == "heatmaply") {
    
    side_colors <- NULL
    side_palette <- NULL
    if (!is.null(ann_colors)) {
      side_colors <- data.frame(ann_colors[, color_by, drop = FALSE])
      colnames(side_colors) <- color_by
      rownames(side_colors) <- ann_colors$sample_id
      
      # Create named color palette for side colors
      unique_vals <- unique(side_colors[[color_by]])
      n_colors <- length(unique_vals)
      side_palette <- setNames(
        colorRampPalette(brewer.pal(min(8, max(3, n_colors)), "Set2"))(n_colors),
        unique_vals
      )
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
      labRow <- NULL
      labCol <- NULL
    } else if (show_names == "x") {
      labRow <- rownames(dist_matrix)
      labCol <- NULL
    } else if (show_names == "y") {
      labRow <- NULL
      labCol <- colnames(dist_matrix)
    } else {
      labRow <- rownames(dist_matrix)
      labCol <- colnames(dist_matrix)
    }
    
    hm <- heatmaply::heatmaply(
      dist_matrix,
      colors = colorRampPalette(brewer.pal(9, color_scheme))(256),
      Rowv = show_rowv,
      Colv = show_colv,
      labRow = labRow,
      labCol = labCol,
      row_side_colors = side_colors,
      col_side_colors = side_colors,
      row_side_palette = side_palette,
      col_side_palette = side_palette,
      main = "Sample Distance Heatmap",
      showticklabels = c(!is.null(labRow), !is.null(labCol)),
      colorbar_title = "Distance",
      width = 900,
      height = 700
    )
    
    hm <- hm %>% plotly::layout(
      legend = list(
        orientation = "v",
        x = 1.3,
        y = 0.5,
        xanchor = "left",
        yanchor = "middle",
        bgcolor = 'rgba(255,255,255,0.9)',
        bordercolor = '#636363',
        borderwidth = 1
      )
    )
    
    return(hm)
  }
}  

# ------------------ GENE EXPRESSION HEATMAP ------------------

prepare_gene_expression_matrix <- function(counts, 
                                           metadata, 
                                           gene_list, 
                                           remove_samples = NULL,
                                           scaling = "none") {
  
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
  
  if (tolower(scaling) == "row") {
    expr_matrix <- t(scale(t(expr_matrix)))
  } else if (tolower(scaling) == "column") {
    expr_matrix <- scale(expr_matrix)
  } else if (tolower(scaling) == "log") {
    expr_matrix <- log2(expr_matrix + 1)
  }
  
  expr_matrix[is.na(expr_matrix)] <- 0
  
  return(expr_matrix)
}

plot_gene_expression_heatmap <- function(expr_matrix,
                                         metadata = NULL,
                                         color_by = NULL,
                                         sidebar_color_scheme = "Set1",
                                         heatmap_title = "Gene Expression Heatmap",
                                         heatmap_color_scheme = "RdYlBu",
                                         scaling = "none",
                                         cluster = "both",
                                         dendrogram = "both",
                                         show_names = "both",
                                         heatmap_type = "ggplot") {
  
  expr_matrix <- as.matrix(expr_matrix)
  
  ann_colors <- NULL
  side_colors <- NULL
  
  if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
    ann_colors <- metadata[, c("sample_id", color_by), drop = FALSE]
    rownames(ann_colors) <- ann_colors$sample_id
    
    common_samples <- intersect(colnames(expr_matrix), ann_colors$sample_id)
    if (length(common_samples) > 0) {
      ann_colors <- ann_colors[common_samples, , drop = FALSE]
      side_colors <- data.frame(ann_colors[, color_by, drop = FALSE])
      colnames(side_colors) <- color_by
      rownames(side_colors) <- ann_colors$sample_id
    }
  }
  
  if (tolower(heatmap_type) == "ggplot") {
    expr_melt <- reshape2::melt(expr_matrix)
    colnames(expr_melt) <- c("Gene", "Sample", "Expression")
    
    p <- ggplot(expr_melt, aes(x = Sample, y = Gene, fill = Expression)) +
      geom_tile() +
      scale_fill_distiller(palette = heatmap_color_scheme, direction = -1, name = "Expression") +
      labs(title = heatmap_title, x = "", y = "") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "right"
      )
    
    if (show_names == "none") {
      p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    } else if (show_names == "column" || show_names == "x") {
      p <- p + theme(axis.text.y = element_blank())
    } else if (show_names == "row" || show_names == "y") {
      p <- p + theme(axis.text.x = element_blank())
    }
    
    return(p)
  }
  
  else if (tolower(heatmap_type) == "heatmaply") {
    
    side_colors <- NULL
    side_palette <- NULL
    if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
      ann_colors <- metadata[, c("sample_id", color_by), drop = FALSE]
      rownames(ann_colors) <- ann_colors$sample_id
      
      common_samples <- intersect(colnames(expr_matrix), ann_colors$sample_id)
      if (length(common_samples) > 0) {
        ann_colors <- ann_colors[common_samples, , drop = FALSE]
        side_colors <- data.frame(ann_colors[, color_by, drop = FALSE])
        colnames(side_colors) <- color_by
        rownames(side_colors) <- ann_colors$sample_id
        
        # Create named color palette for side colors
        unique_vals <- unique(side_colors[[color_by]])
        n_colors <- length(unique_vals)
        side_palette <- setNames(
          colorRampPalette(brewer.pal(min(8, max(3, n_colors)), sidebar_color_scheme))(n_colors),
          unique_vals
        )
      }
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
      labRow <- NULL
      labCol <- NULL
    } else if (show_names == "column" || show_names == "x") {
      labRow <- NULL
      labCol <- colnames(expr_matrix)
    } else if (show_names == "row" || show_names == "y") {
      labRow <- rownames(expr_matrix)
      labCol <- NULL
    } else {
      labRow <- rownames(expr_matrix)
      labCol <- colnames(expr_matrix)
    }
    
    heatmap_colors <- colorRampPalette(rev(brewer.pal(11, heatmap_color_scheme)))(256)

    hm <- heatmaply::heatmaply(
      expr_matrix,
      colors = heatmap_colors,
      Rowv = show_rowv,
      Colv = show_colv,
      labRow = labRow,
      labCol = labCol,
      main = heatmap_title,
      col_side_colors = side_colors,
      col_side_palette = side_palette,
      showticklabels = c(!is.null(labRow), !is.null(labCol)),
      colorbar_title = "Expression",
      width = 900,
      height = 700
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
    
    return(hm)
  }
}
