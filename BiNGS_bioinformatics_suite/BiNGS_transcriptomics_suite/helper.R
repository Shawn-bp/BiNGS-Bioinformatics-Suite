
#------------------- PCA specific variables ------------------- 
pca_color_palette_list = c(
  "Set1 (Bright reds/blues)" = "Set1",
  "Dark2 (Darker tones)" = "Dark2",
  "Paired (Light & dark pairs)" = "Paired",
  "Accent (High contrast)" = "Accent",
  "Set2 (Pastel)" = "Set2",
  "Set3 (Very light)" = "Set3",
  "Pastel1 (Soft)" = "Pastel1",
  "Pastel2 (Muted)" = "Pastel2"
)

pca_metadata_columns_to_remove = c("species", "strandedness", "alignment_bam", "fastq_1", "fastq_2", "fastq_1_trimmed", "fastq_2_trimmed", "bings_sample_id", "bings_case_id", "bings_patient_id","file_name", "bings_block_id")

plotly_marker_symbols <- c(
  # --- Filled ---
  "circle", "square", "diamond", "triangle-up", "triangle-down",
  "triangle-left", "triangle-right", "pentagon", "hexagon", "hexagon2",
  "octagon", "star", "hexagram", "diamond-tall", "diamond-wide",
  "star-triangle-up", "star-triangle-down", "star-square", "star-diamond",
  "hourglass", "bowtie",
  
  # --- Filled with dot ---
  "circle-dot", "square-dot", "diamond-dot", "triangle-up-dot",
  "pentagon-dot", "hexagon-dot", "hexagon2-dot", "octagon-dot",
  "star-dot", "hexagram-dot", "diamond-tall-dot", "diamond-wide-dot",
  "hash-dot",
  
  # --- Diagonal triangles (filled) ---
  "triangle-ne", "triangle-se", "triangle-sw", "triangle-nw",
  
  # --- Open variants  ---
  "circle-open", "square-open", "diamond-open", "triangle-up-open",
  "triangle-down-open", "triangle-left-open", "triangle-right-open",
  "pentagon-open", "hexagon-open", "hexagon2-open", "octagon-open",
  "star-open", "hexagram-open", "diamond-tall-open", "diamond-wide-open",
  "star-triangle-up-open", "star-triangle-down-open",
  "star-square-open", "star-diamond-open",
  "hourglass-open", "bowtie-open",
  
  # --- Open with dot ---
  "circle-open-dot", "square-open-dot", "diamond-open-dot",
  "triangle-up-open-dot", "cross-open-dot", "x-open-dot",
  
  # --- Line/symbol shapes ---
  "cross", "cross-open", "cross-dot",
  "x", "x-open", "x-dot",
  "asterisk", "hash", "hash-open",
  "y-up", "y-up-open", "y-down", "y-down-open",
  "y-left", "y-left-open", "y-right", "y-right-open",
  "line-ew", "line-ns", "line-ne", "line-nw"
)
#------------------- General variables ------------------- 

color_palette_list = c(
  "Dark2", 
  "Set1",
  "Set2",
  "Set3"
)

QC_list = c(
  "Yes" ="yes",
  "No" = "no"
)

Log_list = c(
  "No" = "no",
  "Yes" ="yes"
)

type_list = c(
  "Interactive" = "plotly",
  "Static" ="ggplot"
)

#------------------- Boxplot specific variables -------------------
box_plot_list = c(
  "Boxplot" = "Boxplot",
  "Boxplot with points" = "box_points",
  "Violin plot" = "Violin plot",
  "Violin plot with pionts" = "violin_points"
)

metadata_columns_to_remove = c("species","sequencing_batch", "strandedness", "alignment_bam", "fastq_1", "fastq_2", "fastq_1_trimmed", "fastq_2_trimmed", "bings_sample_id", "bings_case_id", "bings_patient_id","file_name", "bings_block_id")

color_palette = "Set1"

log = "Yes"

QC_check = "No"

gene = "ST7"

fact_var = "condition"


#------------------- Heatmap specific variables ------------------- 
sample_distance_scaling_list = c("None" = "none", "Row" = "row", "Column" = "column", "Log scale" = "log")

gene_expression_scaling_list <- c("None" = "none", "Z-score" = "z-score", "Log2" = "log")

dendrogram_list = c("None" = "none", "Row" = "row", "Column" = "column", "Both" = "both")

sample_distance_show_names_list = c("None" = "none", "X-axis only" = "x", "Y-axis only" = "y", "Both" = "both")

gene_heatmap_show_names_list = c("None" = "none", "Gene Names" = "row", "Sample IDs" = "column", "Both" = "both")

heatmap_type_list = c("Interactive" = "heatmaply", "Static" = "ggplot")

heatmap_color_scheme_list = c(
  "Red-White-Blue" = "RdBu",
  "Red-Yellow-Blue" = "RdYlBu",
  "Spectral" = "Spectral",
  "Yellow-Blue-Purple" = "viridis"
)
