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

gene_expression_scaling_list = c("None" = "none", "Row" = "row", "Log scale" = "log")

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
