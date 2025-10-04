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
  "Plotly" = "plotly",
  "GGplot" ="ggplot"
)

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


#Heatmap specific variables
scaling_list = c("None" = "none", "Row" = "row", "Column" = "column", "Log scale" = "log")

clustering_list = c("None" = "none", "Row" = "row", "Column" = "column", "Both" = "both")

dendrogram_list = c("None" = "none", "Row" = "row", "Column" = "column", "Both" = "both")

show_names_list = c("None" = "none", "Row" = "row", "Column" = "column", "Both" = "both")

heatmap_type_list = c("Interactive" = "heatmaply", "Static" = "ggplot")
