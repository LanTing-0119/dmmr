library(maftools)
library(ggplot2)
library(pheatmap)
library(caret)
library(cluster)
library(dplyr)
library(tidyr)

source('scripts/tools/02_calculate_mutation_type_scores.R')

# Import maf file
dmmr_all = read.maf(maf = 'data/vep_merge_oncokb.maf', clinicalData = 'outputs/dmmr_clinical_hg38.csv')
dmmr_all_filtered = subsetMaf(maf=dmmr_all, query = "ONCOGENIC %in% c('Inconclusive','Likely Neutral','Likely Oncogenic','Oncogenic','Resistance')")
response_info = dmmr_all@clinical.data %>%
  select(Tumor_Sample_Barcode, ICI_response) %>%
  filter(!is.na(ICI_response))

# define dmmr_resist & dmmr_sen
ICI_resistant_id = response_info %>% filter(ICI_response=="ICI-resistant") %>% pull(Tumor_Sample_Barcode) %>% as.character()
dmmr_resist = subsetMaf(maf = dmmr_all, tsb = ICI_resistant_id)
ICI_sensitive_id = response_info %>% filter(ICI_response=="ICI-sensitive") %>% pull(Tumor_Sample_Barcode) %>% as.character()
dmmr_sen = subsetMaf(maf = dmmr_all, tsb = ICI_sensitive_id)

dmmr=dmmr_all
# Generate Mutation Matrix
mutation_scores <- calculate_mutation_scores(dmmr@data)
# Create a binary matrix: 1 for mutation, 0 for no mutation
mutation_matrix <- mutation_scores %>%
  # Reshape the data: pivot_wider for reshaping
  pivot_wider(
    names_from = Hugo_Symbol,            # The new columns (genes)
    values_from = Total_Mutation_Score,  # The values (mutation scores)
    values_fill = list(Total_Mutation_Score = 0)  # Replace NA with 0
  )

# Convert mutation matrix to binary (1 if > 0, 0 if == 0)
mutation_matrix_binary = mutation_matrix %>%
  mutate(across(-Tumor_Sample_Barcode, ~ ifelse(. > 0, 1, 0)))

mutation_matrix_dataset <- mutation_matrix_binary %>%
  left_join(response_info, by = "Tumor_Sample_Barcode") %>%
  rename(Group=ICI_response) %>%
  select(-Tumor_Sample_Barcode) %>%
  select(Group, everything()) %>%
  na.omit()
mutation_matrix_dataset$Group <- as.factor(mutation_matrix_dataset$Group)





# Perform hierarchical clustering
dist_matrix <- dist(mutation_matrix[, -1])  # Exclude the 'Tumor_Sample_Barcode' column for clustering
hclust_result <- hclust(dist_matrix, method = "ward.D2")
# Plot hierarchical clustering dendrogram
plot(hclust_result, main = "Hierarchical Clustering of Samples Based on Mutations", 
     xlab = "", sub = "", cex = 0.9)


# Load necessary libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Step 1: Preprocess the data
# Extract the features (exclude the Group column) and convert to a matrix
data_matrix <- mutation_matrix_dataset %>%
  select(-Group) %>%
  as.matrix()

data_matrix_t = t(data_matrix)

# Extract the Group column as a factor to use for grouping
group_info <- mutation_matrix_dataset$Group

# Step 2: Check that the number of columns in data_matrix and length of group_info match
if (length(group_info) != ncol(data_matrix_t)) {
  stop("Number of samples (columns) in the data matrix and the group information do not match.")
}

# Step 3: Ensure that the group information is a factor with correct levels
group_info <- factor(group_info, levels = c("ICI-resistant", "ICI-sensitive"))

# Step 4: Create a color palette for the Group variable (for annotations)
group_colors <- c("ICI-resistant" = "red", "ICI-sensitive" = "blue")

# Step 5: Create the heatmap with annotations
Heatmap(
  data_matrix,                           # The matrix of data to plot
  name = "Mutation Status",               # Label for the heatmap
  row_title = "Samples",                    # Title for rows (genes)
  column_title = "Genes",               # Title for columns (samples)
  left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:3),
                                                   labels = c("group1", "group2"), 
                                                   labels_gp = gpar(col = "white", fontsize = 10))),
  row_km = 2,
  show_column_names = TRUE,               # Show column names (sample names)
  show_row_names = FALSE,                 # Optionally hide row names (gene names)
  clustering_distance_rows = "euclidean", # Row clustering method
  clustering_distance_columns = "euclidean", # Column clustering method
  clustering_method_rows = "complete",    # Row clustering method
  clustering_method_columns = "complete", # Column clustering method
  heatmap_legend_param = list(title = "Mutation", at = c(0, 1), labels = c("No Mutation", "Mutation")), # Legend details
  column_names_gp = gpar(fontsize = 10),  # Font size for column names
  row_names_gp = gpar(fontsize = 10)      # Font size for row names
)


