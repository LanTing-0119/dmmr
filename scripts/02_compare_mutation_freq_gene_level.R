# explor if there is no significant difference in the mutation rate of a given gene between the resistant and sensitive groups.
# Steps:
# 1. Generate a mutation matrix based on the MAF (Mutation Annotation Format) file.
# 2. Perform a t-test to compare the mutation frequencies between the primary resistant and sensitive groups.

#### Step 1: Load the Required Libraries ####
library(maftools)  # for handling MAF objects
library(dplyr)     # for data manipulation
library(tidyr)     # for reshaping data
library(ggplot2)   # for plotting

#### Step 2: Load the MAF File and Clinical Data ####
dmmr = read.maf(maf = 'outputs/dmmr.maf.gz', clinicalData = 'outputs/dmmr_clinical.csv')


#### Step 3: Convert MAF Data into a Mutation Matrix ####
# Extract the relevant columns: Patient_ID and Gene
mutation_data <- dmmr@data %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  mutate(Mutation_Present = 1)  # Mark mutation present with 1

# Step 3.1: Group by Sample_ID and Hugo_Symbol, and summarize (count mutations)
mutation_matrix <- mutation_data %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  summarize(Mutation_Present = n(), .groups = 'drop') %>%
  
  # Step 3.2: Spread to wide format (pivot_wider)
  pivot_wider(names_from = Hugo_Symbol, values_from = Mutation_Present, values_fill = list(Mutation_Present = 0))

# Step 3.3: Convert Mutation_Present to binary (1 if greater than 0, else 0)
mutation_matrix_binary <- mutation_matrix %>%
  mutate(across(-Tumor_Sample_Barcode, ~ ifelse(. > 0, 1, 0))) # across(-Sample_ID, ...) ensures that only the mutation columns (i.e., TP53, BRCA1, etc.) are transformed, not the Sample_ID column.

mutation_matrix
mutation_matrix_binary

#### Step 4: add clinical information
response_info = dmmr@clinical.data %>%
  select(Tumor_Sample_Barcode, ICI_response)
response_info$Tumor_Sample_Barcode = as.factor(response_info$Tumor_Sample_Barcode)
full_data = mutation_matrix %>% 
  left_join(response_info, by = "Tumor_Sample_Barcode") %>%
  select(Tumor_Sample_Barcode, ICI_response, sort(names(.)[!(names(.) %in% c("Tumor_Sample_Barcode", "ICI_response"))]))

full_data$ICI_response %>% table()

#### Step 5L Perform t-test/Wilcoxon
# Get the list of genes from the mutation matrix (excluding Tumor_Sample_Barcode and Response columns)
genes <- colnames(full_data)[-c(1,2)]

# Initialize a vector to store p-values for each gene
p_values <- numeric(length(genes))

# Perform t-test for each gene
for (gene in genes) {
  print(paste0(which(gene==genes),' / ',length(genes)))
  # Subset data for the current gene and perform t-test (resistant vs sensitive)
  gene_data <- full_data %>%
    select(Tumor_Sample_Barcode, ICI_response, gene) %>%
    filter(!is.na(ICI_response))  # Remove any missing responses
  # Perform t-test to compare mutation frequency between resistant and sensitive
  t_test_result <- t.test(gene_data[[gene]] ~ gene_data$ICI_response)
  # Store the p-value for the current gene
  p_values[gene] <- t_test_result$p.value
}
# Create a data frame to store p-values
t_test_results <- data.frame(Gene = genes, P_Value = p_values)

# Adjust p-values for multiple testing (Benjamini-Hochberg FDR)
t_test_results$FDR <- p.adjust(t_test_results$P_Value, method = "BH")
# View the results with FDR
View(t_test_results)




# Extract gene mutation data
gene_summary <- getGeneSummary(dmmr)

# View the gene mutation summary (this will give you mutation counts for each gene)
gene_summary
# Now, we need to create a binary matrix for each patient and gene.
# The `gene_summary` contains information on each mutation for each gene. We can transform this into a matrix.
