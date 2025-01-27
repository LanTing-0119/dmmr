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
  select(Group, everything()) %>% # reorder
  na.omit() %>%
  mutate(Group = as.factor(Group))


### 1.Feature Selection ###

library(Boruta)
# Perform feature selection
boruta_output <- Boruta(Group ~ ., data = mutation_matrix_dataset, doTrace = 2)
print(boruta_output)

# Get the final features selected by Boruta
confirmed_features <- getSelectedAttributes(boruta_output)
print(confirmed_features)

### 2.Building Classification Model ###
library(randomForest)
# Train a random forest model using selected features
rf_model <- randomForest(Group ~ ., data = mutation_matrix_dataset[, c(confirmed_features, "Group")])
print(rf_model)


### SVM
library(e1071)
# Train an SVM model
svm_model <- svm(Group ~ ., data = mutation_matrix_dataset[, c(confirmed_features, "Group")])
print(svm_model)

# 10-fold cross-validation using caret
library(caret)
train_control <- trainControl(method = "cv", number = 10)
model <- train(Group ~ ., data = mutation_matrix_dataset[, c(confirmed_features, "Group")], method = "rf", trControl = train_control)
print(model)
