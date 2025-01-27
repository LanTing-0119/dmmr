# Load necessary libraries
library(paradox)
library(mlr3)
library(mlr3tuning)
library(mlr3proba)

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

colnames(mutation_matrix_dataset) = sub("-", "_", colnames(mutation_matrix_dataset))


# Assuming 'train_set' contains your features and target, with 'mutation' as the target column
# Replace 'train_set' with the actual dataset variable

# Define the classification task
task_classif <- TaskClassif$new(id = "ClassifTask", backend = mutation_matrix_dataset, target = "Group")

# Define the hyperparameter search space for classic.ranger
search_space_classif <- ps(
  num.trees = p_int(lower = 100, upper = 1000),
  mtry = p_int(lower = 1, upper = floor(sqrt(ncol(mutation_matrix_dataset)))),
  min.node.size = p_int(lower = 1, upper = 20),
  max.depth = p_int(lower = 1, upper = 10),
  sample.fraction = p_dbl(lower = 0.5, upper = 1),
  replace = p_lgl(),
  splitrule = p_fct(levels = c("gini", "extratrees", "hellinger")),
  importance = p_fct(levels = c("none", "impurity", "permutation")),
  num.threads = p_int(lower = 1, upper = parallel::detectCores(logical = FALSE)),
  verbose = p_lgl()
)

# Create the learner (ranger classifier)
learner_classif <- lrn("classif.ranger", predict_type = "response")

# Set resampling strategy: 5-fold cross-validation
resampling_classif <- rsmp("cv", folds = 5)

# Set performance measure: Accuracy or AUC (depending on your task)
#measure_classif <- msr("classif.acc")  # You can also use msr("classif.auc") for AUC
measure_classif <- msr("classif.recall") # Set performance measure: prioritize recall (sensitivity)
# Define the termination criterion: stop after 50 evaluations
termination_classif <- trm("evals", n_evals = 150)

# Set up the tuning instance for classification
tuning_instance_classif <- TuningInstanceSingleCrit$new(
  task = task_classif,
  learner = learner_classif,
  resampling = resampling_classif,
  measure = measure_classif,
  search_space = search_space_classif,
  terminator = termination_classif
)

# Choose a tuner: Random Search
tuner_classif <- tnr("random_search")

# Tune the model
tuner_classif$optimize(tuning_instance_classif)

# Extract best hyperparameter settings
best_params_classif <- tuning_instance_classif$result_learner_param_vals

# Print the best model parameters
print(best_params_classif)

# Train the model on the entire dataset using the best parameters
learner_classif$param_set$values <- best_params_classif
learner_classif$train(task_classif)

# Make predictions on the training set
predictions <- learner_classif$predict(task_classif)

# Compute the confusion matrix
conf_matrix <- predictions$confusion
print(conf_matrix)












# Define the classification task
task_classif <- TaskClassif$new(id = "ClassifTask", backend = mutation_matrix_dataset, target = "Group")

# Define the hyperparameter search space for ranger
search_space_classif <- ps(
  num.trees = p_int(lower = 100, upper = 1000),
  mtry = p_int(lower = 1, upper = floor(sqrt(ncol(mutation_matrix_dataset)))),
  min.node.size = p_int(lower = 1, upper = 20),
  max.depth = p_int(lower = 1, upper = 10),
  sample.fraction = p_dbl(lower = 0.5, upper = 1),
  replace = p_lgl(),
  splitrule = p_fct(levels = c("gini", "extratrees", "hellinger")),
  importance = p_fct(levels = c("none", "impurity", "permutation")),
  num.threads = p_int(lower = 1, upper = parallel::detectCores(logical = FALSE)),
  verbose = p_lgl()
)

# Create the learner (ranger classifier)
learner_classif <- lrn("classif.ranger", predict_type = "response")

# Set resampling strategy: 5-fold cross-validation
resampling_classif <- rsmp("cv", folds = 5)

# Set performance measure: Recall (sensitivity)
measure_classif <- msr("classif.recall")

# Define the termination criterion: stop after 150 evaluations
termination_classif <- trm("evals", n_evals = 20)

# Set up the tuning instance for classification
tuning_instance_classif <- TuningInstanceSingleCrit$new(
  task = task_classif,
  learner = learner_classif,
  resampling = resampling_classif,
  measure = measure_classif,
  search_space = search_space_classif,
  terminator = termination_classif
)

# Choose a tuner: Random Search
tuner_classif <- tnr("random_search")

# Perform the optimization (tuning)
tuner_classif$optimize(tuning_instance_classif)

# Extract the best hyperparameter settings
best_params_classif <- tuning_instance_classif$result_learner_param_vals

# Print the best model parameters
print(best_params_classif)

# Manually apply class weights to the learner
# Define the class weights for ICI-resistant (higher weight) and ICI-sensitive
class_weights <- c("ICI-resistant" = 15, "ICI-sensitive" = 1)  # Adjust weights as needed

# Set the best parameters and class weights
learner_classif$param_set$values <- best_params_classif
learner_classif$param_set$values$class.weights <- class_weights  # Apply class weights here

# Train the model on the entire dataset using the best parameters and class weights
learner_classif$train(task_classif)

# Make predictions on the training set
predictions <- learner_classif$predict(task_classif)

# Compute the confusion matrix
conf_matrix <- predictions$confusion
print(conf_matrix)

