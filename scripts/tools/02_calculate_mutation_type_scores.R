# get mutation score of a gene of a sample. 
# input: maf@data
# output: a dataframe with 3 columns: Tumor_Sample_Barcode, Hugo_Symbol, Total_Mutation_Score
# Example usage:
# Assuming `dmmr@data` is your MAF data
# mutation_scores <- calculate_mutation_scores(dmmr@data)
# head(mutation_scores)

# Load necessary libraries
library(dplyr)
library(tidyr)

# Define the function
calculate_mutation_scores <- function(maf_data) {
  
  # Define the mapping for Variant_Classification to scores
  mutation_mapping <- c(
    "Missense_Mutation" = 1,
    "Frame_Shift_Del" = 2,
    "Splice_Site" = 4,
    "In_Frame_Del" = 8,
    "Frame_Shift_Ins" = 16,
    "Nonsense_Mutation" = 32,
    "Nonstop_Mutation" = 64,
    "In_Frame_Ins" = 128,
    "Translation_Start_Site" = 256
  )
  
  # Apply the mapping and create the Variant_Classification_Score column
  maf_data <- maf_data %>%
    mutate(
      Variant_Classification_Score = recode(
        Variant_Classification, 
        !!!mutation_mapping
      )
    )
  
  # Calculate the total mutation score per gene and tumor sample
  mutation_scores <- maf_data %>%
    group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    summarise(Total_Mutation_Score = sum(Variant_Classification_Score, na.rm = TRUE)) %>%
    ungroup()  # Remove grouping
  
  # Return the reshaped data
  return(mutation_scores)
}


