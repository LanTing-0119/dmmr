# Sample data as a character vector
data = read.table(file = 'data/data_external/scRNA/2023_Deng/sample_info.txt', sep = '\t')

# Read the data into a data frame
df <- data.frame(
  datasets = character(),
  labels = character(),
  stringsAsFactors = FALSE
)

n_rows = dim(data)[1]

# Process each line
for (i in 1:n_rows) {
  # Extract the dataset
  dataset <- data[i,1]
  
  # Split the second part by colon and comma
  details <- strsplit(data[i,2], ": ")[[1]][2]
  details <- strsplit(details, ", ")[[1]]
  
  # Extract the required information
  patient <- details[1]
  tissue <- details[2]
  condition <- details[3]
  treatment <- details[4]
  
  # Map the conditions and treatments to the specified labels
  tissue_label <- ifelse(condition == "normal", "N", "T")
  treatment_label <- ifelse(
    treatment == "untreated", "trt-",
    ifelse(treatment == "anti-PD-1", "trt+",
           ifelse(treatment == "Anti-PD-1+celecoxib", "trt++", treatment))
  )
  
  # Create the label
  label <- paste0(patient, tissue_label, "_", treatment_label)
  
  # Append to the data frame
  df <- rbind(df, data.frame(datasets = dataset, labels = label, stringsAsFactors = FALSE))
}

write.csv(df, file = 'data/data_external/scRNA/2023_Deng/all_samples.csv')
