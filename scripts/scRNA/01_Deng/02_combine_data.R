
library(Seurat)
datasets=c("GSM6213956","GSM6213957","GSM6213958","GSM6213959","GSM6213960","GSM6213961","GSM6213962","GSM6213963","GSM6213964","GSM6213965","GSM6213966","GSM6213967","GSM6213968","GSM6213969","GSM6213970","GSM6213971","GSM6213972","GSM6213973","GSM6213974","GSM6213975","GSM6213976","GSM6213977","GSM6213978","GSM6213979","GSM6213980","GSM6213981","GSM6213982","GSM6213983","GSM6213984","GSM6213985","GSM6213986","GSM6213987","GSM6213988","GSM6213989","GSM6213990","GSM6213991","GSM6213992","GSM6213993","GSM6213994","GSM6213995")
folder_input='/Users/lanting/Project/dmmr/data/data_external/scRNA/2023_Deng/01_clean_data'

# Create an empty list to store the Seurat objects
seurat_objects <- list()

# Loop over each dataset name in the 'datasets' list
for (dataset in datasets) {
  print(paste0(which(datasets==dataset),'/', length(datasets)))
  # Read in the data using Read10X
  data_path <- paste0(folder_input,'/', dataset, '/')
  dataset_data <- Read10X(data.dir = data_path)
  
  # Create the Seurat object
  seurat_object <- CreateSeuratObject(counts = dataset_data, project = dataset, min.cells = 3, min.features = 200)
  
  # Store the Seurat object in the list with the dataset name as the key
  seurat_objects[[dataset]] <- seurat_object
}

dmmr=merge(seurat_objects[[1]], y=c(seurat_objects[2:length(seurat_objects)]), project='merged_dmmr')
saveRDS(dmmr, file = 'data/data_external/scRNA/2023_Deng/02_dmmr_merged.RDS')
rm(seurat_objects)
# Now you can access the Seurat object for a specific dataset
# Example: Access Seurat object for GSM6213956
seurat_objects[["GSM6213956"]]
