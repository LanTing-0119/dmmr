
library(Seurat)

df = read.csv("data/data_external/scRNA/2023_Deng/paired_pcr_npcr.csv")
folder_input='/Users/lanting/Project/dmmr/data/data_external/scRNA/2023_Deng/01_clean_data'
list.dirs(folder_input)

# Create an empty list to store the Seurat objects
seurat_objects <- list()

# import data and createseurat object
for (dataset in df$datasets) {
  print(paste0(which(datasets==dataset),'/', length(datasets)))
  # Read in the data using Read10X
  label = df$labels[df$datasets == dataset]
  data_path <- paste0(folder_input,'/', dataset, '/')
  dataset_data <- Read10X(data.dir = data_path)
  
  # Create the Seurat object
  seurat_object <- CreateSeuratObject(counts = dataset_data, project = label, min.cells = 3, min.features = 200) # project 这个参数会影响orig.ident选项
  
  # Store the Seurat object in the list with the dataset name as the key
  seurat_objects[[label]] <- seurat_object
}

# quality control
# quality control
seurat_objects = lapply(X = seurat_objects, FUN = function(x){
  x = subset(x, 
             subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & 
             nCount_RNA > 400 & nCount_RNA < 25000
             )
})

dmmr=merge(seurat_objects[[1]], y=seurat_objects[-1], project='merged_dmmr') # add.cell.ids 会在所有细胞名称前添加对应所属组别的label(有时候Barcode可能会重复: 不用添加add.cell.ids, 会强制自动转换)
saveRDS(dmmr, file = 'data/data_external/scRNA/2023_Deng/03_dmmr_merged_paired_pcr_npcr.RDS')


# rm(seurat_objects)
# Now you can access the Seurat object for a specific dataset
# Example: Access Seurat object for GSM6213956
# seurat_objects[["GSM6213956"]]
