
library(Seurat)
library(dplyr)
options(future.globals.maxSize = 3e+11)


Epi = subset(dmmr, idents = 'Epithelial')

dmmr_Epi = CreateSeuratObject(counts = GetAssayData(Epi, assay = 'SCT', layer = 'counts'), 
                              meta.data = Epi@meta.data)

dmmr_Epi = SCTransform(dmmr_Epi) %>% 
  RunPCA(npcs = 30, verbose = F) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2, cluster.name = "pca_clusters") %>%
  RunUMAP(reduction = "pca", dims = 1:30, reduction.name = "umap")

DimPlot(dmmr_Epi, reduction = 'umap', label = T, raster = FALSE)
DimPlot(dmmr_Epi, reduction = 'umap', split.by = 'orig.ident', label = T, raster = FALSE)

# saveRDS(dmmr_Epi, file = 'outputs/Deng/dmmr_epi.RDS')
