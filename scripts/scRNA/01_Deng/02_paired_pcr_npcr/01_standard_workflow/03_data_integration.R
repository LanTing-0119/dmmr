rm(list = ls())
library(Seurat)
library(BPCells)
options(future.globals.maxSize = 3e+11)

dmmr=readRDS('data/data_external/scRNA/2023_Deng/03_dmmr_merged_paired_pcr_npcr.RDS')

dmmr <- SCTransform(dmmr)
dmmr <- RunPCA(dmmr, npcs = 30, verbose = F)

dmmr <- IntegrateLayers(
  object = dmmr,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

# saveRDS(dmmr, file = 'outputs/Deng/dmmr_merged_paired_pcr_npcr.RDS')
