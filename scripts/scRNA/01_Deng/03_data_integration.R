rm(list = ls())

dmmr=readRDS('data/data_external/scRNA/2023_Deng/03_dmmr_merged_paired_pcr_npcr.RDS')
library(Seurat)
library(BPCells)
options(future.globals.maxSize = 3e+11)
dmmr <- SCTransform(dmmr)
dmmr <- RunPCA(dmmr, npcs = 30, verbose = F)
dmmr <- IntegrateLayers(
  object = dmmr,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
dmmr <- FindNeighbors(dmmr, dims = 1:30, reduction = "integrated.dr")
dmmr <- FindClusters(dmmr, resolution = 2)