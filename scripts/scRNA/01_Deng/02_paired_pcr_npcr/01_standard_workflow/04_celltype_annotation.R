
library(Seurat)

dmmr = readRDS('outputs/Deng/dmmr_merged_paired_pcr_npcr.RDS')

dmmr <- FindNeighbors(dmmr, reduction = "integrated.cca", dims = 1:30)
dmmr <- FindClusters(dmmr, resolution = 0.2, cluster.name = "cca_clusters")

dmmr <- RunUMAP(dmmr, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

p1 <- DimPlot(
  dmmr,
  reduction = "umap.cca",
  group.by = c("Method", "predicted.celltype.l0.2", "cca_clusters"),
  combine = FALSE, label.size = 2
)

DimPlot(dmmr, reduction = "umap.cca", split.by = "orig.ident")

dmmr = PrepSCTFindMarkers(dmmr)
dmmr.markers <- FindAllMarkers(dmmr, only.pos = TRUE)

n=30
dmmr.markers_top_n = dmmr.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>% 
  dplyr::arrange(desc(avg_log2FC)) %>%
  top_n(n=n, wt = avg_log2FC)

formatted_genes <- dmmr.markers_top_n %>%
  group_by(cluster) %>%                                      # Group by 'cluster'
  summarise(genes = paste(gene, collapse = ",")) %>%         # Concatenate genes with commas
  arrange(as.numeric(as.character(cluster))) %>%              # Sort clusters numerically
  mutate(formatted = paste0(cluster, ":\n", genes)) %>%      # Create the desired format
  pull(formatted)                                            # Extract the formatted strings

writeLines(formatted_genes, 'outputs/Deng/cluster_top_n.txt')


# this markers are from paper: Remodeling of the immune and stromal cell compartment by PD-1 blockade in mismatch repair- deficient colorectal cancer
markers = c('CD3D', 'CD3E', 'TRAC', 'TRBC1') # T cells 3,5
markers = c('CD79A', 'CD79B', 'MS4A1', 'TNFRSF17', 'MZB1') # B cells 7,8
markers = c('CD14', 'CD68') # Myeloid cells 9
markers = c('EPCAM', 'CD24') # Epithelial cells 0,1,2,6
markers = c('COL1A2', 'COL3A1', 'MYH11', 'ACTA2') # Fibroblast 10
markers = c('VWF', 'PECAM1') # Endothelial cells 4

markers = c('CD3D', 'CD3E', 'TRAC', 'TRBC1', 'CD79A', 'CD79B', 'MS4A1', 'TNFRSF17', 'MZB1', 'CD14', 'CD68', 'EPCAM', 'CD24', 'COL1A2', 'COL3A1', 'MYH11', 'ACTA2', 'VWF', 'PECAM1')

# use dotplot to help annotate the celltype
DotPlot(dmmr, features = markers)  + RotatedAxis()
# VlnPlot(dmmr, features = markers)
# FeaturePlot(dmmr, features = markers)

# create a datafrane to record the relationship
celltype = data.frame(ClusterID = 0:10,
                      Celltype = 'Unknown')
celltype[celltype$ClusterID %in% c(3,5), 2] = 'Tcell'
celltype[celltype$ClusterID %in% c(7,8), 2] = 'Bcell'
celltype[celltype$ClusterID %in% c(9), 2] = 'Myeloid'
celltype[celltype$ClusterID %in% c(0,1,2,6), 2] = 'Epithelial'
celltype[celltype$ClusterID %in% c(10), 2] = 'Fibroblast'
celltype[celltype$ClusterID %in% c(4), 2] = 'Endothelial'

# update the cluster.id in seurat
new.cluster.ids <- celltype$Celltype
names(new.cluster.ids) <- levels(dmmr)
dmmr <- RenameIdents(dmmr, new.cluster.ids)
DimPlot(dmmr, reduction = "umap.cca", label = TRUE, pt.size = 0.5) + NoLegend()

