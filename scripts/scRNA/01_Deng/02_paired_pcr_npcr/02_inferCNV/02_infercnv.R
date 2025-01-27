library(Seurat)

# dmmr_normal = readRDS('data/data_external/scRNA/2023_Deng/03_dmmr_merged_paired_pcr_npcr_normal.RDS')
# dmmr_Epi = readRDS('outputs/Deng/dmmr_epi.RDS')
samples = c('30','31')

for (sample in samples){
  
  # 创建逻辑向量 & 使用逻辑向量提取亚组
  cells <- grepl(paste0('^P',sample), dmmr_Epi@meta.data$orig.ident)
  subset <- subset(dmmr_Epi, cells = colnames(dmmr_Epi)[cells])
  
  # P30T_trt P30T_untrt   P31T_trt P31T_untrt 
  # 1094       3131       7390       2291 
  
  #### 1. Extract raw counts matrix
  # 提取正常组织的原始计数矩阵
  JoinLayers(dmmr_normal)
  normal_counts <- GetAssayData(dmmr_normal, assay = "RNA", layer = paste0('counts.P',sample,'N'))
  
  # 提取肿瘤组织的原始计数矩阵
  tumor_counts <- GetAssayData(subset, assay = "RNA", layer = "counts")
  
  # 检查基因是否一致
  # 确保肿瘤矩阵的基因顺序与正常矩阵一致
  genes = intersect(rownames(tumor_counts), rownames(normal_counts)) %>% sort()
  tumor_counts <- tumor_counts[genes,]
  normal_counts = normal_counts[genes,]
  
  # 合并矩阵
  combined_counts <- cbind(normal_counts, tumor_counts)
  
  
  
  # 提取肿瘤样本的簇信息
  tumor_clusters <- subset@meta.data$pca_clusters  # 假设聚类结果存储在 seurat_clusters 列
  tumor_cell_ids <- colnames(subset)  # 肿瘤样本的细胞名称
  
  # 创建细胞注释
  cell_annotations <- data.frame(
    cell_id = colnames(combined_counts),
    cell_type = ifelse(colnames(combined_counts) %in% colnames(normal_counts), "normal", "tumor")
  )
  
  # 为肿瘤细胞添加簇信息
  cell_annotations$cell_type[cell_annotations$cell_type == "tumor"] <- paste0("tumor_", tumor_clusters)
  write.table(cell_annotations, file = 'data/infercnv/cell_annotations.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # infercnv
  # create the infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=combined_counts,
                                      annotations_file='data/infercnv/cell_annotations.txt',
                                      delim="\t",
                                      gene_order_file="data/infercnv/hg38_gencode_v27.txt",
                                      ref_group_names=c("normal"))
  
  # perform infercnv operations to reveal cnv signal
  output_dir = paste0('outputs/Deng/infercnv/',sample)
  dir.create(output_dir, recursive = T)
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=output_dir,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T,
                               HMM=T
  )
}