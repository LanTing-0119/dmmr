# load the CytoTRACE 2 package
library(CytoTRACE2) 

# download the .rds file (this will download the file to your working directory)
# download.file("https://drive.google.com/uc?export=download&id=1TYdQsMoDIJjoeuiTD5EO_kZgNJUyfRY2", "Pancreas_10x_downsampled.rds")

# load rds
# data <- readRDS("Pancreas_10x_downsampled.rds")
dmmr_Epi = readRDS('outputs/Deng/dmmr_epi.RDS')

samples = c('30','31')

for (sample in samples){
  
  # 创建逻辑向量 & 使用逻辑向量提取亚组
  cells <- grepl(paste0('^P',sample), dmmr_Epi@meta.data$orig.ident)
  subset <- subset(dmmr_Epi, cells = colnames(dmmr_Epi)[cells])
  
  expression_data <- GetAssayData(subset, assay = "SCT", layer = "counts")
  
  # running CytoTRACE 2 main function - cytotrace2 - with default parameters
  cytotrace2_result <- cytotrace2(expression_data, ncores = 6, species = 'human')
  saveRDS(cytotrace2_result, file = paste0('outputs/Deng/cytotrace2_result/',sample,'.RDS'))
  
  # extract annotation data(need a dataframe output)
  subset@meta.data$clusters.group=paste0(subset@meta.data$pca_clusters, '_', subset@meta.data$orig.ident)
  annotation <- subset@meta.data['clusters.group']
  
  # generate prediction and phenotype association plots with plotData function
  plots <- plotData(cytotrace2_result = cytotrace2_result, 
                    annotation = annotation,
                    expression_data = expression_data
  )
}



