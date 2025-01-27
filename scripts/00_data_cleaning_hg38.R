library(maftools)
library(data.table)
library(readxl)


#### 01. Import Data ####
maf_data <- read.maf(maf = 'data/vep_merge.maf')


# Import the Excel file
maf_clinical <- read_excel("data/merged_output.xlsx")
colnames(maf_clinical)[colnames(maf_clinical) == "分子编号"] <- "Tumor_Sample_Barcode"
maf_clinical$Tumor_Sample_Barcode <- paste0(maf_clinical$Tumor_Sample_Barcode, "_HUM_C")
write.csv(maf_clinical, file = 'outputs/dmmr_clinical_hg38.csv', row.names = FALSE)


# get df for oncokb annotation
dmmr_all = read.maf(maf = 'data/vep_merge.maf', clinicalData = 'outputs/dmmr_clinical_hg38.csv')
# oncokb_clinical_df <- data.frame(Column1 = dmmr_all@clinical.data$Tumor_Sample_Barcode, Column2 = rep('COAD', length(dmmr_all@clinical.data$Tumor_Sample_Barcode)))
oncokb_clinical_df <- data.frame(Column1 = dmmr_all@clinical.data$Tumor_Sample_Barcode, Column2 = rep('COADREAD', length(dmmr_all@clinical.data$Tumor_Sample_Barcode)))
colnames(oncokb_clinical_df) = c('Tumor_Sample_Barcode','ONCOTREE_CODE') 
write.csv(oncokb_clinical_df, file = 'outputs/oncokb_clinical.csv', row.names = FALSE)

# run oncokb annotation
# /scripts/01_oncokb_annotation.sh

# update the maf file: filter out mutations which is unknown
