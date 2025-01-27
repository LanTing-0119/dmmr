library(maftools)
library(data.table)
library(readxl)


#### 01. Import Data ####
annovar.laml <- annovarToMaf(annovar = 'data/annovar_merge.vcf',
                             refBuild = 'hg38',
                             tsbCol = 'Tumor_Sample_Barcode',
                             table = 'refGene',
                             MAFobj = T)
# Transform the object to a maf
saveRDS(annovar.laml,file = 'data/merged_maf.RDS')

maf_data <- readRDS(file = 'data/merged_maf.RDS')


# Clean the table
  
  # Step 1: Convert the factor to a character vector
  maf_data@data[["Tumor_Sample_Barcode"]] <- as.character(maf_data@data[["Tumor_Sample_Barcode"]])
  
  # Step 2.1: Remove the part between '-' and '_', inclusive
  maf_data@data[["Tumor_Sample_Barcode"]] <- sub("-[^_]+_", "_", maf_data@data[["Tumor_Sample_Barcode"]])
  # Step 2.2: Remove the suffix '_HUM'
  maf_data@data[["Tumor_Sample_Barcode"]] <- sub("_HUM$", "", maf_data@data[["Tumor_Sample_Barcode"]])
  
  # Step 3: Optionally, convert back to factor (if you need it to be a factor)
  maf_data@data[["Tumor_Sample_Barcode"]] <- factor(maf_data@data[["Tumor_Sample_Barcode"]])
  
  # Check the result
  head(maf_data@data[["Tumor_Sample_Barcode"]])

# Export the Maf Object  
maf_df = maf_data@data # Extract the MAF Data (from a`MAF`object)
output_file = 'outputs/dmmr.maf.gz' # Path to save the compressed MAF file
write.table(maf_df, gzfile(output_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# clinical_id = maf_data@clinical.data %>%
#   separate(Tumor_Sample_Barcode, into = c("First_Part", "Second_Part"), sep = "_", remove = TRUE) %>%
#   pull(First_Part) %>%
#   sub("-\\d+$", "", .) %>%
#   as.data.table() # 需要转变为table才能合并
# 
# maf_data@clinical.data = clinical_id
# 
# saveRDS(maf_data,file = 'outputs/dmmr_maf.RDS')


# Import the Excel file
maf_clinical <- read_excel("data/merged_output.xlsx")
colnames(maf_clinical)[colnames(maf_clinical) == "分子编号"] <- "Tumor_Sample_Barcode"
write.csv(maf_clinical, file = 'outputs/dmmr_clinical.csv', row.names = FALSE)



