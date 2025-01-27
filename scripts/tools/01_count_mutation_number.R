dmmr_all = read.maf(maf = 'data/vep_merge.maf', clinicalData = 'outputs/dmmr_clinical_hg38.csv')

# count mutation gene number
n=0
dmmr=dmmr_resist
count=0
gene = 'SMARCA4'
for (item in dmmr_all@data[["Gene.refGene"]]){
  n=n+1
  if (item == gene){
    count=count+1
    print(item)
  }
}
# n
count


# count mutation gene with specific mutation site
n=0
dmmr=dmmr_all
count=0
gene = 'KRAS'
site = 'p.A146T' # protein change
for (item in dmmr_all@data[["Hugo_Symbol"]]){
  n=n+1
  if (item == gene){
    print(item)
    print(dmmr@data[["HGVSp_Short"]][n])
    if (dmmr@data[["HGVSp_Short"]][n] == site){
      count=count+1
      # print(dmmr_all@data[["Tumor_Sample_Barcode"]][n])
    }
  }
}
# n
count





