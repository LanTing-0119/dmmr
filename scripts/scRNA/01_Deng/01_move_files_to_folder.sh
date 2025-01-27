#!/bin/bash

folder_input=/Users/lanting/Project/dmmr/data/data_external/scRNA/2023_Deng/GSE205506_RAW
folder_output=/Users/lanting/Project/dmmr/data/data_external/scRNA/2023_Deng/01_clean_data
# Loop through each file in the directory
for file in ${folder_input}/*; do
  
  echo ${file} # /Users/lanting/Project/dmmr/data/data_external/scRNA/2023_Deng/GSE205506_RAW/GSM6213992_Colon_MGI_321_v3_barcodes.tsv.gz
  
  # Extract the prefix before the first underscore
  prefix=$(basename $file | cut -d'_' -f1) # GSM6213992
  echo ${prefix}
  echo ${folder_output}/"$prefix"

  # Create a directory named after the prefix
  # mkdir -p ${folder_output}/"$prefix"

  # Get the base filename (to make sure no directories are included)
  basefile=$(basename "$file") # GSE205506_RAW/GSM6213992_Colon_MGI_321_v3_barcodes.tsv.gz
  echo ${basefile}
  # Create symbolic links for the relevant files
  case "$basefile" in
    *barcodes.tsv.gz)
      ln -s "$file" ${folder_output}/$prefix/barcodes.tsv.gz
      ;;
    *features.tsv.gz)
      ln -s "$file" ${folder_output}/$prefix/features.tsv.gz
      ;;
    *matrix.mtx.gz)
      ln -s "$file" ${folder_output}/$prefix/matrix.mtx.gz
      ;;
  esac
done
