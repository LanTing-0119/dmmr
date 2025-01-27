library(maftools)  # for handling MAF objects
library(dplyr)     # for data manipulation
library(tidyr)     # for reshaping data
library(ggplot2)   # for plotting

dmmr_all = read.maf(maf = 'outputs/dmmr.maf.gz', clinicalData = 'outputs/dmmr_clinical.csv')

response_info = dmmr_all@clinical.data %>%
  select(Tumor_Sample_Barcode, ICI_response)

ICI_resistant_id = response_info %>% filter(ICI_response=="ICI-resistant") %>% pull(Tumor_Sample_Barcode) %>% as.character()
dmmr_resist = subsetMaf(maf = dmmr_all, tsb = ICI_resistant_id)
ICI_sensitive_id = response_info %>% filter(ICI_response=="ICI-sensitive") %>% pull(Tumor_Sample_Barcode) %>% as.character()
dmmr_sen = subsetMaf(maf = dmmr_all, tsb = ICI_sensitive_id)

# combine subgroups to a list for easy analysis(3 subgroups)
dmmr_list = list(total=dmmr_all, resistant=dmmr_resist, sensitive=dmmr_sen)

for (subgroup in names(dmmr_list)){
  
  # Create a directory for figures
  dir.create(paste0("figures/",subgroup), recursive = TRUE)
  dmmr = dmmr_list[[subgroup]]
  
  #### Part1: Visualization ####
  
  ##### 1.1 MAF summary ####
  maf_summary_plot = plotmafSummary(maf = dmmr, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  pdf(paste0("figures/", subgroup, "_01_01_maf_summary.pdf"))
  print(maf_summary_plot)  # Print plot to the PDF
  dev.off()  # Close the PDF device
  
  ##### 1.2 Oncoplots ####
  n = 10
  pdf(paste0("figures/", subgroup, "_01_02_oncoplot.pdf"), width = 10, height = 8)
  oncoplot(maf = dmmr, top = n)  # Oncoplot for the top 10 mutated genes
  dev.off()  # Close the PDF device
  
  ##### 1.3 Transition and Transversions ####
  dmmr_titv = titv(maf = dmmr, plot = FALSE, useSyn = TRUE)
  pdf(paste0("figures/", subgroup, "_01_03_titv.pdf"))
  plotTiTv(res = dmmr_titv)  # Plot Transition and Transversions
  dev.off()  # Close the PDF device
  
  ##### 1.6 Compare mutation load against TCGA cohorts ####
  dmmr_mutload = tcgaCompare(maf = dmmr, cohortName = subgroup, logscale = TRUE, capture_size = 50)
  pdf(paste0("figures/", subgroup, "_01_04_mutation_load_comparison.pdf"))
  print(dmmr_mutload)  # Print mutation load comparison
  dev.off()  # Close the PDF device
  
  #### Part2: Analysis ####
  
  ##### 2.1 Somatic Interaction ####
  pdf(paste0("figures/", subgroup, "_02_01_somatic_interaction.pdf"))
  somaticInteractions(maf = dmmr, top = 25, pvalue = c(0.05, 0.1))
  dev.off()  # Close the PDF device
  
  ##### 2.2 Detecting cancer driver ####
  dmmr_sig = oncodrive(maf = dmmr, AACol = 'aaChange', minMut = 5, pvalMethod = 'zscore')
  pdf(paste0("figures/", subgroup, "_02_02_cancer_driver_detection.pdf"))
  plotOncodrive(res = dmmr_sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
  dev.off()  # Close the PDF device
  
  ##### 2.3 Adding and summarizing pfam (protein family) domains ####
  pdf(paste0("figures/", subgroup, "_02_03_pfam_domains.pdf"))
  dmmr_pfam = pfamDomains(maf = dmmr, AACol = 'aaChange', top = 10)
  print(dmmr_pfam$proteinSummary[, 1:7, with = FALSE])  # Protein summary
  dev.off()  # Close the PDF device
  
  ##### 2.4 Survival analysis ####

  ###### 2.4.1 Differences between mutation in any gene ####
  genes = dmmr@data$Hugo_Symbol
  genes = 'ARID1A'
  for (gene in genes){
    pdf(paste0("figures/", subgroup, "_02_04_01_survival_analysis_", gene, ".pdf"))
    mafSurvival(maf = dmmr, genes = gene, time = 'PFS.month', Status = 'Status', isTCGA = FALSE)
    dev.off()  # Close the PDF device
  }
  
  ###### 2.4.2 Predict genesets associated with survival ####
  #Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
  pdf(paste0("figures/", subgroup, "_02_04_01_survival_analysis_", gene, ".pdf"))
  prog_geneset = survGroup(maf = dmmr, top = 100, geneSetSize = 2, time = "PFS.month", Status = "Status", verbose = FALSE)
  print(prog_geneset)
  prog_geneset$P_value %>% hist()
  mafSurvGroup(maf = dmmr, geneSet = c("PMS1", "ERCC3"), time = "PFS.month", Status = "Status")
  
  ##### 2.5 Drug-Gene Interactions ####
  dgi = drugInteractions(maf = dmmr, fontSize = 0.75)
  # get known/reported drugs to interact with DNMT3A
  genes = c("TP53","POLE")
  genes.dgi = drugInteractions(genes = genes, drugs = TRUE)
  genes.dgi[,.(genes, interaction_types, drug_name, drug_claim_name)]
  
  ##### 2.6 Oncogenic Signaling Pathways ####
  pws = pathways(maf = dmmr, plotType = 'treemap')
  plotPathways(maf = dmmr, pathlist = pws)
  
  ##### 2.7 Tumor heterogeneity and MATH scores ####
  # 待完善
  ###### 2.7.1 Heterogeneity in tumor samples ####
  library("mclust")
  
  ##### 2.8 mutational signatures ####
  ###### 2.8.1 APOBEC Enrichment estimation and Differences ####
  library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
  dmmr.tnm = trinucleotideMatrix(maf = dmmr, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
  plotApobecDiff(tnm = dmmr.tnm, maf = dmmr, pVal = 0.2)
  ###### 2.8.2 Signature analysis ####
  library('NMF')
  dmmr.sign = estimateSignatures(mat = dmmr.tnm, nTry = 6)
  plotCophenetic(res = dmmr.sign) # choose the best possible value
  dmmr.sig = extractSignatures(mat = dmmr.tnm, n = 3)
  # Compare detected signatures to COSMIC Legacy or SBS signature database.
  dmmr.og30.cosm = compareSignatures(nmfRes = dmmr.sig, sig_db = "legacy")
  # Compate against updated version3 60 signatures 
  dmmr.v3.cosm = compareSignatures(nmfRes = dmmr.sig, sig_db = "SBS")
  library('pheatmap')
  pheatmap::pheatmap(mat = dmmr.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
  maftools::plotSignatures(nmfRes = dmmr.sig, title_size = 1.2, sig_db = "SBS")
  # Visualize first signature for fancy 3D plot
  # 待完善
  library("barplot3d")
  sig1 = dmmr.sig$signatures[,1]
  barplot3d::legoplot3d(contextdata = sig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)

}
#### Part3: Comparing 2 cohorts ####

##### 3.1 Forest plots ####
sen_vs_resist <- mafCompare(m1 = dmmr_sen, m2 = dmmr_resist, m1Name = 'Sensitive', m2Name = 'Resistant', minMut = 5)
pdf(paste0("figures/", subgroup, "_09_forest_plot.pdf"))
forestPlot(mafCompareRes = sen_vs_resist, pVal = 0.05)
dev.off()  # Close the PDF device

##### 3.2 Co-onco plots ####
top = 20
genes = sen_vs_resist$results$Hugo_Symbol[1:top]
pdf(paste0("figures/", subgroup, "_10_co_onco_plot.pdf"), width = 10, height = 8)
coOncoplot(m1 = dmmr_sen, m2 = dmmr_resist, m1Name = 'Sensitive', m2Name = 'Resistant', genes = genes, removeNonMutated = TRUE)
dev.off()  # Close the PDF device

##### 3.3 Co-bar plots ####
pdf(paste0("figures/", subgroup, "_11_co_bar_plot.pdf"), width = 10, height = 8)
coBarplot(m1 = dmmr_sen, m2 = dmmr_resist, m1Name = "Sensitive", m2Name = "Resistant")
dev.off()  # Close the PDF device

##### 3.4 Lollipop plot-2 ####
gene = "TP53"
pdf(paste0("figures/", subgroup, "_12_lollipop_plot_2.pdf"))
lollipopPlot2(m1 = dmmr_sen, m2 = dmmr_resist, gene = gene, AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Sensitive", m2_name = "Resistant")
dev.off()  # Close the PDF device

##### 3.5 Clinical enrichment analysis ####
response_ce = clinicalEnrichment(maf = dmmr_all, clinicalFeature = 'ICI_response', minMut = 3, useCNV = FALSE, pathways = FALSE)
pdf(paste0("figures/", subgroup, "_13_clinical_enrichment.pdf"))
print(response_ce$groupwise_comparision$p_value %>% hist())
dev.off()  # Close the PDF device

