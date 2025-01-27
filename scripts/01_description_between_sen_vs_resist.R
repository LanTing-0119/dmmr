library(maftools)  # for handling MAF objects
library(dplyr)     # for data manipulation
library(tidyr)     # for reshaping data
library(ggplot2)   # for plotting

dmmr_all = read.maf(maf = 'data/vep_merge_oncokb.maf', clinicalData = 'outputs/dmmr_clinical_hg38.csv')
dmmr_all_filtered = subsetMaf(maf=dmmr_all, query = "ONCOGENIC %in% c('Inconclusive','Likely Neutral','Likely Oncogenic','Oncogenic','Resistance')")
response_info = dmmr_all@clinical.data %>%
  select(Tumor_Sample_Barcode, ICI_response) %>%
  filter(!is.na(ICI_response))

# write.maf

ICI_resistant_id = response_info %>% filter(ICI_response=="ICI-resistant") %>% pull(Tumor_Sample_Barcode) %>% as.character()
dmmr_resist = subsetMaf(maf = dmmr_all, tsb = ICI_resistant_id)
ICI_sensitive_id = response_info %>% filter(ICI_response=="ICI-sensitive") %>% pull(Tumor_Sample_Barcode) %>% as.character()
dmmr_sen = subsetMaf(maf = dmmr_all, tsb = ICI_sensitive_id)

# combine subgroups to a list for easy analysis(3 subgroups)
dmmr_list = list(total=dmmr_all, resistant=dmmr_resist, sensitive=dmmr_sen)
dmmr=dmmr_all_filtered
for (subgroup in names(dmmr_list)){
  
  dir.create(paste0("figures/",subgroup), recursive = TRUE)
  dmmr = dmmr_list[[subgroup]]
  
  #### Part1: Visualization ####
  
  ##### 1.1 MAF summary ####
  maf_summary_plot=plotmafSummary(maf = dmmr, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  ##### 1.2 Oncoplots ####
  n=40
  oncoplot(maf = dmmr, top = n)
  # oncoplot(maf = dmmr_resist, top = 801)
  
  ##### 1.3 Transition and Transversions ####
  dmmr.titv = titv(maf = dmmr, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = dmmr.titv)
  
  # ## --------------------------
  # ## --------------------------
  # 
  # ##### 1.4 Lollipop plots for amino acid changes ####
  # # 待完善
  # genes = dmmr@data$Hugo_Symbol
  # genes = 'TP53'
  # for (gene in genes){
  #   tryCatch(
  #     {
  #       lollipopPlot(
  #       maf = dmmr,
  #       gene = gene,
  #       AACol = 'aaChange', # 还要修改
  #       showMutationRate = TRUE,
  #       labelPos = 882
  #     )
  #       dev.off},
  #     error = function(e){
  #       print(gene)
  #       }
  #   )
  # }
  # 
  # ##### 1.5 Rainfall plots ####
  # # 待完善
  # patients = dmmr@clinical.data$Tumor_Sample_Barcode %>% as.character()
  # patients =  dmmr@clinical.data$Tumor_Sample_Barcode[1] %>% as.character()
  # for (patient in patients){
  #   print(paste0(which(patient==patients),' / ',length(patients))) 
  #   tryCatch({
  #     dmmr.patient = subsetMaf(maf = dmmr, tsb = patient)
  #     rainfallPlot(maf = dmmr.patient, detectChangePoints = TRUE, pointSize = 0.4)
  #     dev.off},
  #     error = function(e){
  #     print(patient)
  #     }
  #   )
  # }
  # 
  # ## --------------------------
  # ## --------------------------
  
  ##### 1.6 Compare mutation load against TCGA cohorts ####
  dmmr.mutload = tcgaCompare(maf = dmmr, cohortName = subgroup, logscale = TRUE, capture_size = 50)
  dmmr.mutload$
  dev.off
  
  ##### 1.7 Plotting VAF ####
  # lack of vaf information
  # plotVaf(maf = dmmr, vafCol = 't_vaf')
  
  
  #### Part2: Analysis ####
  
  ##### 2.1 Somatic Interaction ####
  somaticInteractions(maf = dmmr, top = 25, pvalue = c(0.05, 0.1))
  
  ##### 2.2 Detecting cancer driver ####
  dmmr.sig = oncodrive(maf = dmmr, AACol = 'aaChange', minMut = 5, pvalMethod = 'zscore')
  plotOncodrive(res = dmmr.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
  
  ##### 2.3 Adding and summarizing pfam(protein family) domains ####
  dmmr.pfam = pfamDomains(maf = dmmr, AACol = 'aaChange', top = 10)
  dmmr.pfam$proteinSummary[,1:7, with = FALSE] # Protein summary
  dmmr.pfam$domainSummary[,1:3, with = FALSE] # Domain summary
  
  ##### 2.4 Survival analysis ####
  
  ###### 2.4.1 Differences between mutation in any gene ####
  genes = dmmr@data$Hugo_Symbol
  genes = 'KRAS'
  for (gene in genes){
  mafSurvival(maf = dmmr, genes = gene, time = 'PFS.month', Status = 'Status', isTCGA = FALSE)
  }
  
  ###### 2.4.2 Predict genesets associated with survival ####
  #Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
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
sen.vs.resist <- mafCompare(m1 = dmmr_sen, m2 = dmmr_resist, m1Name = 'Sensitive', m2Name = 'Resistant', minMut = 5)
print(sen.vs.resist)
forestPlot(mafCompareRes = sen.vs.resist, pVal = 0.05)

  # plot common genes with pCR after neoadjuvant
  .genes1 = sen.vs.resist$results %>% filter(pval<0.05) %>% pull(Hugo_Symbol)
  .sen.vs.resist=sen.vs.resist
  genes_common=intersect(.genes1,.genes2)
  .sen.vs.resist$results = .sen.vs.resist$results %>% filter(Hugo_Symbol %in% genes_common)
  forestPlot(mafCompareRes = .sen.vs.resist)

##### 3.2 Co-onco plots ####
top=20
genes=sen.vs.resist$results$Hugo_Symbol[1:top]
coOncoplot(m1 = dmmr_sen, m2 = dmmr_resist, m1Name = 'Sensitive', m2Name = 'Resistant', genes = genes, removeNonMutated = TRUE)

##### 3.3 Co-bar plots ####
coBarplot(m1 = dmmr_sen, m2 = dmmr_resist, m1Name = "Sensitive", m2Name = "Resistant")

##### 3.4 Lollipop plot-2 ####
gene="B2M"
lollipopPlot2(m1 = dmmr_sen, m2 = dmmr_resist, gene = gene, AACol1 = "aaChange", AACol2 = "aaChange", m1_name = "Sensitive", m2_name = "Resistant")

##### 3.5 Clinical enrichment analysis ####
response.ce = clinicalEnrichment(maf = dmmr_all, clinicalFeature = 'ICI_response', minMut = 3, useCNV = FALSE, pathways = FALSE)
response.ce$groupwise_comparision$p_value %>% hist()
response.ce$groupwise_comparision[p_value < 0.1]

##### 3.6 Oncogenic Signaling Pathways ####
pws = pathways(maf = dmmr, plotType = 'treemap')

# This will produce the figure with correlation of top 25 mutated genes
somaticInteractions(
  maf = laml,
  top = 25,  # Number of top mutated genes to include
  pvalue = c(0.05, 0.1) # Significance levels to display
)



