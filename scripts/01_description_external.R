#### Step 1: Load the Required Libraries ####
library(maftools)  # for handling MAF objects
library(dplyr)     # for data manipulation
library(tidyr)     # for reshaping data
library(ggplot2)   # for plotting

#### Step 2: Load the MAF File and Clinical Data ####
maf_data = 'data/data_external/healthy_bowel/SYSUCC_CRC_Mutation_withSilent.maf'
changkang_all = read.maf(maf = maf_data)

changkang

# combine subgroups to a list for easy analysis(3 subgroups)
changkang_list = list(total=changkang_all)
# changkang_list = list(total=changkang_all, resistant=changkang_resist, sensitive=changkang_sen)

for (subgroup in names(changkang_list)){
  
  dir.create(paste0("figures/",subgroup), recursive = TRUE)
  changkang = changkang_list[[subgroup]]
  
  #### Part1: Visualization ####
  
  ##### 1.1 MAF summary ####
  maf_summary_plot=plotmafSummary(maf = changkang, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  ##### 1.2 Oncoplots ####
  n=10
  oncoplot(maf = changkang, top = n)
  # oncoplot(maf = changkang_resist, top = 801)
  
  ##### 1.3 Transition and Transversions ####
  changkang.titv = titv(maf = changkang, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = changkang.titv)
  
  # ## --------------------------
  # ## --------------------------
  # 
  # ##### 1.4 Lollipop plots for amino acid changes ####
  # # 待完善
  # genes = changkang@data$Hugo_Symbol
  # genes = 'TP53'
  # for (gene in genes){
  #   tryCatch(
  #     {
  #       lollipopPlot(
  #       maf = changkang,
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
  # patients = changkang@clinical.data$Tumor_Sample_Barcode %>% as.character()
  # patients =  changkang@clinical.data$Tumor_Sample_Barcode[1] %>% as.character()
  # for (patient in patients){
  #   print(paste0(which(patient==patients),' / ',length(patients))) 
  #   tryCatch({
  #     changkang.patient = subsetMaf(maf = changkang, tsb = patient)
  #     rainfallPlot(maf = changkang.patient, detectChangePoints = TRUE, pointSize = 0.4)
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
  changkang.mutload = tcgaCompare(maf = changkang, cohortName = subgroup, logscale = TRUE, capture_size = 50)
  changkang.mutload$
    dev.off
  
  ##### 1.7 Plotting VAF ####
  # lack of vaf information
  # plotVaf(maf = changkang, vafCol = 't_vaf')
  
  
  #### Part2: Analysis ####
  
  ##### 2.1 Somatic Interaction ####
  somaticInteractions(maf = changkang, top = 25, pvalue = c(0.05, 0.1))
  
  ##### 2.2 Detecting cancer driver ####
  changkang.sig = oncodrive(maf = changkang, AACol = 'aaChange', minMut = 5, pvalMethod = 'zscore')
  plotOncodrive(res = changkang.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
  
  ##### 2.3 Adding and summarizing pfam(protein family) domains ####
  changkang.pfam = pfamDomains(maf = changkang, AACol = 'aaChange', top = 10)
  changkang.pfam$proteinSummary[,1:7, with = FALSE] # Protein summary
  changkang.pfam$domainSummary[,1:3, with = FALSE] # Domain summary
  
  ##### 2.4 Survival analysis ####
  
  ###### 2.4.1 Differences between mutation in any gene ####
  genes = changkang@data$Hugo_Symbol
  genes = 'KRAS'
  for (gene in genes){
    mafSurvival(maf = changkang, genes = gene, time = 'PFS.month', Status = 'Status', isTCGA = FALSE)
  }
  
  ###### 2.4.2 Predict genesets associated with survival ####
  #Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
  prog_geneset = survGroup(maf = changkang, top = 100, geneSetSize = 2, time = "PFS.month", Status = "Status", verbose = FALSE)
  print(prog_geneset)
  prog_geneset$P_value %>% hist()
  mafSurvGroup(maf = changkang, geneSet = c("PMS1", "ERCC3"), time = "PFS.month", Status = "Status")
  
  ##### 2.5 Drug-Gene Interactions ####
  dgi = drugInteractions(maf = changkang, fontSize = 0.75)
  # get known/reported drugs to interact with DNMT3A
  genes = c("TP53","POLE")
  genes.dgi = drugInteractions(genes = genes, drugs = TRUE)
  genes.dgi[,.(genes, interaction_types, drug_name, drug_claim_name)]
  
  ##### 2.6 Oncogenic Signaling Pathways ####
  pws = pathways(maf = changkang, plotType = 'treemap')
  plotPathways(maf = changkang, pathlist = pws)
  
  ##### 2.7 Tumor heterogeneity and MATH scores ####
  # 待完善
  ###### 2.7.1 Heterogeneity in tumor samples ####
  library("mclust")
  
  ##### 2.8 mutational signatures ####
  ###### 2.8.1 APOBEC Enrichment estimation and Differences ####
  library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
  changkang.tnm = trinucleotideMatrix(maf = changkang, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
  plotApobecDiff(tnm = changkang.tnm, maf = changkang, pVal = 0.2)
  ###### 2.8.2 Signature analysis ####
  library('NMF')
  changkang.sign = estimateSignatures(mat = changkang.tnm, nTry = 6)
  plotCophenetic(res = changkang.sign) # choose the best possible value
  changkang.sig = extractSignatures(mat = changkang.tnm, n = 3)
  # Compare detected signatures to COSMIC Legacy or SBS signature database.
  changkang.og30.cosm = compareSignatures(nmfRes = changkang.sig, sig_db = "legacy")
  # Compate against updated version3 60 signatures 
  changkang.v3.cosm = compareSignatures(nmfRes = changkang.sig, sig_db = "SBS")
  library('pheatmap')
  pheatmap::pheatmap(mat = changkang.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
  maftools::plotSignatures(nmfRes = changkang.sig, title_size = 1.2, sig_db = "SBS")
  # Visualize first signature for fancy 3D plot
  # 待完善
  library("barplot3d")
  sig1 = changkang.sig$signatures[,1]
  barplot3d::legoplot3d(contextdata = sig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)
}

#### Part3: Comparing 2 cohorts ####
##### 3.1 Forest plots ####
pCR.vs.npCR <- mafCompare(m1 = changkang_pCR, m2 = changkang_npCR, m1Name = 'pCR', m2Name = 'npCR', minMut = 5)
print(pCR.vs.npCR)
forestPlot(mafCompareRes = pCR.vs.npCR, pVal = 0.05)

# plot common genes with primary resistant
.genes2 = pCR.vs.npCR$results %>% filter(pval<0.05) %>% pull(Hugo_Symbol)
.pCR.vs.npCR=pCR.vs.npCR
genes_common=intersect(.genes1,.genes2)
.pCR.vs.npCR$results = .pCR.vs.npCR$results %>% filter(Hugo_Symbol %in% genes_common)
forestPlot(mafCompareRes = .pCR.vs.npCR)

##### 3.2 Co-onco plots ####
top=20
genes=sen.vs.resist$results$Hugo_Symbol[1:top]
coOncoplot(m1 = changkang_pCR, m2 = changkang_npCR, m1Name = 'pCR', m2Name = 'npCR', genes = genes, removeNonMutated = TRUE)

##### 3.3 Co-bar plots ####
coBarplot(m1 = changkang_pCR, m2 = changkang_npCR, m1Name = 'pCR', m2Name = 'npCR')

##### 3.4 Lollipop plot-2 ####
gene="SMARCA4"
lollipopPlot2(m1 = changkang_pCR, m2 = changkang_npCR, gene = gene, AACol1 = "aaChange", AACol2 = "aaChange", m1Name = 'pCR', m2Name = 'npCR')

##### 3.5 Clinical enrichment analysis ####
response.ce = clinicalEnrichment(maf = changkang_all, clinicalFeature = 'pCR', minMut = 3, useCNV = FALSE, pathways = FALSE)
response.ce$groupwise_comparision$p_value %>% hist()
response.ce$groupwise_comparision[p_value < 0.1]

##### 3.6 Oncogenic Signaling Pathways ####
pws = pathways(maf = changkang, plotType = 'treemap')

# This will produce the figure with correlation of top 25 mutated genes
somaticInteractions(
  maf = laml,
  top = 25,  # Number of top mutated genes to include
  pvalue = c(0.05, 0.1) # Significance levels to display
)



