# Tau-pathology-drives-dementia-risk-associated-gene-networks-towards-chronic-inflammatory-states-and

This code is companion to (Rexach et al. 2020) "Tau pathology drives dementia risk-associated gene networks towards chronic inflammatory states and immunosuppression"

The code file provides code in R for 1. a consensus WGCNA network analysis (Rexach et al, 2020) that combines Tg4510 Tau mouse purified microglia (data source Wang et al. 2018; PMID 30558641) and TPR50 Tau mouse brain dataset (cortex; 6 months old; data source Swarup et al. 2019; PMID 30510257) and 2. a secondary clustering analysis of selected modules to explore for secondary biological associations  (Rexach et al. 2020).

The repository has one code file (rexach_github_cWGCNA_verified.R) and all necesseary formatted input data files required to run the code in the form of .RData and .csv files.  This includes the following:

1. input sample metadata = input_targets.rda (this code-ready formatted. The published source data is from Swarup et al. 2019; PMID 30510257, Wang et al. 2018; PMID 30558641)

2. input gene expression data = input_datExpr.rda (this code-ready formatted. The published source data is from Swarup et al 2019 PMID 30510257, Wang et al. 2018; PMID 30558641)

3. PPI data downloaded from Biogrid (Stark et al. 2014; PMID 16381927) and Inweb (Rossin et al. 2011; PMID 21249183) and combined into one matrix =  PPI_matrix_BioGrid_Inweb_combined_2014.rda (as decribed in Methods of Swarup et al. 2019; PMID 30510257)

4. consensus kME table = consensusKMEtable.csv (similar to file published in Rexach et al 2020 as Supplementary Table 1; but formatted ready for this code)
