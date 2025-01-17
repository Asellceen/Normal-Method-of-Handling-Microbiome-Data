##############################
## Working Directory Set Up ##
current_wd<- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_wd)

###############
## Libraries ##
library(ape)
library(mbImpute)
library(phyloMDA)
library(phyloseq)
library(zCompositions)
library(ZIPPCAlnm)
library(mvtnorm)
library(eBay)

##################
## load the otu ##
load("simulated_true_abundance.RData")
load("otu_counts_tree2.RData")

######################
## Fitting function ##
###################################
## Fitting Values of four methods: mbImpute, phyloMDA, mbImputeosition, ZIPPCAlnm
## otu: OTUs count dataset
## D_phylo: phylogenetic distance matrix
## phylo_tree: phylogenetic tree (binary)
## ni: sum of OTU counts for each sample
fitting_ZIPPCA<- function(otu){
  
  ## ZIPPCAlnm package
  ## All the setting is default, the number of latent environmental factor is two
  time7<- Sys.time()
  fitting_value<- ZIPPCAlnm(otu)
  if (is.atomic(fitting_value)){
    imputed_value_ZIPPCA<- fitting_value
  } else{
    imputed_value_ZIPPCA<- fitting_value$Q
  }
  time8<- Sys.time()
  time_ZIPPCA<- difftime(time8, time7, units= "secs")
  otu_count_noColZeros<- otu[, which(colSums(otu) != 0)]
  column_remove<- ncol(otu)- ncol(otu_count_noColZeros)
  prop_zeros_noColZeros<- mean(otu_count_noColZeros == 0)
  
  ## Return the results 
  fitting_results<- list(imputed_value_ZIPPCA= imputed_value_ZIPPCA, 
                         prop_zeros_noColZeros= prop_zeros_noColZeros, 
                         column_remove= column_remove, time_ZIPPCA= time_ZIPPCA)
  return(fitting_results)
}

imputed_ZIPPCA_tree2<- list()
prop_zeros_noColZeros_ZIPPCA_tree2<- list()
column_remove_ZIPPCA_tree2<- list()
time_ZIPPCA_tree2<- list()
for (i in 1:500){
  temp_value<- fitting_ZIPPCA(otu_list_tree2[[i]])
  imputed_ZIPPCA_tree2[[i]]<- temp_value$imputed_value_ZIPPCA
  prop_zeros_noColZeros_ZIPPCA_tree2[[i]]<- temp_value$prop_zeros_noColZeros
  column_remove_ZIPPCA_tree2[[i]]<- temp_value$column_remove
  time_ZIPPCA_tree2[[i]]<- temp_value$time_ZIPPCA
}
save(imputed_ZIPPCA_tree2, file= "imputed_ZIPPCA_tree2.RData")
save(prop_zeros_noColZeros_ZIPPCA_tree2, file= "prop_zeros_noColZeros_ZIPPCA_tree2.RData")
save(time_ZIPPCA_tree2, file= "time_ZIPPCA_tree2.RData")
save(column_remove_ZIPPCA_tree2, file= "column_remove_ZIPPCA_tree2.RData")
