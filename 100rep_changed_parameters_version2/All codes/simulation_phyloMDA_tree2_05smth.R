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
## Fitting Values of four methods: mbImpute, phyloMDA, zComposition, ZIPPCAlnm
## otu: OTUs count dataset
## D_phylo: phylogenetic distance matrix
## phylo_tree: phylogenetic tree (binary)
## ni: sum of OTU counts for each sample
fitting_phylo<- function(otu, phylo_tree){
  
  ## phyloMDA package: Zero-Inflated Dirichlet Tree Multinomial Model
  time3<- Sys.time()
  imputed_value_ZIDTM<- eBay_comps(otu_tab= otu, prior = "zero-inflated-Dirichlet-tree",
                                   tree= phylo_tree, model= 0.5)
  time4<- Sys.time()
  time_ZIDTM<- difftime(time4, time3, units= "secs")
  otu_count_noColZeros<- otu[, which(colSums(otu) != 0)]
  column_remove<- ncol(otu)- ncol(otu_count_noColZeros)
  prop_zeros_noColZeros<- mean(otu_count_noColZeros == 0)
  
  ## Return the results 
  fitting_results<- list(imputed_value_ZIDTM= imputed_value_ZIDTM, 
                         prop_zeros_noColZeros= prop_zeros_noColZeros, 
                         column_remove= column_remove, time_ZIDTM= time_ZIDTM)
  return(fitting_results)
}

#############
## Evaluation
evaluation_phylo<- function(otu, true_abundance, imputed_value_ZIDTM){
  
  ## dimensions
  n<- nrow(true_abundance)
  ## index
  index1<- which(colSums(otu) != 0)
  true_abundance_noColZeros<- true_abundance[, index1]
  
  ## Frobenius norm error
  Frobenius_ZIDTM<- sqrt(sum((imputed_value_ZIDTM[, index1]- true_abundance_noColZeros)^2))
  
  ## Simpson's index mean squared error
  Simpson_ZIDTM<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums((imputed_value_ZIDTM[, index1])^2))^2))/ n
  
  ## Wasserstein distance
  ## mean of each taxon
  true_abund_noColZeros_mean<- (apply(true_abundance_noColZeros, 2, sum))/n
  ZIDTM_mean<- (apply(imputed_value_ZIDTM[, index1], 2, sum))/n
  ## standard deviation of each taxon
  true_abund_noColZeros_sd<- sqrt(apply((true_abundance_noColZeros- true_abund_noColZeros_mean)^2, 2, sum)/ (n-1))
  ZIDTM_sd<- sqrt(apply((imputed_value_ZIDTM[, index1]- ZIDTM_mean)^2, 2, sum)/ (n-1))
  ## mean/ sd in order statstics
  true_abund_noColZeros_mean_sd<- sort(true_abund_noColZeros_mean/ true_abund_noColZeros_sd)
  ZIDTM_mean_sd<- sort(ZIDTM_mean/ ZIDTM_sd)
  ## Mean of the Wasserstein Distance Error between the distribution of mean/sd
  ZIDTM_WassD<- sum(abs(true_abund_noColZeros_mean_sd- ZIDTM_mean_sd)) / ncol(true_abundance_noColZeros)
  
  ## Results
  evaluation_results<- list(Frobenius_ZIDTM= Frobenius_ZIDTM, 
                            Simpson_ZIDTM= Simpson_ZIDTM, 
                            ZIDTM_WassD= ZIDTM_WassD)
}


###################################################
## Replicate Function: x means number of replicates
replicate_function<- function(x, otu_list, true_abundance, phylo_tree){
  ## Starting Time
  time1<- Sys.time()
  Frobenius_error_ZIDTM<- array(NA, c(x,1))
  colnames(Frobenius_error_ZIDTM)<- c("phyloMDA-ZIDTM")
  Simpson_index_ZIDTM<- array(NA, c(x, 1))
  colnames(Simpson_index_ZIDTM)<- c("phyloMDA-ZIDTM")
  Wasserstein_Distance_ZIDTM<- array(NA, c(x, 1))
  colnames(Wasserstein_Distance_ZIDTM)<- c("phyloMDA-ZIDTM")
  zero_proportion_imupted_ZIDTM<- array(NA, c(x, 1))
  colnames(zero_proportion_imupted_ZIDTM)<- c("phyloMDA-ZIDTM")
  zero_proportion_otu_noColZeros<- array(NA, c(x, 1))
  colnames(zero_proportion_otu_noColZeros)<- c("zero_proportion_noColumnZeros_otu")
  column_remove_ZIDTM<- array(NA, c(x, 1))
  colnames(column_remove_ZIDTM)<- c("column_removes_ZIDTM")
  ZIDTM_time<- array(NA, c(x, 1))
  colnames(ZIDTM_time)<- c("ZIDTM_time")
  for (i in 1:x) {
    otu<- otu_list[[i]]
    fitting_ZIDTM<- fitting_phylo(otu, phylo_tree)
    imputed_value_ZIDTM<- fitting_ZIDTM$imputed_value_ZIDTM
    results<- evaluation_phylo(otu, true_abundance, imputed_value_ZIDTM)
    Frobenius_error_ZIDTM[i, 1]<- results$Frobenius_ZIDTM
    Simpson_index_ZIDTM[i, 1]<- results$Simpson_ZIDTM
    Wasserstein_Distance_ZIDTM[i, 1]<- results$ZIDTM_WassD
    zero_proportion_otu_noColZeros[i, 1]<- fitting_ZIDTM$prop_zeros_noColZeros
    zero_proportion_imupted_ZIDTM[i, 1]<- mean(imputed_value_ZIDTM == 0)
    column_remove_ZIDTM[i, 1]<- fitting_ZIDTM$column_remove
    ZIDTM_time[i, 1]<- fitting_ZIDTM$time_ZIDTM
  }
  ## End Time
  time2<- Sys.time()
  ## Processing Time
  time<- time2-time1
  output<- list(Frobenius_error_ZIDTM= Frobenius_error_ZIDTM, 
                Simpson_index_ZIDTM= Simpson_index_ZIDTM, 
                Wasserstein_Distance_ZIDTM= Wasserstein_Distance_ZIDTM,
                zero_proportion_otu_noColZeros= zero_proportion_otu_noColZeros,
                zero_proportion_imupted_ZIDTM= zero_proportion_imupted_ZIDTM,
                column_remove_ZIDTM= column_remove_ZIDTM, 
                ZIDTM_time= ZIDTM_time, time= time)
  return(output)
}


results_phylo_tree2_05smth<- replicate_function(100, otu_list= otu_list_tree2,
                                         true_abundance= true_abundance_combo2, 
                                         phylo_tree = COMBO_tree2)
save(results_phylo_tree2_05smth, file= "results_phylo_tree2_05smth.RData")
