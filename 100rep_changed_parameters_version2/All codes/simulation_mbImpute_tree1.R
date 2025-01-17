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
load("otu_counts_tree1.RData")

######################

## Fitting function ##
###################################
## Fitting Values of four methods: mbImpute, phyloMDA, mbImputeosition, ZIPPCAlnm
## otu: OTUs count dataset
## D_phylo: phylogenetic distance matrix
## phylo_tree: phylogenetic tree (binary)
## ni: sum of OTU counts for each sample
fitting_mbImpute<- function(n, k, ni, otu, D_phylo){
  
  ## mbImpute package: Gamma-Normal Mixed Model
  time1<- Sys.time()
  otu_mbImpute<- mbImpute(otu_tab= otu, D= D_phylo)$imp_count_mat_origlibsize
  time2<- Sys.time()
  time_mbImpute<- difftime(time2, time1, units= "secs")
  ## Change the otu counts into proportion
  imputed_value_mbImpute<- array(NA, c(n, k))
  for (i in 1:n) {
    for (j in 1:k) {
      imputed_value_mbImpute[i, j]<- otu_mbImpute[i, j] / ni[i]
    }
  }
  otu_count_noColZeros<- otu[, which(colSums(otu) != 0)]
  column_remove<- ncol(otu)- ncol(otu_count_noColZeros)
  prop_zeros_noColZeros<- mean(otu_count_noColZeros == 0)
  
  ## Return the results 
  fitting_results<- list(imputed_value_mbImpute= imputed_value_mbImpute, 
                         prop_zeros_noColZeros= prop_zeros_noColZeros, 
                         column_remove= column_remove, time_mbImpute= time_mbImpute)
  return(fitting_results)
}

#############
## Evaluation
evaluation_mbImpute<- function(otu, true_abundance, imputed_value_mbImpute){
  
  ## dimensions
  n<- nrow(true_abundance)
  ## index
  index1<- which(colSums(otu) != 0)
  true_abundance_noColZeros<- true_abundance[, index1]
  
  ## Frobenius norm error
  Frobenius_mbImpute<- sqrt(sum((imputed_value_mbImpute[, index1]- true_abundance_noColZeros)^2))
  
  ## Simpson's index mean squared error
  Simpson_mbImpute<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums((imputed_value_mbImpute[, index1])^2))^2))/ n
  
  ## Wasserstein distance
  ## mean of each taxon
  true_abund_noColZeros_mean<- (apply(true_abundance_noColZeros, 2, sum))/n
  mbImpute_mean<- (apply(imputed_value_mbImpute[, index1], 2, sum))/n
  ## standard deviation of each taxon
  true_abund_noColZeros_sd<- sqrt(apply((true_abundance_noColZeros- true_abund_noColZeros_mean)^2, 2, sum)/ (n-1))
  mbImpute_sd<- sqrt(apply((imputed_value_mbImpute[, index1]- mbImpute_mean)^2, 2, sum)/ (n-1))
  ## mean/ sd in order statstics
  true_abund_noColZeros_mean_sd<- sort(true_abund_noColZeros_mean/ true_abund_noColZeros_sd)
  mbImpute_mean_sd<- sort(mbImpute_mean/ mbImpute_sd)
  ## Mean of the Wasserstein Distance Error between the distribution of mean/sd
  mbImpute_WassD<- sum(abs(true_abund_noColZeros_mean_sd- mbImpute_mean_sd)) / ncol(true_abundance_noColZeros)
  
  ## Results
  evaluation_results<- list(Frobenius_mbImpute= Frobenius_mbImpute, 
                            Simpson_mbImpute= Simpson_mbImpute, 
                            mbImpute_WassD= mbImpute_WassD)
}

###################################################
## Replicate Function: x means number of replicates
replicate_function<- function(x, n, k, ni, otu_list, true_abundance, D_phylo){
  ## Starting Time
  time1<- Sys.time()
  Frobenius_error_mbImpute<- array(NA, c(x,1))
  colnames(Frobenius_error_mbImpute)<- c("mbImpute")
  Simpson_index_mbImpute<- array(NA, c(x, 1))
  colnames(Simpson_index_mbImpute)<- c("mbImpute")
  Wasserstein_Distance_mbImpute<- array(NA, c(x, 1))
  colnames(Wasserstein_Distance_mbImpute)<- c("mbImpute")
  zero_proportion_imupted_mbImpute<- array(NA, c(x, 1))
  colnames(zero_proportion_imupted_mbImpute)<- c("mbImpute")
  zero_proportion_otu_noColZeros<- array(NA, c(x, 1))
  colnames(zero_proportion_otu_noColZeros)<- c("zero_proportion_noColumnZeros_otu")
  column_remove_mbImpute<- array(NA, c(x, 1))
  colnames(column_remove_mbImpute)<- c("column_removes_mbImpute")
  mbImpute_time<- array(NA, c(x, 1))
  colnames(mbImpute_time)<- c("mbImpute_time")
  for (i in 1:x) {
    otu<- otu_list[[i]]
    fitting_mbImpute<- fitting_mbImpute(n, k, ni, otu, D_phylo)
    imputed_value_mbImpute<- fitting_mbImpute$imputed_value_mbImpute
    results<- evaluation_mbImpute(otu, true_abundance, imputed_value_mbImpute)
    Frobenius_error_mbImpute[i, 1]<- results$Frobenius_mbImpute
    Simpson_index_mbImpute[i, 1]<- results$Simpson_mbImpute
    Wasserstein_Distance_mbImpute[i, 1]<- results$mbImpute_WassD
    zero_proportion_otu_noColZeros[i, 1]<- fitting_mbImpute$prop_zeros_noColZeros
    zero_proportion_imupted_mbImpute[i, 1]<- mean(imputed_value_mbImpute == 0)
    column_remove_mbImpute[i, 1]<- fitting_mbImpute$column_remove
    mbImpute_time[i, 1]<- fitting_mbImpute$time_mbImpute
  }
  ## End Time
  time2<- Sys.time()
  ## Processing Time
  time<- time2-time1
  output<- list(Frobenius_error_mbImpute= Frobenius_error_mbImpute, 
                Simpson_index_mbImpute= Simpson_index_mbImpute, 
                Wasserstein_Distance_mbImpute= Wasserstein_Distance_mbImpute,
                zero_proportion_otu_noColZeros= zero_proportion_otu_noColZeros,
                zero_proportion_imupted_mbImpute= zero_proportion_imupted_mbImpute,
                column_remove_mbImpute= column_remove_mbImpute, 
                mbImpute_time= mbImpute_time, time= time)
  return(output)
}

results_mbImpute_tree1<- replicate_function(100, n=98, k= 62, ni=ni_combo, otu_list= otu_list_tree1,
                                            true_abundance= true_abundance_combo1, D_phylo= D_phylo_combo1)
save(results_mbImpute_tree1, file= "results_mbImpute_tree1.RData")
