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
## Fitting Values of four methods: mbImpute, phyloMDA, zComposition, ZIPPCAlnm
## otu: OTUs count dataset
## D_phylo: phylogenetic distance matrix
## phylo_tree: phylogenetic tree (binary)
## ni: sum of OTU counts for each sample
fitting_zComp<- function(otu){
  
  ## zComposition package
  ## Geometric Bayesian multiplicative (SQ)
  ## Remove the columns with all zeros
  otu_count_noColZeros<- otu[, which(colSums(otu) != 0)]
  time5<- Sys.time()
  ## Updated Fix: default z.warning is 0.8, columns and rows with more than 80% zeros will 
  ## automatically remove. 
  imputed_value_SQ<- as.matrix(cmultRepl(otu_count_noColZeros, method= "SQ", z.warning = 1))
  time6<- Sys.time()
  ## Updated Fix: default z.warning is 0.8, columns and rows with more than 80% zeros will 
  ## automatically remove. 
  imputed_value_GBM<- as.matrix(cmultRepl(otu_count_noColZeros, z.warning = 1))
  time7<- Sys.time()
  time_SQ<- difftime(time6, time5, units= "secs")
  time_GBM<- difftime(time7, time6, units= "secs")
  column_remove<- ncol(otu)- ncol(otu_count_noColZeros)
  prop_zeros_noColZeros<- mean(otu_count_noColZeros == 0)
  
  ## Return the results 
  fitting_results<- list(imputed_value_SQ= imputed_value_SQ, 
                         imputed_value_GBM= imputed_value_GBM,
                         prop_zeros_noColZeros= prop_zeros_noColZeros, 
                         column_remove= column_remove, time_SQ= time_SQ, time_GBM= time_GBM)
  return(fitting_results)
}

## Testing case ##
# fit_test<- fitting_zComp(otu_list_tree1[[1]])

#############
## Evaluation
evaluation_zComp<- function(otu, true_abundance, imputed_value_SQ, imputed_value_GBM){
  ## dimensions
  n<- nrow(true_abundance)
  ## index
  index1<- which(colSums(otu) != 0)
  true_abundance_noColZeros<- true_abundance[, index1]
  
  ## Frobenius norm error
  Frobenius_SQ<- sqrt(sum((imputed_value_SQ- true_abundance_noColZeros)^2))
  Frobenius_GBM<- sqrt(sum((imputed_value_GBM- true_abundance_noColZeros)^2))
  
  ## Simpson's index mean squared error
  Simpson_SQ<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums(imputed_value_SQ^2))^2))/ n
  Simpson_GBM<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums(imputed_value_GBM^2))^2))/ n
  
  ## Wasserstein distance
  ## mean of each taxon
  true_abund_noColZeros_mean<- (apply(true_abundance_noColZeros, 2, sum))/n
  SQ_mean<- (apply(imputed_value_SQ, 2, sum))/n
  GBM_mean<- (apply(imputed_value_GBM, 2, sum))/n
  ## standard deviation of each taxon
  true_abund_noColZeros_sd<- sqrt(apply((true_abundance_noColZeros- true_abund_noColZeros_mean)^2, 2, sum)/ (n-1))
  SQ_sd<- sqrt(apply((imputed_value_SQ- SQ_mean)^2, 2, sum)/ (n-1))
  GBM_sd<- sqrt(apply((imputed_value_GBM- GBM_mean)^2, 2, sum)/ (n-1))
  ## mean/ sd in order statstics
  true_abund_noColZeros_mean_sd<- sort(true_abund_noColZeros_mean/ true_abund_noColZeros_sd)
  SQ_mean_sd<- sort(SQ_mean/ SQ_sd)
  GBM_mean_sd<- sort(GBM_mean/ GBM_sd)
  ## Mean of the Wasserstein Distance Error between the distribution of mean/sd
  SQ_WassD<- sum(abs(true_abund_noColZeros_mean_sd- SQ_mean_sd)) / ncol(true_abundance_noColZeros)
  GBM_WassD<- sum(abs(true_abund_noColZeros_mean_sd- GBM_mean_sd)) / ncol(true_abundance_noColZeros)
  
  ## Results
  evaluation_results<- list(Frobenius_SQ= Frobenius_SQ, Frobenius_GBM= Frobenius_GBM, 
                            Simpson_SQ= Simpson_SQ, Simpson_GBM= Simpson_GBM, 
                            SQ_WassD= SQ_WassD, GBM_WassD= GBM_WassD)
}

## Testing case ##
# evaluation_test<- evaluation_zComp(otu= otu_list_tree1[[1]], 
#                                    true_abundance= true_abundance_combo1, 
#                                    imputed_value_SQ= fit_test$imputed_value_SQ, 
#                                    imputed_value_GBM= fit_test$imputed_value_GBM)
# imputed_value_SQ= fit_test$imputed_value_SQ
# imputed_value_GBM= fit_test$imputed_value_GBM
# index1<- which(colSums(otu_list_tree1[[1]]) != 0)
# true_abundance_noColZeros<- true_abundance_combo1[, index1]
# d_test= imputed_value_SQ-true_abundance_noColZeros

###################################################
## Replicate Function: x means number of replicates
replicate_function<- function(x, otu_list, true_abundance){
  ## Starting Time
  time1<- Sys.time()
  Frobenius_error_zComp<- array(NA, c(x,2))
  colnames(Frobenius_error_zComp)<- c("zComposition-SQ", "zComposition-GBM")
  Simpson_index_zComp<- array(NA, c(x, 2))
  colnames(Simpson_index_zComp)<- c("zComposition-SQ", "zComposition-GBM")
  Wasserstein_Distance_zComp<- array(NA, c(x, 2))
  colnames(Wasserstein_Distance_zComp)<- c("zComposition-SQ", "zComposition-GBM")
  zero_proportion_imupted_zComp<- array(NA, c(x, 2))
  colnames(zero_proportion_imupted_zComp)<- c("zComposition-SQ", "zComposition-GBM")
  zero_proportion_otu_noColZeros<- array(NA, c(x, 1))
  colnames(zero_proportion_otu_noColZeros)<- c("zero_proportion_noColumnZeros_otu")
  column_remove_zComp<- array(NA, c(x, 1))
  colnames(column_remove_zComp)<- c("column_removes")
  zComp_time<- array(NA, c(x, 2))
  colnames(zComp_time)<- c("zComposition-SQ", "zComposition-GBM")
  for (i in 1:x) {
    otu<- otu_list[[i]]
    fitting_zCompstin<- fitting_zComp(otu)
    imputed_value_SQ<- fitting_zCompstin$imputed_value_SQ
    imputed_value_GBM<- fitting_zCompstin$imputed_value_GBM
    results<- evaluation_zComp(otu, true_abundance, imputed_value_SQ, imputed_value_GBM)
    Frobenius_error_zComp[i, 1]<- results$Frobenius_SQ
    Frobenius_error_zComp[i, 2]<- results$Frobenius_GBM
    Simpson_index_zComp[i, 1]<- results$Simpson_SQ
    Simpson_index_zComp[i, 2]<- results$Simpson_GBM
    Wasserstein_Distance_zComp[i, 1]<- results$SQ_WassD
    Wasserstein_Distance_zComp[i, 2]<- results$GBM_WassD
    zero_proportion_otu_noColZeros[i, 1]<- fitting_zCompstin$prop_zeros_noColZeros
    zero_proportion_imupted_zComp[i, 1]<- mean(imputed_value_SQ == 0)
    zero_proportion_imupted_zComp[i, 2]<- mean(imputed_value_GBM == 0)
    column_remove_zComp[i, 1]<- fitting_zCompstin$column_remove
    zComp_time[i, 1]<- fitting_zCompstin$time_SQ
    zComp_time[i, 2]<- fitting_zCompstin$time_GBM
  }
  ## End Time
  time2<- Sys.time()
  ## Processing Time
  time<- time2-time1
  output<- list(Frobenius_error_zComp= Frobenius_error_zComp, 
                Simpson_index_zComp= Simpson_index_zComp, 
                Wasserstein_Distance_zComp= Wasserstein_Distance_zComp,
                zero_proportion_otu_noColZeros= zero_proportion_otu_noColZeros,
                zero_proportion_imupted_zComp= zero_proportion_imupted_zComp,
                column_remove_zComp= column_remove_zComp, 
                zComp_time= zComp_time, time= time)
  return(output)
}


results_zComp_tree1<- replicate_function(100, otu_list= otu_list_tree1,
                                         true_abundance= true_abundance_combo1)
save(results_zComp_tree1, file= "results_zComp_tree1.RData")
