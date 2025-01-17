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

######################
## load the dataset ##
load("simulated_true_abundance.RData")
load("otu_counts_tree1.RData")
load("imputed_ZIDTM_tree1.RData")
load("prop_zeros_noColZeros_ZIDTM_tree1.RData")
load("time_ZIDTM_tree1.RData")


################
## Evaluation ##
Frobenius_error_ZIDTM_tree1<- list()
Simpson_index_ZIDTM_tree1<- list()
Wasserstein_Distance_ZIDTM_tree1<- list()
zero_proportion_imupted_ZIDTM_tree1<- list()
for (i in 1:100){
  n<- nrow(imputed_ZIDTM_tree1[[i]])
  index1<- which(colSums(otu_list_tree1[[i]]) != 0)
  true_abundance_noColZeros<- true_abundance_combo1[, index1]
  if (any(is.nan(imputed_ZIDTM_tree1[[i]]))) {
    Frobenius_error_ZIDTM_tree1[[i]]<- NaN
    Simpson_index_ZIDTM_tree1[[i]]<- NaN
    Wasserstein_Distance_ZIDTM_tree1[[i]]<- NaN
    zero_proportion_imupted_ZIDTM_tree1[[i]]<- NaN
  } else {
    imputed_value_ZIDTM<- imputed_ZIDTM_tree1[[i]]
    zero_proportion_imupted_ZIDTM_tree1[[i]]<- mean(imputed_value_ZIDTM == 0)
    ## Frobenius norm error
    Frobenius_error_ZIDTM_tree1[[i]]<- sqrt(sum((imputed_value_ZIDTM[, index1]- true_abundance_noColZeros)^2))
    ## Simpson's index mean squared error
    Simpson_index_ZIDTM_tree1[[i]]<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums((imputed_value_ZIDTM[, index1])^2))^2))/ n
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
    Wasserstein_Distance_ZIDTM_tree1[[i]]<- sum(abs(true_abund_noColZeros_mean_sd- ZIDTM_mean_sd)) / ncol(true_abundance_noColZeros)
  }
}

save(Frobenius_error_ZIDTM_tree1, file= "Frobenius_error_ZIDTM_tree1.RData")
save(Simpson_index_ZIDTM_tree1, file= "Simpson_index_ZIDTM_tree1.RData")
save(Wasserstein_Distance_ZIDTM_tree1, file= "Wasserstein_Distance_ZIDTM_tree1.RData")
save(zero_proportion_imupted_ZIDTM_tree1, file= "zero_proportion_imupted_ZIDTM_tree1.RData")



