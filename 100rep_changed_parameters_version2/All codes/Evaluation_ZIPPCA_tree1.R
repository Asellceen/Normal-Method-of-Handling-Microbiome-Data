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
load("imputed_ZIPPCA_tree1.RData")
load("prop_zeros_noColZeros_ZIPPCA_tree1.RData")
load("time_ZIPPCA_tree1.RData")


################
## Evaluation ##
Frobenius_error_ZIPPCA_tree1<- list()
Simpson_index_ZIPPCA_tree1<- list()
Wasserstein_Distance_ZIPPCA_tree1<- list()
zero_proportion_imupted_ZIPPCA_tree1<- list()
for (i in 1:100){
  n<- nrow(imputed_ZIPPCA_tree1[[i]])
  index1<- which(colSums(otu_list_tree1[[i]]) != 0)
  true_abundance_noColZeros<- true_abundance_combo1[, index1]
  if (any(is.nan(imputed_ZIPPCA_tree1[[i]]))) {
    Frobenius_error_ZIPPCA_tree1[[i]]<- NaN
    Simpson_index_ZIPPCA_tree1[[i]]<- NaN
    Wasserstein_Distance_ZIPPCA_tree1[[i]]<- NaN
    zero_proportion_imupted_ZIPPCA_tree1[[i]]<- NaN
  } else {
    imputed_value_ZIPPCA<- imputed_ZIPPCA_tree1[[i]]
    zero_proportion_imupted_ZIPPCA_tree1[[i]]<- mean(imputed_value_ZIPPCA == 0)
    ## Frobenius norm error
    Frobenius_error_ZIPPCA_tree1[[i]]<- sqrt(sum((imputed_value_ZIPPCA[, index1]- true_abundance_noColZeros)^2))
    ## Simpson's index mean squared error
    Simpson_index_ZIPPCA_tree1[[i]]<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums((imputed_value_ZIPPCA[, index1])^2))^2))/ n
    ## Wasserstein distance
    ## mean of each taxon
    true_abund_noColZeros_mean<- (apply(true_abundance_noColZeros, 2, sum))/n
    ZIPPCA_mean<- (apply(imputed_value_ZIPPCA[, index1], 2, sum))/n
    ## standard deviation of each taxon
    true_abund_noColZeros_sd<- sqrt(apply((true_abundance_noColZeros- true_abund_noColZeros_mean)^2, 2, sum)/ (n-1))
    ZIPPCA_sd<- sqrt(apply((imputed_value_ZIPPCA[, index1]- ZIPPCA_mean)^2, 2, sum)/ (n-1))
    ## mean/ sd in order statstics
    true_abund_noColZeros_mean_sd<- sort(true_abund_noColZeros_mean/ true_abund_noColZeros_sd)
    ZIPPCA_mean_sd<- sort(ZIPPCA_mean/ ZIPPCA_sd)
    ## Mean of the Wasserstein Distance Error between the distribution of mean/sd
    Wasserstein_Distance_ZIPPCA_tree1[[i]]<- sum(abs(true_abund_noColZeros_mean_sd- ZIPPCA_mean_sd)) / ncol(true_abundance_noColZeros)
  }
}

save(Frobenius_error_ZIPPCA_tree1, file= "Frobenius_error_ZIPPCA_tree1.RData")
save(Simpson_index_ZIPPCA_tree1, file= "Simpson_index_ZIPPCA_tree1.RData")
save(Wasserstein_Distance_ZIPPCA_tree1, file= "Wasserstein_Distance_ZIPPCA_tree1.RData")
save(zero_proportion_imupted_ZIPPCA_tree1, file= "zero_proportion_imupted_ZIPPCA_tree1.RData")



