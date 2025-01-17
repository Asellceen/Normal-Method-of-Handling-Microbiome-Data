##############################
## Working Directory Set Up ##
current_wd<- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_wd)

###############
## Libraries ##
library(ggplot2)
library(dplyr)

###################
## Load the data ##

## original otu
load("otu_zero_prop_tree1.RData")
load("otu_zero_prop_tree2.RData")

## zComposition package
load("results_zComp_tree1.RData")
load("results_zComp_tree2.RData")

## mbImpute package
load("results_mbImpute_tree1.RData")
load("results_mbImpute_tree2.RData")

## phyloMDA package
### 0.5 smoothing method
load("results_phylo_tree1_05smth.RData")
load("results_phylo_tree2_05smth.RData")
### MIX smoothing method
load("Frobenius_error_ZIDTM_tree1.RData")
load("Frobenius_error_ZIDTM_tree2.RData")
load("Simpson_index_ZIDTM_tree1.RData")
load("Simpson_index_ZIDTM_tree2.RData")
load("Wasserstein_Distance_ZIDTM_tree1.RData")
load("Wasserstein_Distance_ZIDTM_tree2.RData")
load("zero_proportion_imupted_ZIDTM_tree1.RData")
load("zero_proportion_imupted_ZIDTM_tree2.RData")
load("prop_zeros_noColZeros_ZIDTM_tree1.RData")
load("prop_zeros_noColZeros_ZIDTM_tree2.RData")
load("column_remove_ZIDTM_tree1.RData")
load("column_remove_ZIDTM_tree2.RData")
load("time_ZIDTM_tree1.RData")
load("time_ZIDTM_tree2.RData")

## ZIPPCA package
load("Frobenius_error_ZIPPCA_tree1.RData")
load("Frobenius_error_ZIPPCA_tree2.RData")
load("Simpson_index_ZIPPCA_tree1.RData")
load("Simpson_index_ZIPPCA_tree2.RData")
load("Wasserstein_Distance_ZIPPCA_tree1.RData")
load("Wasserstein_Distance_ZIPPCA_tree2.RData")
load("zero_proportion_imupted_ZIPPCA_tree1.RData")
load("zero_proportion_imupted_ZIPPCA_tree2.RData")
load("prop_zeros_noColZeros_ZIPPCA_tree1.RData")
load("prop_zeros_noColZeros_ZIPPCA_tree2.RData")
load("column_remove_ZIPPCA_tree1.RData")
load("column_remove_ZIPPCA_tree2.RData")
load("time_ZIPPCA_tree1.RData")
load("time_ZIPPCA_tree2.RData")


#################################################################################
## Comparison: 200 replicates
## Since we use different tree distance calculation, we get two distance matrix. 
## Tree1: sum of true edge's length; Tree2: number of edges connect two taxa. 
#################################################################################

## Mean for Frobenius Error 
mean_zComp_Frb_SQ_tr1<- mean(results_zComp_tree1$Frobenius_error_zComp[1:200,1]) ## zComp-SQ 
mean_zComp_Frb_GBM_tr1<- mean(results_zComp_tree1$Frobenius_error_zComp[1:200,2]) ## zComp-GBM
mean_zComp_Frb_SQ_tr2<- mean(results_zComp_tree2$Frobenius_error_zComp[1:200,1]) ## zComp-SQ 
mean_zComp_Frb_GBM_tr2<- mean(results_zComp_tree2$Frobenius_error_zComp[1:200,2]) ## zComp-GBM
mean_mbImpute_Frb_tr1<- mean(results_mbImpute_tree1$Frobenius_error_mbImpute[1:200,1]) ## mbImpute
mean_mbImpute_Frb_tr2<- mean(results_mbImpute_tree2$Frobenius_error_mbImpute[1:200,1]) ## mbImpute
mean_phyloMDA_05_Frb_tr1<- mean(results_phylo_tree1_05smth$Frobenius_error_ZIDTM[1:200,1]) ## phyloMDA with 0.5 smoothing
mean_phyloMDA_05_Frb_tr2<- mean(results_phylo_tree2_05smth$Frobenius_error_ZIDTM[1:200,1]) ## phyloMDA with 0.5 smoothing
mean_phyloMDA_MIX_Frb_tr1<- mean(unlist(Frobenius_error_ZIDTM_tree1)[1:200], na.rm = T) ## phyloMDA with MIX smoothing
mean_phyloMDA_MIX_Frb_tr2<- mean(unlist(Frobenius_error_ZIDTM_tree2)[1:200], na.rm = T) ## phyloMDA with MIX smoothing
mean_ZIPPCA_Frb_tr1<- mean(unlist(Frobenius_error_ZIPPCA_tree1)[1:200], na.rm = T) ## ZIPPCAlnm 
mean_ZIPPCA_Frb_tr2<- mean(unlist(Frobenius_error_ZIPPCA_tree2)[1:200], na.rm = T) ## ZIPPCAlnm 
Frb_mean_tr1<- round(c(mean_zComp_Frb_SQ_tr1, mean_zComp_Frb_GBM_tr1, mean_mbImpute_Frb_tr1, 
                       mean_phyloMDA_05_Frb_tr1, mean_phyloMDA_MIX_Frb_tr1, mean_ZIPPCA_Frb_tr1), 6); Frb_mean_tr1
Frb_mean_tr2<- round(c(mean_zComp_Frb_SQ_tr2, mean_zComp_Frb_GBM_tr2, mean_mbImpute_Frb_tr2, 
                       mean_phyloMDA_05_Frb_tr2, mean_phyloMDA_MIX_Frb_tr2, mean_ZIPPCA_Frb_tr2), 6); Frb_mean_tr2

## Mean for mean squared error of Simpson's index 
mean_zComp_SmpI_SQ_tr1<- mean(results_zComp_tree1$Simpson_index_zComp[1:200,1]) ## zComp-SQ 
mean_zComp_SmpI_GBM_tr1<- mean(results_zComp_tree1$Simpson_index_zComp[1:200,2]) ## zComp-GBM
mean_zComp_SmpI_SQ_tr2<- mean(results_zComp_tree2$Simpson_index_zComp[1:200,1]) ## zComp-SQ 
mean_zComp_SmpI_GBM_tr2<- mean(results_zComp_tree2$Simpson_index_zComp[1:200,2]) ## zComp-GBM
mean_mbImpute_SmpI_tr1<- mean(results_mbImpute_tree1$Simpson_index_mbImpute[1:200,1]) ## mbImpute
mean_mbImpute_SmpI_tr2<- mean(results_mbImpute_tree2$Simpson_index_mbImpute[1:200,1]) ## mbImpute
mean_phyloMDA_05_SmpI_tr1<- mean(results_phylo_tree1_05smth$Simpson_index_ZIDTM[1:200,1]) ## phyloMDA with 0.5 smoothing
mean_phyloMDA_05_SmpI_tr2<- mean(results_phylo_tree2_05smth$Simpson_index_ZIDTM[1:200,1]) ## phyloMDA with 0.5 smoothing
mean_phyloMDA_MIX_SmpI_tr1<- mean(unlist(Simpson_index_ZIDTM_tree1)[1:200], na.rm = T) ## phyloMDA with MIX smoothing
mean_phyloMDA_MIX_SmpI_tr2<- mean(unlist(Simpson_index_ZIDTM_tree2)[1:200], na.rm = T) ## phyloMDA with MIX smoothing
mean_ZIPPCA_SmpI_tr1<- mean(unlist(Simpson_index_ZIPPCA_tree1)[1:200], na.rm = T) ## ZIPPCAlnm 
mean_ZIPPCA_SmpI_tr2<- mean(unlist(Simpson_index_ZIPPCA_tree2)[1:200], na.rm = T) ## ZIPPCAlnm 
SmpI_mean_tr1<- round(c(mean_zComp_SmpI_SQ_tr1, mean_zComp_SmpI_GBM_tr1, mean_mbImpute_SmpI_tr1, 
                        mean_phyloMDA_05_SmpI_tr1, mean_phyloMDA_MIX_SmpI_tr1, mean_ZIPPCA_SmpI_tr1), 6); SmpI_mean_tr1
SmpI_mean_tr2<- round(c(mean_zComp_SmpI_SQ_tr2, mean_zComp_SmpI_GBM_tr2, mean_mbImpute_SmpI_tr2, 
                        mean_phyloMDA_05_SmpI_tr2, mean_phyloMDA_MIX_SmpI_tr2, mean_ZIPPCA_SmpI_tr2), 6); SmpI_mean_tr2

## Mean for Wasserstein Distance Error
mean_zComp_WasD_SQ_tr1<- mean(results_zComp_tree1$Wasserstein_Distance_zComp[1:200,1]) ## zComp-SQ 
mean_zComp_WasD_GBM_tr1<- mean(results_zComp_tree1$Wasserstein_Distance_zComp[1:200,2]) ## zComp-GBM
mean_zComp_WasD_SQ_tr2<- mean(results_zComp_tree2$Wasserstein_Distance_zComp[1:200,1]) ## zComp-SQ 
mean_zComp_WasD_GBM_tr2<- mean(results_zComp_tree2$Wasserstein_Distance_zComp[1:200,2]) ## zComp-GBM
mean_mbImpute_WasD_tr1<- mean(results_mbImpute_tree1$Wasserstein_Distance_mbImpute[1:200,1]) ## mbImpute
mean_mbImpute_WasD_tr2<- mean(results_mbImpute_tree2$Wasserstein_Distance_mbImpute[1:200,1]) ## mbImpute
mean_phyloMDA_05_WasD_tr1<- mean(results_phylo_tree1_05smth$Wasserstein_Distance_ZIDTM[1:200,1]) ## phyloMDA with 0.5 smoothing
mean_phyloMDA_05_WasD_tr2<- mean(results_phylo_tree2_05smth$Wasserstein_Distance_ZIDTM[1:200,1]) ## phyloMDA with 0.5 smoothing
mean_phyloMDA_MIX_WasD_tr1<- mean(unlist(Wasserstein_Distance_ZIDTM_tree1)[1:200], na.rm = T) ## phyloMDA with MIX smoothing
mean_phyloMDA_MIX_WasD_tr2<- mean(unlist(Wasserstein_Distance_ZIDTM_tree2)[1:200], na.rm = T) ## phyloMDA with MIX smoothing
mean_ZIPPCA_WasD_tr1<- mean(unlist(Wasserstein_Distance_ZIPPCA_tree1)[1:200], na.rm = T) ## ZIPPCAlnm 
mean_ZIPPCA_WasD_tr2<- mean(unlist(Wasserstein_Distance_ZIPPCA_tree2)[1:200], na.rm = T) ## ZIPPCAlnm 
WasD_mean_tr1<- round(c(mean_zComp_WasD_SQ_tr1, mean_zComp_WasD_GBM_tr1, mean_mbImpute_WasD_tr1, 
                        mean_phyloMDA_05_WasD_tr1, mean_phyloMDA_MIX_WasD_tr1, mean_ZIPPCA_WasD_tr1), 6); WasD_mean_tr1
WasD_mean_tr2<- round(c(mean_zComp_WasD_SQ_tr2, mean_zComp_WasD_GBM_tr2, mean_mbImpute_WasD_tr2, 
                        mean_phyloMDA_05_WasD_tr2, mean_phyloMDA_MIX_WasD_tr2, mean_ZIPPCA_WasD_tr2), 6); WasD_mean_tr2

## Processing time for 200 replicates
time_zComp_SQ_tr1<- sum(results_zComp_tree1$zComp_time[1:200, 1])/60 ## zComp-SQ in minutes
time_zComp_GBM_tr1<- sum(results_zComp_tree1$zComp_time[1:200, 2])/60 ## zComp-GBM in minutes
time_zComp_SQ_tr2<- sum(results_zComp_tree2$zComp_time[1:200, 1])/60 ## zComp-SQ in minutes
time_zComp_GBM_tr2<- sum(results_zComp_tree2$zComp_time[1:200, 2])/60 ## zComp-GBM in minutes
time_phyloMDA_05_tr1<- sum(results_phylo_tree1_05smth$ZIDTM_time[1:200,1])/60 ## phyloMDA with 0.5 smoothing in minutes
time_phyloMDA_05_tr2<- sum(results_phylo_tree2_05smth$ZIDTM_time[1:200,1])/60 ## phyloMDA with 0.5 smoothing in minutes
time_phyloMDA_MIX_tr1<- sum(unlist(time_ZIDTM_tree1)[1:200], na.rm = T)/60 ## phyloMDA with MIX smoothing in minutes
time_phyloMDA_MIX_tr2<- sum(unlist(time_ZIDTM_tree2)[1:200], na.rm = T)/60 ## phyloMDA with MIX smoothing in minutes
time_mbImpute_tr1<- sum(results_mbImpute_tree1$mbImpute_time[1:200, 1])/60 ## mbImpute in minutes
time_mbImpute_tr2<- sum(results_mbImpute_tree2$mbImpute_time[1:200, 1])/60 ## mbImpute in minutes
time_ZIPPCA_tr11<- sum(unlist(time_ZIPPCA_tree1)[1:200], na.rm = T)/60 ## ZIPPCAlnm in minutes
time_ZIPPCA_tr22<- sum(unlist(time_ZIPPCA_tree2)[1:200], na.rm = T)/60 ## ZIPPCAlnm in minutes
time_tr1<- round(c(time_zComp_SQ_tr1, time_zComp_GBM_tr1, time_mbImpute_tr1, 
                   time_phyloMDA_05_tr1, time_phyloMDA_MIX_tr1, time_ZIPPCA_tr11), 2); time_tr1
time_tr2<- round(c(time_zComp_SQ_tr2, time_zComp_GBM_tr2, time_mbImpute_tr2, 
                   time_phyloMDA_05_tr2, time_phyloMDA_MIX_tr2, time_ZIPPCA_tr22), 2); time_tr2


###############################
## Output for 200 replicates ##
results_200<- list(Frb_mean_tr1, Frb_mean_tr2, SmpI_mean_tr1, SmpI_mean_tr2, 
                   WasD_mean_tr1, WasD_mean_tr2, time_tr1, time_tr2)
save(results_200, file= "results_200.RData")


######################################################
# Plot of the Mean for different evaluation metrics ##
## Plot for Frobenius Norm Error w/o ZIPPCA method
df <- data.frame(
  Method = rep(c('zComp-SQ', 'zComp-GBM', 'mbImpute', 'phyloMDA-0.5', 'phyloMDA-MIX'), each = 2),
  Tree = rep(c('Tree1', 'Tree2'), times = 5),
  Mean_Frobenius_Error = c(mean_zComp_Frb_SQ_tr1, mean_zComp_Frb_SQ_tr2, 
                           mean_zComp_Frb_GBM_tr1, mean_zComp_Frb_GBM_tr2, 
                           mean_mbImpute_Frb_tr1, mean_mbImpute_Frb_tr2, 
                           mean_phyloMDA_05_Frb_tr1, mean_phyloMDA_05_Frb_tr2,
                           mean_phyloMDA_MIX_Frb_tr1, mean_phyloMDA_MIX_Frb_tr2)
)
ggplot(df, aes(x=Method, y=Mean_Frobenius_Error, fill=Tree)) + 
  geom_bar(stat='identity', position=position_dodge()) +
  geom_text(aes(label=round(Mean_Frobenius_Error, 6)), position=position_dodge(width=1), vjust=-0.5) +
  theme_minimal() +
  labs(title="Mean of Frobenius Norm Error (200 replicates)", x="Methods", y="", fill="Tree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(hjust = 0.5), axis.title.y = element_text(hjust = 0.5))


## Plot for MSE of Simpson's Index w/o ZIPPCA method
df <- data.frame(
  Method = rep(c('zComp-SQ', 'zComp-GBM', 'mbImpute', 'phyloMDA-0.5', 'phyloMDA-MIX'), each = 2),
  Tree = rep(c('Tree1', 'Tree2'), times = 5),
  Mean_MSE_SmpI_Error = c(mean_zComp_SmpI_SQ_tr1, mean_zComp_SmpI_SQ_tr2,
                           mean_zComp_SmpI_GBM_tr1, mean_zComp_SmpI_GBM_tr2,
                           mean_mbImpute_SmpI_tr1, mean_mbImpute_SmpI_tr2,
                           mean_phyloMDA_05_SmpI_tr1, mean_phyloMDA_05_SmpI_tr2, 
                           mean_phyloMDA_MIX_SmpI_tr1, mean_phyloMDA_MIX_SmpI_tr2)
)
ggplot(df, aes(x=Method, y= Mean_MSE_SmpI_Error, fill=Tree)) + 
  geom_bar(stat='identity', position=position_dodge()) +
  geom_text(aes(label=round(Mean_MSE_SmpI_Error, 6)), position=position_dodge(width=1), vjust=-0.5) +
  theme_minimal() +
  labs(title="Mean of the Mean Squared Errors for Simpson's Index (200 replicates)", x="Methods", y="", fill="Tree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(hjust = 0.5), axis.title.y = element_text(hjust = 0.5))


## Plot for Mean for Wasserstein Distance Error w/o ZIPPCA method
df <- data.frame(
  Method = rep(c('zComp-SQ', 'zComp-GBM', 'mbImpute', 'phyloMDA-0.5', 'phyloMDA-MIX'), each = 2),
  Tree = rep(c('Tree1', 'Tree2'), times = 5),
  Mean_WasD_Error = c(mean_zComp_WasD_SQ_tr1, mean_zComp_WasD_SQ_tr2,
                          mean_zComp_WasD_GBM_tr1, mean_zComp_WasD_GBM_tr2,
                          mean_mbImpute_WasD_tr1, mean_mbImpute_WasD_tr2,
                          mean_phyloMDA_05_WasD_tr1, mean_phyloMDA_05_WasD_tr2, 
                          mean_phyloMDA_MIX_WasD_tr1, mean_phyloMDA_MIX_WasD_tr2)
)
ggplot(df, aes(x=Method, y= Mean_WasD_Error, fill=Tree)) + 
  geom_bar(stat='identity', position=position_dodge()) +
  geom_text(aes(label=round(Mean_WasD_Error, 6)), position=position_dodge(width=1), vjust=-0.5) +
  theme_minimal() +
  labs(title="Mean for Wasserstein Distance Error (200 replicates)", x="Methods", y="", fill="Tree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(hjust = 0.5), axis.title.y = element_text(hjust = 0.5))

## Plot for Time for the 200 replicates w/o ZIPPCA method
df <- data.frame(
  Method = rep(c('zComp-SQ', 'zComp-GBM', 'mbImpute', 'phyloMDA-0.5', 'phyloMDA-MIX'), each = 2),
  Tree = rep(c('Tree1', 'Tree2'), times = 5),
  Mean_time = c(time_zComp_SQ_tr1, time_zComp_SQ_tr2,
                      time_zComp_GBM_tr1, time_zComp_GBM_tr2,
                      time_mbImpute_tr1, time_mbImpute_tr2,
                      time_phyloMDA_05_tr1, time_phyloMDA_05_tr2,
                      time_phyloMDA_MIX_tr1, time_phyloMDA_MIX_tr2)
)
ggplot(df, aes(x=Method, y=  Mean_time, fill=Tree)) + 
  geom_bar(stat='identity', position=position_dodge()) +
  geom_text(aes(label=round(Mean_time, 2)), position=position_dodge(width=1), vjust=-0.5) +
  theme_minimal() +
  labs(title="Time for 200 Replications in minutes", x="Methods", y="", fill="Tree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(hjust = 0.5), axis.title.y = element_text(hjust = 0.5))

## Plot for ZIPPCA method with 200 replicates
ZIPPCA_tr1<- round(c(mean_ZIPPCA_Frb_tr1, mean_ZIPPCA_SmpI_tr1, mean_ZIPPCA_WasD_tr1, time_ZIPPCA_tr11), 5); ZIPPCA_tr1
ZIPPCA_tr2<- round(c(mean_ZIPPCA_Frb_tr2, mean_ZIPPCA_SmpI_tr2, mean_ZIPPCA_WasD_tr2, time_ZIPPCA_tr22), 5); ZIPPCA_tr2





#################################################################################
sd_tree1_ZIDTM_MIX_Frob<- sd(unlist(Frobenius_error_ZIDTM_tree1)[1:200], na.rm = T)
sd_tree2_ZIDTM_MIX_Frob<- sd(unlist(Frobenius_error_ZIDTM_tree2)[1:200], na.rm = T)
sd_tree1_ZIDTM_MIX_SI<- sd(unlist(Simpson_index_ZIDTM_tree1)[1:200], na.rm = T)
sd_tree2_ZIDTM_MIX_SI<- sd(unlist(Simpson_index_ZIDTM_tree2)[1:200], na.rm = T)
sd_tree1_ZIDTM_MIX_WasD<- sd(unlist(Wasserstein_Distance_ZIDTM_tree1)[1:200], na.rm = T)
sd_tree2_ZIDTM_MIX_WasD<- sd(unlist(Wasserstein_Distance_ZIDTM_tree2)[1:200], na.rm = T)













