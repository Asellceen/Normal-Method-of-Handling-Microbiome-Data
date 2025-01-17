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

########################
## Parameters Setting ##
# Real dataset: COMBO from phyloMDA package
COMBO_OTU<- as.matrix(t(combo.phyloseq.obj@otu_table@.Data))
COMBO_tree1<- phy_tree(combo.phyloseq.obj) ## Original edge length
COMBO_tree2<- COMBO_tree1
COMBO_tree2$edge.length<- rep(1, 122) ## Edge length equal to 1
# D: phylogenetic distance matrix
D_phylo_combo1<- cophenetic.phylo(COMBO_tree1) ## sum of all true edge lengths
D_phylo_combo2<- cophenetic.phylo(COMBO_tree2) ## sum of the number of edges
# alpha0_j and alpha1: vary the proportion of structural zeros
set.seed(1)
alpha_combo<- rnorm(ncol(COMBO_OTU), mean = 0.5, sd = 1.5) ##### CHANGE MEAN FOMR 0 TO 0.5, SD FROM 1 TO 1.5
# evolution rate and sigma value
rho<- 0.5 ################# FROM 0.25 TO 0.5, HIGHER EVOLUTION RATE BETWEEN TWO TAXA
sigma_phy<- 1.27  ############ FROM 1 TO 1.27
# mean vector for the multivariate normal distribution
# affect the proportion of sampling zero
K1_combo <- floor(62* 0.85) 
set.seed(2)
theta_combo <- rnorm(62, mean = c(rep(0.3, K1_combo), rep(7, 62- K1_combo)), sd=1) ####### CHANGE FROM 0.1 TO 0.3
# total count for each sample
ni_combo<- array(apply(COMBO_OTU, 1, sum), c(98, 1))

##############################
## True Abundance Generator ##
true_abundance<- function(alpha, theta, sigma, rho, D_phylo, n, k, Ni){
  
  ###################
  ## Starting time ##
  ###################
  time1<- Sys.time()
  
  ##############################################################
  ## 2.1: Generate the probability (p_{ij}) that delta_{ij}=1 ##
  ##############################################################
  p_delta1<- array(NA, c(1, k))
  for (j in 1:k){
    p_delta1[1, j]<- 1/ (1+ exp(-(alpha[j])))
  }
  
  ###########################################################
  ## 2.2: Generate the delta_param set                     ##
  ## delta_{ij}= 1: Structural Zero (taxa does not exist); ##
  ## delta_{ij}= 0: Otherwise                              ##
  ###########################################################
  delta_param<- array(NA, c(n, k))
  for (i in 1:n){
    for (j in 1:k){
      delta_param[i,j]<- rbinom(1, 1, p_delta1[1,j])
    }
  }
  
  ##############################################################################
  ## 2.3: Generate the var-cov matrix with the phylogenentic tree information ##
  ## rho: evolution rate; D_phy: phylogenentic tree distance matrix           ##
  ##############################################################################
  var_cov_matrix<- sigma^2* exp(-2* rho* D_phylo)
  
  ###############################################################################
  ## 2.4: Logit-normal-distribution function based on log-ratio transformation ##
  ###############################################################################
  mu_prob<- function(x){
    pl<- length(x)
    tmp<- sum(exp(x[1: (pl-1)]))
    ans<- c(exp(x[1: (pl-1)])/ (1+ tmp), 1/ (1+ tmp))
  }
  
  ######################################
  ## 2.5: Generate the true abundance ##
  ######################################
  true_abundance<- array(0, c(n, k))
  for (i in 1:n) {
    ind<- which(delta_param[i, ]== 0)
    ## y_ind<- log(\pi_{ij} / \pi_{ip}); where taxa p is the reference taxa
    ## True abundance with logit scale
    y_ind<- rmvnorm(1, mean= theta[ind], sigma= var_cov_matrix[ind, ind])
    true_abundance[i, ind]<- mu_prob(y_ind)
  }
  
  ##############
  ## End time ##
  ##############
  time<- Sys.time()-time1
  
  #############################
  ## 2.6: Return the results ##
  #############################
  results<- list(p_delta1= p_delta1, delta_param= delta_param, time= time,
                 var_cov_matrix= var_cov_matrix, true_abundance= true_abundance)
  return(results)
}

####################
## True Abundance ##
####################
set.seed(3)
simulated_combo_data1<- true_abundance(alpha= alpha_combo, theta= theta_combo, 
                                 sigma= sigma_phy, rho= rho, D_phylo= D_phylo_combo1, 
                                 n= 98, k=62, Ni=ni_combo)
simulated_combo_data2<- true_abundance(alpha= alpha_combo, theta= theta_combo, 
                                 sigma= sigma_phy, rho= rho, D_phylo= D_phylo_combo2, 
                                 n= 98, k=62, Ni=ni_combo)
true_abundance_combo1<- simulated_combo_data1$true_abundance
colnames(true_abundance_combo1)<- colnames(COMBO_OTU)
true_abundance_combo2<- simulated_combo_data2$true_abundance
colnames(true_abundance_combo2)<- colnames(COMBO_OTU)

###############
## Save Data ##
###############
save.image("simulated_true_abundance.RData")

