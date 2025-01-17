##############################
## Working Directory Set Up ##
current_wd<- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_wd)

#################################
## Load the relative abundance ##
#################################
load("simulated_true_abundance.RData")

#####################
## Generate OTU Count
## true_abundance: a fixed true_abundance that generate from previous function
## ni: sum of OTU counts for each sample
otu_generator<- function(n, k, true_abundance, ni){
  otu<- array(0, c(n, k))
  colnames(otu)<- colnames(true_abundance)
  for (i in 1:n) {
    ind<- which(true_abundance[i, ] != 0)
    otu[i, ind]<- rmultinom(1, ni[i], true_abundance[i, ind])
  }
  prop_zeros<- mean(otu == 0)
  return(list(otu= otu, prop_zeros= prop_zeros))
}
otu_list_tree1<- list()
otu_zero_prop_tree1<- list()
for (i in 1:100){
  set.seed(i)
  otu_list_tree1[[i]]<- otu_generator(n= 98, k= 62, true_abundance = true_abundance_combo1, ni= ni_combo)$otu
  otu_zero_prop_tree1[[i]]<- otu_generator(n= 98, k= 62, true_abundance = true_abundance_combo1, ni= ni_combo)$prop_zeros
}
save(otu_list_tree1, file= "otu_counts_tree1.RData")
save(otu_zero_prop_tree1, file= "otu_zero_prop_tree1.RData")

