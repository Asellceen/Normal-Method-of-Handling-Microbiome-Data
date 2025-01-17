load(file = "results_100.RData")

###########################################################
## Mean of Frobenius Norm Error ##
Frb_mean_tr1 <- results_100[[1]]
Frb_mean_tr2 <- results_100[[2]]
cFrb <- c(Frb_mean_tr1,Frb_mean_tr2)
Frb_Matrix <- matrix(cFrb, nrow = 6)

###########################################################
## Mean of the Mean Squared Errors for Simpson's Index ##
SmpI_mean_tr1 <- results_100[[3]]
SmpI_mean_tr2 <- results_100[[4]]
cSmpI <- c(SmpI_mean_tr1,SmpI_mean_tr2)
SmpI_Matrix <- matrix(cSmpI, nrow = 6)

###########################################################
## Mean of Wasserstein Distance Error ##
WasD_mean_tr1 <- results_100[[5]]
WasD_mean_tr2 <- results_100[[6]]
cWasD <- c(WasD_mean_tr1,WasD_mean_tr2)
WasD_Matrix <- matrix(cWasD, nrow = 6)

###########################################################
## Time for 1 Replications in minutes##
time_tr1 <- results_100[[7]]
time_tr2 <- results_100[[8]]
ctime <- c(time_tr1,time_tr2)
time_Matrix <- matrix(ctime, nrow = 6)

###########################################################
## Combination of all the matrix ##

Matrix_100_rep = cbind(Frb_Matrix, SmpI_Matrix, WasD_Matrix, time_Matrix)

rownames(Matrix_100_rep) = c("zComp_SQ", "zComp_GBM", "mbImpute", "phyloMDA_05", "phyloMDA_MIX", "ZIPPCA")
colnames(Matrix_100_rep) = c("Frb_mean_tr1", "Frb_mean_tr2", "SmpI_mean_tr1", "SmpI_mean_tr2", 
                         "WasD_mean_tr1", "WasD_mean_tr2", "time_tr1", "time_tr2")

rep100 <- Matrix_100_rep; rep100
save(Matrix_100_rep, file= "Matrix_100_rep.RData")
