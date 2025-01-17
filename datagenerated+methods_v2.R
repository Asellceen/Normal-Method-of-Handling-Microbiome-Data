## Libraries ##
library(ape)
library(mbImpute)
library(phyloMDA)
library(phyloseq)
library(zCompositions)
library(ZIPPCAlnm)
library(eBay)
library(mvtnorm)

########################
## Parameters Setting ##
# Real dataset: COMBO from phyloMDA package
COMBO_OTU<- as.matrix(t(combo.phyloseq.obj@otu_table@.Data))
COMBO_tree1<- phy_tree(combo.phyloseq.obj) ## Original edge length
COMBO_tree2<- COMBO_tree1
COMBO_tree2$edge.length<- rep(1, 122) ## Edge length equal to 1
# D: phylogenetic distance matrix
D_phylo_combo2<- cophenetic.phylo(COMBO_tree2) ## sum of the number of edges
# alpha should be the scale number
set.seed(1)
alpha_combo<- 5; alpha = alpha_combo
# evolution rate and sigma value
rho<- 0.8
sigma_phy<- 1; sigma = sigma_phy
# mean vector for the multivariate normal distribution
# affect the proportion of sampling zero
K1_combo <- 31
set.seed(2)
theta_combo <- rnorm(62, mean = c(rep(1, K1_combo), rep(5, 62- K1_combo)), sd=1); theta = theta_combo
# total count for each sample
ni_combo<- array(apply(COMBO_OTU, 1, sum), c(98, 1)); Ni = ni_combo; n = 98; k = 62
# number of generated otu table
rep_otu<- 10

##############################
## True Abundance Generator ##
true_abundance<- function(alpha, theta, sigma, rho, D_phylo, n, k, Ni, rep_otu){
  
  ###################
  ## Starting time ##
  ###################
  time1<- Sys.time()
  
  ##############################################################
  ## 2.1: Generate the probability (p_{ij}) that delta_{ij}=1 ##
  ##############################################################
  p_delta1<- array(NA, c(1, k))
  for (j in 1:k){
    p_delta1[1, j]<- 1-exp((-D_phylo[j,k])/alpha)
  }
  
  ###########################################################
  ## 2.2: Generate the delta_param set                     ##
  ## delta_{ij}= 1: Biological Zero (taxa does not exist); ##
  ## delta_{ij}= 0: Otherwise                              ##
  ###########################################################
  delta_param<- array(NA, c(1, k)); set.seed(123)
  for (j in 1:k){
    delta_param[1,j]<- rbinom(1, 1, p_delta1[1,j])
  }
  # proportion of biological zero
  prop.bio<- mean(delta_param == 1)
  
  ##############################################################################
  ## 2.3: Generate the var-cov matrix with the phylogenentic tree information ##
  ## rho: evolution rate; D_phy: phylogenentic tree distance matrix           ##
  ##############################################################################
  var_cov_matrix<- sigma^2* exp(-2* rho* (D_phylo^2))
  
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
  ## testing purpose: check whether the true abundance remains the same
  ## true_abundance_list<- list() 
  ind<- which(delta_param == 0)
  set.seed(123)
  for (i in 1:n){
    ## y_ind<- log(\pi_{ij} / \pi_{ip}); where taxa p is the reference taxa
    ## True abundance with logit scale
    y_ind<- rmvnorm(1, mean= theta[ind], sigma= var_cov_matrix[ind, ind])
    true_abundance[i, ind]<- mu_prob(y_ind)
    colnames(true_abundance)<- colnames(D_phylo)
  }

  otu<- array(0, c(n, k))
  otu_list<- list()
  prop.zero<- c()
  pro.smplng<- c()
  for (l in 1:rep_otu){
    set.seed(l)
    for (i in 1:n){
      ## log(True abundance ratio) ~ multivariate normal
      otu[i, ind]<- rmultinom(1, Ni[i], true_abundance[i, ind])
      colnames(otu)<- colnames(D_phylo)
    }
    otu_list[[l]]<- otu
    
    ## testing purpose: check whether the true abundance remains the same
    ## true_abundance_list[[l]]<- true_abundance
    
    # proportion of zero
    prop.zero[l]<- mean(otu == 0)
    # proportion of sampling zero
    pro.smplng[l]<- prop.zero[l]- prop.bio
  }

  ##############
  ## End time ##
  ##############
  time<- Sys.time()-time1
  
  #############################
  ## 2.6: Return the results ##
  #############################
  results<- list(p_delta1= p_delta1, delta_param= delta_param, time= time, otu_list= otu_list,
                 var_cov_matrix= var_cov_matrix, true_abundance= true_abundance, 
                 prop.zero= prop.zero, prop.bio= prop.bio, pro.smplng= pro.smplng)
  return(results)
}

####################
## True Abundance ##
####################
simulated_combo_data2<- true_abundance(alpha= alpha_combo, theta= theta_combo, 
                                       sigma= sigma_phy, rho= rho, D_phylo= D_phylo_combo2, 
                                       n= 98, k=62, Ni=ni_combo, rep_otu = rep_otu)
true_abundance_combo2<- simulated_combo_data2$true_abundance
otu_list_tree<- simulated_combo_data2$otu_list

###########################################################################
## Error Debug for phyloMDA: zidm methods 
## change the while loop condition from && to |
## remove the the if-else condition and remove case for number of taxa=2
## Error message will comes up when the number of taxa is large
## Modified function: 
est.zidm.EM_modified <- function(X, init.a=NULL, init.pi=NULL, iter=100, conv=1e-6) {
  print("modified function used")
  # for only binary
  Y <- X
  Y <- Y[rowSums(Y)>0, ]
  p <- ncol(Y)
  if(p > 2) {
    est.zidm <- function(Y, convergence=1e-4, iter = 1000) {
      f_t2 <- function(x,child,delta_e,A_e,B_e) {
        J=length(child)-1
        b<-c()
        for(j in 1:(J+1)){
          b<-c(b,x[j])
        }
        fx=0
        I=nrow(delta_e)
        for(i in 1:I){
          for(j in 1:J){
            fx=fx+(1-delta_e[i,child[j]])*(-log(beta(b[j],sum(b[(j+1):(J+1)])))+(b[j]-1)*A_e[i,child[j]]+(sum(b[(j+1):(J+1)])-1)*B_e[i,child[j]])
            #print(fx)
          }
        }
        -fx
      }
      
      g_t2 <- function(x,child,delta_e,A_e,B_e) {
        J=length(child)-1
        b<-c()
        for(j in 1:(J+1)){
          b<-c(b,x[j])
        }
        gx.all<-c()
        
        ### 1
        j=1
        gx=0
        I=nrow(delta_e)
        for(i in 1:I){
          gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(b[j])+A_e[i,child[j]])
        }
        gx.all<-c(gx.all,gx)
        
        ### 2-J
        if(J>1){
          for(k in 2:J){
            gx=0
            for(j in 1:(k-1)){
              for(i in 1:I){
                gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(sum(b[(j+1):(J+1)]))+B_e[i,child[j]])
              }
            }
            j=k
            for(i in 1:I){
              gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(b[j])+A_e[i,child[j]])
            }
            gx.all<-c(gx.all,gx)
          }
        }
        ### J+1
        gx=0
        for(j in 1:J){
          for(i in 1:I){
            gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(sum(b[(j+1):(J+1)]))+B_e[i,child[j]])
          }
        }
        gx.all<-c(gx.all,gx)
        -gx.all
      }
      tree <- ape::rtree(p)
      tree$edge <- cbind(p+1, 1:p)
      tree$Nnode <- 1
      tree$tip.label <- as.character(1:p)
      
      I=nrow(Y)
      M=nrow(tree$edge)
      leaf=length(tree$tip.label)
      
      pai_e=array(rep(0.5,I*M),c(I,M))
      A_e=array(rep(NA,I*M),c(I,M))
      B_e=array(rep(NA,I*M),c(I,M))
      delta_e=array(rep(NA,I*M),c(I,M))
      #delta_e=array(rep(0,I*M),c(I,M))
      a_e=array(rep(5,I*M),c(I,M))
      colnames(pai_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(A_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(B_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(delta_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(a_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      
      Q_a=100
      Q_ab=1000
      num=0
      
      ###########################################################################################
      ## Error found here, change && to or !!!!
      while( (abs(Q_a-Q_ab)/abs(Q_ab))>conv | num<iter){
        
        # print((abs(Q_a-Q_ab)/abs(Q_ab)))
        # num=num+1
        # print(num)
        
        for(i in 1:I){
          for(j in (leaf + 1) : (leaf + tree$Nnode)){
            child=which(tree$edge[, 1] == j)
            for(c in child[1:(length(child)-1)]){
              
              if(Y[i,c]>0){
                delta_e[i,c]=0
              }else{
                mid=beta(a_e[i,c]+Y[i,c],sum(a_e[i,child[(which(child==c)+1):length(child)]]+
                                               Y[i,child[(which(child==c)+1):length(child)]]))/
                  beta(a_e[i,c],sum(a_e[i,child[(which(child==c)+1):length(child)]]))
                delta_e[i,c]=pai_e[i,c]/(pai_e[i,c]+(1-pai_e[i,c])*mid)
              }
            }
          }
        }
        
        for(i in 1:I){
          for(j in (leaf + 1) : (leaf + tree$Nnode)){
            child=which(tree$edge[, 1] == j)
            for(c in child[1:(length(child)-1)]){
              
              A_e[i,c]=digamma(a_e[i,c]+Y[i,c])-
                digamma(sum(a_e[i,child[which(child==c):length(child)]]+Y[i,child[which(child==c):length(child)]]))
              B_e[i,c]=digamma(sum(a_e[i,child[(which(child==c)+1):length(child)]]+Y[i,child[(which(child==c)+1):length(child)]]))-
                digamma(sum(a_e[i,child[which(child==c):length(child)]]+Y[i,child[which(child==c):length(child)]]))
            }
          }
        }
        
        # pai_e[1,]=constrOptim(rep(0.2,M), f_t1, g_t1, ui =rbind(diag(1,M,M),diag(-1,M,M)), ci=c(rep(0.000001,M),rep(-1,M)),
        #                       method="BFGS",tree=tree,delta_e=delta_e,outer.iterations = 100)$par
        
        pai_e[1,]<- c(colMeans(as.matrix(delta_e[,1:(M-1)])),0)
        pai_e[1,]=pai_e[1,]+1e-6
        for(i in 1:I){
          pai_e[i,]=pai_e[1,]
        }
        
        
        for(j in (leaf + 1) : (leaf + tree$Nnode)){
          child=which(tree$edge[, 1] == j)
          H=length(child)
          a_e[1,child]=constrOptim(rep(2,H),f_t2,grad=g_t2,ui = diag(1,H,H),ci=c(rep(0,H)),method="BFGS",
                                   child=child,delta_e=delta_e,A_e=A_e,B_e=B_e)$par
        }
        
        for(i in 2:I){
          a_e[i,]=a_e[1,]
        }
        
        Q_ab=Q_a
        Q_a=0
        
        for(i in 1:I){
          for(j in (leaf + 1) : (leaf + tree$Nnode)){
            child=which(tree$edge[, 1] == j)
            for(c in child[1:(length(child)-1)]){
              Q_a=Q_a+delta_e[i,c]*log(pai_e[i,c])+(1-delta_e[i,c])*log(1-pai_e[i,c])
              +(1-delta_e[i,c])*(-log(beta(a_e[i,c],sum(a_e[i,child[(which(child==c)+1):length(child)]])))+
                                   (a_e[i,c]-1)*A_e[i,c]+
                                   (sum(a_e[i,child[(which(child==c)+1):length(child)]])-1)*B_e[i,c])
            }
          }
        }
      }
      #print("done")
      return(list(alpha=a_e[1,], pi=pai_e[1,], loglik=Q_a))
    }
    res <- est.zidm(Y)
    return(res)
  }else{
    N  <- nrow(Y)
    yl <- Y[, 1]
    yr <- Y[, 2]
    zr <- which(yl==0)
    Ql <- function(a,b,A,B,d) {
      q1 <- sum(d*log(b)+(1-d)*log(1-b))
      q2 <- sum((1-d)*(-lbeta(a[1],a[2]) + (a[1]-1)*A + (a[2]-1)*B))
      q1+q2
    }
    
    # init
    b <- init.pi
    if(is.null(b))  b <-  length(zr)/N
    a <- init.a
    if(is.null(a)) {
      if(sum(yl>0)<2) {
        a <- c(1, 1)
      }else{
        a <-  mom.dm(Y)
      }
    }
    
    if(any(is.infinite(a))) a <- c(1, 1)
    lq <- loglik0 <- loglik_zidm(Y, a, b)
    #lq <- loglik0 <- -Inf
    conv0 <- 1
    delta<- rep(0, N)
    if(b<1e-5) {
      conv0 <- conv/10
      ta <- est.dm.NR(Y)
      a <- ta[[1]]; b <- 0; loglik <- ta[[2]]
    }
    if(any(is.na(a)) | b>0.95) { # can not estimate the value by EM
      conv0 <- conv/10
      ta <- est.dm.NR(Y)
      a <- c(0.5, 0.5); b <- 0; loglik <- -Inf
      if(!any(is.na(ta[[1]]))){
        a<- ta[[1]]; loglik <- ta[[2]]
      }
    }
    Iter <- 0
    ## Error found here, change and to or !!!!
    while(conv0 > conv | Iter < iter) {
      Iter <- Iter+1
      # e-step
      A <- digamma(a[1]+yl) - digamma(sum(a)+rowSums(Y))
      B <- digamma(a[2]+yr) - digamma(sum(a)+rowSums(Y))
      delta[zr] <- 1/(1+(1/b-1)*exp(lbeta(a[1], yr[zr]+a[2])-lbeta(a[1],a[2])))
      
      # m-step
      b <- mean(delta)
      da<- digamma(sum(a))-digamma(a)
      g <- c(N*(1-b)*da[1] + sum((1-delta)*A), N*(1-b)*da[2] + sum((1-delta)*B))
      ta<- trigamma(sum(a))
      H <- N*(1-b)*matrix(c(ta - trigamma(a[1]), ta, ta, ta - trigamma(a[2])),2,2)
      t1 <- try(H1 <- solve(H), silent = T)
      if(class(t1)[[1]]=='try-error') H1 <- matrix(c(1,0,0,1), 2)
      a <- a - as.numeric(H1%*%t(t(g)))
      a[a<0] <- 0.5
      loglik1 <- loglik_zidm(Y, a, b)
      #loglik1 <- Ql(a, b, A, B, delta)
      lq <- c(lq, loglik1)
      conv0 <- abs(loglik1-loglik0)/abs(loglik1)
      loglik0 <-loglik1
    }
    res <- list(alpha=a, pi=b, loglik=loglik0)
    return(res)
  }
}
environment(est.zidm.EM_modified) <- asNamespace("adaANCOM")
assignInNamespace("est.zidm.EM",est.zidm.EM_modified,ns="adaANCOM")


##############
## mbImpute ##
###########################################################################
fitting_mbImpute<- function(otu, D_phylo){
  
  ## mbImpute package: Gamma-Normal Mixed Model
  time1<- Sys.time()
  otu_mbImpute<- mbImpute(otu_tab= otu, D= D_phylo)$imp_count_mat_origlibsize
  time2<- Sys.time()
  time_mbImpute<- difftime(time2, time1, units= "secs")
  ## Change the otu counts into proportion
  imputed_value_mbImpute<- array(NA, c(n, k))
  for (i in 1:n) {
    temp_value<- sum(otu_mbImpute[i,])
    for (j in 1:k) {
      imputed_value_mbImpute[i, j]<- otu_mbImpute[i, j] / temp_value
    }
  }
  ## Return the results 
  mbImpute_results<- list(imputed_value_mbImpute= imputed_value_mbImpute, 
                          time_mbImpute= time_mbImpute)
  return(mbImpute_results)
}

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
  mbImpute_time<- array(NA, c(x, 1))
  colnames(mbImpute_time)<- c("mbImpute_time")
  for (i in 1:x) {
    otu<- simulated_combo_data2$otu_list[[i]]
    fitting_mbImpute<- fitting_mbImpute(otu, D_phylo)
    imputed_value_mbImpute<- fitting_mbImpute$imputed_value_mbImpute
    results<- evaluation_mbImpute(otu, true_abundance, imputed_value_mbImpute)
    Frobenius_error_mbImpute[i, 1]<- results$Frobenius_mbImpute
    Simpson_index_mbImpute[i, 1]<- results$Simpson_mbImpute
    Wasserstein_Distance_mbImpute[i, 1]<- results$mbImpute_WassD
    zero_proportion_imupted_mbImpute[i, 1]<- mean(imputed_value_mbImpute == 0)
    mbImpute_time[i, 1]<- fitting_mbImpute$time_mbImpute
  }
  ## End Time
  time2<- Sys.time()
  ## Processing Time
  time<- time2-time1
  output<- list(Frobenius_error_mbImpute= Frobenius_error_mbImpute, 
                Simpson_index_mbImpute= Simpson_index_mbImpute, 
                Wasserstein_Distance_mbImpute= Wasserstein_Distance_mbImpute,
                mbImpute_time= mbImpute_time, time= time)
  return(output)
}
results_mbImpute<- replicate_function(10, n=98, k= 62, ni=ni_combo, otu_list= otu_list_tree,
                                      true_abundance= true_abundance_combo2, D_phylo= D_phylo_combo2)
save(results_mbImpute, file= "result_mbImpute05.RData")


################
##   ZIPPCA   ##
###########################################################################
fitting_ZIPPCA<- function(otu){
  ## ZIPPCAlnm package
  ## All the setting is default, the number of latent environmental factor is two
  time8<- Sys.time()
  fitting_value<- ZIPPCAlnm(otu)
  if (is.atomic(fitting_value)){
    imputed_value_ZIPPCA<- fitting_value
  } else{
    imputed_value_ZIPPCA<- fitting_value$Q
  }
  time9<- Sys.time()
  time_ZIPPCA<- difftime(time9, time8, units= "secs")
  
  ## Return the results 
  ZIPPCA_results<- list(imputed_value_ZIPPCA= imputed_value_ZIPPCA, 
                        time_ZIPPCA= time_ZIPPCA)
  return(ZIPPCA_results)
}

imputed_ZIPPCA_tree2<- list()
time_ZIPPCA_tree2<- list()
for (i in 1:200){
  temp_value<- fitting_ZIPPCA(simulated_combo_data2$otu_list[[i]])
  imputed_ZIPPCA_tree2[[i]]<- temp_value$imputed_value_ZIPPCA
  time_ZIPPCA_tree2[[i]]<- temp_value$time_ZIPPCA
}

Frobenius_error_ZIPPCA_tree2<- list()
Simpson_index_ZIPPCA_tree2<- list()
Wasserstein_Distance_ZIPPCA_tree2<- list()
zero_proportion_imupted_ZIPPCA_tree2<- list()
for (i in 1:200){
  n<- nrow(imputed_ZIPPCA_tree2[[i]])
  index1<- which(colSums(otu_list_tree[[i]]) != 0)
  true_abundance_noColZeros<- true_abundance_combo2[, index1]
  if (any(is.nan(imputed_ZIPPCA_tree2[[i]]))) {
    Frobenius_error_ZIPPCA_tree2[[i]]<- NaN
    Simpson_index_ZIPPCA_tree2[[i]]<- NaN
    Wasserstein_Distance_ZIPPCA_tree2[[i]]<- NaN
    zero_proportion_imupted_ZIPPCA_tree2[[i]]<- NaN
  } else {
    imputed_value_ZIPPCA<- imputed_ZIPPCA_tree2[[i]]
    zero_proportion_imupted_ZIPPCA_tree2[[i]]<- mean(imputed_value_ZIPPCA == 0)
    ## Frobenius norm error
    Frobenius_error_ZIPPCA_tree2[[i]]<- sqrt(sum((imputed_value_ZIPPCA[, index1]- true_abundance_noColZeros)^2))
    ## Simpson's index mean squared error
    Simpson_index_ZIPPCA_tree2[[i]]<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums((imputed_value_ZIPPCA[, index1])^2))^2))/ n
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
    Wasserstein_Distance_ZIPPCA_tree2[[i]]<- sum(abs(true_abund_noColZeros_mean_sd- ZIPPCA_mean_sd)) / ncol(true_abundance_noColZeros)
  }
}
result_ZIPPCA <- list(imputed_ZIPPCA_tree2= imputed_ZIPPCA_tree2, time_ZIPPCA_tree2= time_ZIPPCA_tree2, 
                      Frobenius_error_ZIPPCA_tree2= Frobenius_error_ZIPPCA_tree2, Simpson_index_ZIPPCA_tree2= Simpson_index_ZIPPCA_tree2,
                      Wasserstein_Distance_ZIPPCA_tree2= Wasserstein_Distance_ZIPPCA_tree2, zero_proportion_imupted_ZIPPCA_tree2= zero_proportion_imupted_ZIPPCA_tree2)
save(result_ZIPPCA, file= "result_ZIPPCA05.RData")



##############
## phyloMDA ##
###########################################################################
fitting_phylo<- function(otu, phylo_tree) {
  
  ## Initialize the result variable as NaN
  imputed_value_ZIDTM <- NaN
  
  time3<- Sys.time()
  ## Use tryCatch() to handle errors
  tryCatch({
    ## phyloMDA package: Zero-Inflated Dirichlet Tree Multinomial Model
    imputed_value_ZIDTM <- eBay_comps(otu_tab= otu, 
                                      prior = "zero-inflated-Dirichlet-tree",
                                      tree= phylo_tree, 
                                      model= "MIX") ## NOTE!!: model means smooth methods, may use MIX or 0.5
  }, error = function(e) {
    ## Print "NaN" and return when an error occurs
    print("NaN")
    return()
  })
  time4<- Sys.time()
  time_ZIDTM <- difftime(time4, time3, units= "secs")
  
  ## Return the results 
  phyloMDA_results <- list(imputed_value_ZIDTM = imputed_value_ZIDTM, 
                           time_ZIDTM = time_ZIDTM)
  
  return(phyloMDA_results)
}

imputed_ZIDTM_tree2<- list()
time_ZIDTM_tree2<- list() 
for (i in 1:10){
  temp_value<- fitting_phylo(simulated_combo_data2$otu_list[[i]], COMBO_tree2)
  imputed_ZIDTM_tree2[[i]]<- temp_value$imputed_value_ZIDTM
  time_ZIDTM_tree2[[i]]<- temp_value$time_ZIDTM
}

Frobenius_error_ZIDTM_tree2<- list()
Simpson_index_ZIDTM_tree2<- list()
Wasserstein_Distance_ZIDTM_tree2<- list()
zero_proportion_imupted_ZIDTM_tree2<- list()
for (i in 1:10){
  n<- nrow(imputed_ZIDTM_tree2[[i]])
  index1<- which(colSums(simulated_combo_data2$otu_list[[i]]) != 0)
  true_abundance_noColZeros<- true_abundance_combo2[, index1]
  if (any(is.nan(imputed_ZIDTM_tree2[[i]]))) {
    Frobenius_error_ZIDTM_tree2[[i]]<- NaN
    Simpson_index_ZIDTM_tree2[[i]]<- NaN
    Wasserstein_Distance_ZIDTM_tree2[[i]]<- NaN
    zero_proportion_imupted_ZIDTM_tree2[[i]]<- NaN
  } else {
    imputed_value_ZIDTM<- imputed_ZIDTM_tree2[[i]]
    zero_proportion_imupted_ZIDTM_tree2[[i]]<- mean(imputed_value_ZIDTM == 0)
    ## Frobenius norm error
    Frobenius_error_ZIDTM_tree2[[i]]<- sqrt(sum((imputed_value_ZIDTM[, index1]- true_abundance_noColZeros)^2))
    ## Simpson's index mean squared error
    Simpson_index_ZIDTM_tree2[[i]]<- (sum((rowSums(true_abundance_noColZeros^2)- rowSums((imputed_value_ZIDTM[, index1])^2))^2))/ n
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
    Wasserstein_Distance_ZIDTM_tree2[[i]]<- sum(abs(true_abund_noColZeros_mean_sd- ZIDTM_mean_sd)) / ncol(true_abundance_noColZeros)
  }
}
result_ZIDTM <- list(imputed_ZIDTM_tree2= imputed_ZIDTM_tree2, time_ZIDTM_tree2= time_ZIDTM_tree2, 
                     Frobenius_error_ZIDTM_tree2= Frobenius_error_ZIDTM_tree2, Simpson_index_ZIDTM_tree2= Simpson_index_ZIDTM_tree2, 
                     Wasserstein_Distance_ZIDTM_tree2= Wasserstein_Distance_ZIDTM_tree2,
                     zero_proportion_imupted_ZIDTM_tree2= zero_proportion_imupted_ZIDTM_tree2)
save(result_ZIDTM, file= "result_ZIDTM05.RData")

###############
##   zComp   ##
###########################################################################
fitting_zComp<- function(otu){
  
  ## zComposition package
  ## Geometric Bayesian multiplicative (SQ)
  ## Remove the columns with all zeros
  otu_count_noColZeros<- otu[, which(colSums(otu) != 0)]
  time5<- Sys.time()
  ## Updated Fix: default z.warning is 0.8, columns and rows with more than 80% zeros will automatically remove. 
  imputed_value_SQ<- as.matrix(cmultRepl(otu_count_noColZeros, method= "SQ", z.warning = 1))
  time6<- Sys.time()
  ## Updated Fix: default z.warning is 0.8, columns and rows with more than 80% zeros will automatically remove. 
  imputed_value_GBM<- as.matrix(cmultRepl(otu_count_noColZeros, z.warning = 1))
  time7<- Sys.time()
  time_SQ<- difftime(time6, time5, units= "secs")
  time_GBM<- difftime(time7, time6, units= "secs")
  
  ## Return the results 
  zComp_results<- list(imputed_value_SQ= imputed_value_SQ, 
                       imputed_value_GBM= imputed_value_GBM,
                       time_SQ= time_SQ, time_GBM= time_GBM)
  return(zComp_results)
}

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
    otu<- simulated_combo_data2$otu_list[[i]]
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
    zero_proportion_imupted_zComp[i, 1]<- mean(imputed_value_SQ == 0)
    zero_proportion_imupted_zComp[i, 2]<- mean(imputed_value_GBM == 0)
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
                zero_proportion_imupted_zComp= zero_proportion_imupted_zComp,
                zComp_time= zComp_time, time= time)
  return(output)
}
results_zComp_tree2<- replicate_function(10, otu_list= otu_list_tree,
                                         true_abundance= true_abundance_combo2)
save(results_zComp_tree2, file= "result_zComp05.RData")







################
## save image ##
################
save.image("method fit1.RData")