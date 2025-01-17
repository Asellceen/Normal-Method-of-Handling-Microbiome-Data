## mbImpute 
install.packages("devtools", dependencies = T)
install.packages("BiocManager", dependencies = T)
install.packages("glmnet", dependencies = T)
install.packages("Matrix", dependencies = T)
devtools::install_github("ruochenj/mbImpute/mbImpute R package")
## phyloMDA, adaANCOM, eBay
install.packages("MGLM", dependencies = T)
install.packages("plyr", dependencies = T)
install.packages("caper", dependencies = T)
install.packages("genlasso", dependencies = T)
install.packages("magrittr", dependencies = T)
install.packages("foreach", dependencies = T)
install.packages("miLineage", dependencies = T)
install.packages("ape", dependencies = T)
install.packages("ggplot2", dependencies = T)
install.packages("dplyr", dependencies = T)
install.packages("readxl", dependencies = T)
BiocManager::install("phyloseq")
BiocManager::install("ggtree")
devtools::install_github("ZRChao/adaANCOM")
devtools::install_github("liudoubletian/phyloMDA")
devtools::install_github("liudoubletian/eBay")
## zCompositions
install.packages("zCompositions", dependencies = T)
## parallel computing for mbImpute package
install.packages("parallel", dependencies = T)


devtools::install_github("YanyZeng/ZIPPCAlnm") 
