##install required packages
requiredPackages = c('parallel','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")


## Download file from https://github.com/shaoqiangzhang/scRNAseq_Datasets

library(SCENA)

#read file
Express=read.table("./Ting.GSE51372_readCounts.txt", header = T,row.names = 1)

##preprocess input data
Express=datapreprocess(Express,log=T) #"log=T" is to do log-transformation

#do clustering using 5 CPUs
cc=scena_cpu(Express)

#compute ARI
library(mclust)
presetlabel=substring(colnames(Express),1,3)## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.7788999 

##plot scatter graph with PCA
plotPCA(Express,cc)