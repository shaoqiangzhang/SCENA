## Download file from https://github.com/shaoqiangzhang/scRNAseq_Datasets

library(SCENA)

#read file
Express=read.table("./Ting.GSE51372_readCounts.txt", header = T,row.names = 1)


#do clustering using 5 CPUs
cc=scena_cpu(Express,log=T) ## log=T is to do log-tranformation.  


#compute ARI
library(mclust)
presetlabel=substring(colnames(Express),1,3)## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.7788999 

##plot scatter graph with PCA
plotPCA(Express,cc)