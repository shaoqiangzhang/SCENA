
library(SCENA)

#read file
Express=read.table("./Ting.GSE51372_readCounts.txt", header = T,row.names = 1)


#do clustering
cc=scena(Express,log=T, gpu=F) ## log=T is to do log-tranformation,"gpu=F" means to only call CPU.  


#compute ARI
library(mclust)
presetlabel=substring(colnames(Express),1,3)## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.7788999 

##plot scatter graph with PCA
plotPCA(Express,cc)