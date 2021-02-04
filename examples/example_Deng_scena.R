
## Download Deng.zip from https://github.com/shaoqiangzhang/scRNAseq_Datasets
## unzip file to your working directory

library(SCENA)

#read file
Express=read.table("./Deng.txt", header = T,row.names = 1)

#do clustering using 5 CPUs
cc=scena_cpu(Express,log=T) ## log=T is to do log-tranformation  


#compute ARI
library(mclust)
presetlabel=c(substring(colnames(Express),1,6)[1:264],"zy","zy","zy","zy")## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.4587225

##plot scatter graph with PCA
plotPCA(Express,cc)