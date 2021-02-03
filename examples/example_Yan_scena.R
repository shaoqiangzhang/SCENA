###An example for running SCENA on Yan's dataset 

#install.packages("parallel")
#install.packages("SNFtool")
#install.packages("apcluster")
#install.packages("mclust")
#install.packages("devtools")
#devtools::install_github("shaoqiangzhang/SCENA")


library(SCENA)

##download data and read file
furl<-"https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/yan/nsmb.2660-S2.csv"
download.file(furl,destfile="./nsmb.2660-S2.csv")

Express=read.csv("./nsmb.2660-S2.csv", header = T,row.names = 1)


Express=Express[,2:91] #select 90 cells,you can omit it if you use full data

##do SCENA clustering
cc=scena(Express,log=T,gpu=F)  
#log=T is to do log-transformation, log=F is no log-transformation
#gpu=F means calling CPU only 

##plot scatter graph with PCA
plotPCA(Express,cc) 

##compute ARI as follows:
library(mclust)
presetlabel=substring(colnames(Express),1,4) ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.8880762
