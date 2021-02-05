###An example for running SCENA on Yan's dataset 

##install required packages
requiredPackages = c('parallel','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

##download data and read file
furl<-"https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/yan/nsmb.2660-S2.csv"
download.file(furl,destfile="./nsmb.2660-S2.csv")

Express=read.csv("./nsmb.2660-S2.csv", header = T,row.names = 1)

Express=Express[,2:91] #select 90 cells,you can omit it if you use full data

##preprocess input data
Express=datapreprocess(Express,log=T) #"log=T" is to do log-transformation, "log=F" is no log-transformation

##do SCENA clustering uing 5 CPUs
cc=scena_cpu(Express) 
#cc=scena_cpu(Express,  T=20) # or set the number of matrix iterations "T" = 20
#cc=scena_cpu(Express,  T=20, num=6) # or set the number of clusters "num=6"

##plot scatter graph with PCA
plotPCA(Express,cc) 

##compute ARI as follows:
library(mclust)
presetlabel=substring(colnames(Express),1,4) ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.8880762
