###An example for running SCENA on Zeisel's dataset 
### This example should be run in a Linux system

#install.packages("parallel")
#install.packages("SNFtool")
#install.packages("apcluster")
#install.packages("mclust")
#install.packages("devtools")
#devtools::install_github("shaoqiangzhang/SCENA")

install.packages("gpuR")

library(SCENA)

## download data and read file
furl<-"https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/zeisel/expression_mRNA_17-Aug-2014.txt"
download.file(furl,destfile="./expression_mRNA_17-Aug-2014.txt")
Expr=read.table("./expression_mRNA_17-Aug-2014.txt",skip = 7, header = T, sep="\t")

## preprocess data
Express=Expr[4:nrow(Expr),3:ncol(Expr)] #select 3005 cells,delete other annotations
Express=datapreprocess(Express,log=T) 

## do clustering using 5 CPUs
cc=scena_cpu(Express, log=T)

##plot scatter graph with PCA
plotPCA(Express,cc) 

##compute ARI as follows:
library(mclust)
presetlabel=as.vector(t(Expr[1,3:3007])) ## read the line of cell types as a vector
adjustedRandIndex(presetlabel,as.vector(cc)) ## computer ARI
