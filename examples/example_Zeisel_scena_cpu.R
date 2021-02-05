###An example for running SCENA on Zeisel's dataset 

##install required packages
requiredPackages = c('parallel','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

## download data and read file
furl<-"https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/zeisel/expression_mRNA_17-Aug-2014.txt"
download.file(furl,destfile="./expression_mRNA_17-Aug-2014.txt")
Expr=read.table("./expression_mRNA_17-Aug-2014.txt",skip = 7, header = T, sep="\t")

## preprocess data
Express=Expr[4:nrow(Expr),3:ncol(Expr)] #select cells,delete other annotations
Express=datapreprocess(Express,log=T) 

## do clustering using 5 CPUs
cc=scena_cpu(Express, T=20) "T" is number of iterations   

##plot scatter graph with PCA
plotPCA(Express,cc) 

##compute ARI as follows:
library(mclust)
presetlabel=as.vector(t(Expr[1,3:3007])) ## read the line of cell types as a vector
adjustedRandIndex(presetlabel,as.vector(cc)) ##  ARI=0.7290842 
