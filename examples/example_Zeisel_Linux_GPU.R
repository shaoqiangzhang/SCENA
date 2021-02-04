###An example for running SCENA on Zeisel's dataset with GPU
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

## do clustering using GPU
library(gpuR)
source('./ApSpe_GPU.R') ## download the "ApSpe_GPU.R" file to the working directory, and source it 
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=20,T=300,X1=50,X2=100,X3=150,X4=200,X5=250,Express=Express,select_features_GPU)##see Note3
stopCluster(cl)

##do consensus clustering
cc=consClust() #no parameters if using the predicted number of clusters

##plot scatter graph with PCA
plotPCA(Express,cc) 

##compute ARI as follows:
library(mclust)
presetlabel=as.vector(t(Expr[1,3:3007])) ## read the line of cell types as a vector
adjustedRandIndex(presetlabel,as.vector(cc)) ## computer ARI
