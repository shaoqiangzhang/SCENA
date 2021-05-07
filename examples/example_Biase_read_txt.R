###An example for running SCENA on the Biase's dataset 

##install packages
requiredPackages = c('parallel','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

## you need download file from https://github.com/shaoqiangzhang/scRNAseq_Datasets
##read scRNA-seq expression file 
Express=read.table("./Biase3celltypes.txt",header = T,row.names = 1)

##data preprocessing
Express=datapreprocess(Express,log=T)  #log=T is to do log-transformation, log=F is no log-transformation

## do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, Express=Express,select_features)
##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.
stopCluster(cl)

##do consensus clustering
cc=consClust(Express) #no parameters if using the predicted number of clusters
#cc=consClust(3) #set the number of clusters = 3

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=rep(c(1:3),c(9,20,20)) ## preset 3 cell types each containing 9 cells, 20 cells, 20 cells,respectively.
#presetlabel=substring(colnames(Express),1,4) ## or read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=1

