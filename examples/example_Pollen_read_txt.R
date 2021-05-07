###An example for running SCENA on Pollen's dataset 

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
furl<-"https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/pollen/NBT_hiseq_linear_tpm_values.txt"
download.file(furl,destfile="./NBT_hiseq_linear_tpm_values.txt")
Express=read.table("./NBT_hiseq_linear_tpm_values.txt", header = T,row.names = 1)

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
#cc=consClust(Express,11) #set the number of clusters = 11

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=substring(colnames(Express),1,5)## read the column names as preset cell types
#presetlabel=rep(c(1:11),c(22,17,11,37,31,54,24,40,24,15,26)) ## preset 11 cell types 
adjustedRandIndex(presetlabel,as.vector(cc)) ## 'cc' is predicted label ##ARI=0.9510208
