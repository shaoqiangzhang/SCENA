####An example for running SCENA on the Patel's dataset (rds file)

requiredPackages = c('parallel','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

##download data
furl<-"https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/patel.rds"
download.file(furl,destfile="./patel.rds")

#read rds file
library(SingleCellExperiment)
Express=readRDS("patel.rds")
Express=Express@assays$data$logcounts

##data preprocessing
Express=datapreprocess(Express,log=F) #log=F is no log-transformation


## do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=50,X2=100,X3=150,X4=200,X5=250, Express=Express,select_features)
stopCluster(cl)

##do consensus clustering
cc=consClust(9) #set the number of clusters = 9
#cc=consClust() #no parameters if using the predicted number of clusters

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=substring(colnames(Express),1,5) ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.8164146 

