####An example for running SCENA on the Treutlein's dataset 

##install required packages
requiredPackages = c('parallel','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

##download data
furl<-"https://static-content.springer.com/esm/art%3A10.1038%2Fnature13173/MediaObjects/41586_2014_BFnature13173_MOESM31_ESM.txt"
download.file(furl,destfile="./41586_2014_BFnature13173_MOESM31_ESM.txt")

#read file
expr=read.table("./41586_2014_BFnature13173_MOESM31_ESM.txt")
expr=t(expr)
##data preprocessing
Express=expr[5:23275,2:81]
Express=datapreprocess(Express,log=F) 

## do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=50,X2=100,X3=150,X4=200,X5=250, Express=Express,select_features)
##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.
stopCluster(cl)

##do consensus clustering
cc=consClust(Express,3)

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=expr[4:4,2:81] ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## 'cc' is predicted label # ARI=0.7918102

