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

##data preprocessing
Express=Express[,2:91] #select 90 cells,you can omit it if you use full data
Express=datapreprocess(Express,log=T) #"log=T" is to do log-transformation, "log=F" is no log-transformation

##do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=50,X2=100,X3=150,X4=200,X5=250, Express=Express,select_features)
##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.
stopCluster(cl)

##do consensus clustering
cc=consClust() #no parameters if using the predicted number of clusters
#cc=consClust(6) #set the number of clusters = 6

##plot scatter graph with PCA
plotPCA(Express,cc) 

##compute ARI as follows:
library(mclust)
presetlabel=substring(colnames(Express),1,4) ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=0.8880762
