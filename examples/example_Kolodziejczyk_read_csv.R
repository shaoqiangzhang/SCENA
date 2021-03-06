###The example for running SCENA on Kolodziejczyk's dataset 

##install required packages
requiredPackages = c('parallel','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")


library(SCENA)

##download and read file
furl<-"https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/kolodziejczyk/counttable_es.csv"
download.file(furl,destfile="./counttable_es.csv")
Express=read.csv("./counttable_es.csv", sep=" ", header = T,row.names = 1)


##data preprocessing
Express=Express[1:38561,]#select all genes 
Express=Express[apply(Express,1,function(x) mean(x)>10),] #delete low expression value
Express=datapreprocess(Express,log=T)  #log=T is to do log-transformation, log=F is no log-transformation

## do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=20,T=100,X1=50,X2=100,X3=150,X4=200,X5=250, Express=Express,select_features)
##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.
stopCluster(cl)

##do consensus clustering
cc=consClust(Express) #no parameters if using the predicted number of clusters
#cc=consClust(Express,3) #set the number of clusters = 3

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=substring(colnames(Express),9,10) ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## ARI=1
