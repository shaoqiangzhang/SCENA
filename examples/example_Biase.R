###An example for running SCENA on the Biase's dataset (download from NCBI GEO)

#install.packages("parallel")
#install.packages("SNFtool")
#install.packages("apcluster")
#install.packages("mclust")
#install.packages("devtools")
#devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

##download data from GEO and read file
furl<-"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57249&format=file&file=GSE57249%5Ffpkm%2Etxt%2Egz"
download.file(furl,destfile="./GSE57249_fpkm.txt.gz")
Express=read.table(gzfile("GSE57249_fpkm.txt.gz"),header = T,row.names = 1)

##read scRNA-seq expression file if you have downloaded and unzipped the file in advance
#Express=read.table("./GSE57249_fpkm.txt",header = T,row.names = 1)

##data preprocessing
Express=datapreprocess(Express,log=T)  #log=T is to do log-transformation, log=F is no log-transformation

## do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, Express=Express,select_features)
##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.
stopCluster(cl)

##do consensus clustering
cc=consClust() #no parameters if using the predicted number of clusters
#cc=consClust(3) #set the number of clusters = 3

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=rep(c(1:3),c(9,20,27)) ## preset 3 cell types each containing 9 cells, 20 cells, 27 cells,respectively.
adjustedRandIndex(presetlabel,as.vector(cc)) ## 'cc' is predicted label

