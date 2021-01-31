####An example for running SCENA on the Patel's dataset (download from NCBI GEO)

#install.packages("parallel")
#install.packages("SNFtool")
#install.packages("apcluster")
#install.packages("mclust")
#install.packages("devtools")
#devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

##download data from GEO and read file
furl<-"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57872&format=file&file=GSE57872%5FGBM%5Fdata%5Fmatrix%2Etxt%2Egz"
download.file(furl,destfile="./GSE57872_GBM_data_matrix.txt.gz")
Express=read.table(gzfile("GSE57872_GBM_data_matrix.txt.gz"),header = T,row.names = 1)

##read scRNA-seq expression file if you have downloaded and unzipped the file in advance
#Express=read.table("./GSE57872_GBM_data_matrix.txt",header = T,row.names = 1)

##data preprocessing
Express=Express[,1:430] #select 430 cells
Express=datapreprocess(Express,log=F)


## do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=50,X2=100,X3=150,X4=200,X5=250, Express=Express,select_features)
##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.
stopCluster(cl)

##do consensus clustering
cc=consClust(9) #set the number of clusters = 9
#cc=consClust() #no parameters if using the predicted number of clusters

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=substring(colnames(Express),1,5) ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## 'cc' is predicted label # ARI=0.8164146 or 0.6851536

