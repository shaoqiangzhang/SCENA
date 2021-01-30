###An example for running SCENA on the Goolam's dataset (download from EBI ArrayExpress)
install.packages("parallel")
install.packages("SNFtool")
install.packages("apcluster")
install.packages("mclust")
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")

library(SCENA)

##download data from EBI and read file
furl<-"https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3321/E-MTAB-3321.processed.1.zip"
download.file(furl,destfile="./E-MTAB-3321.processed.1.zip")
Express=read.table(unz("E-MTAB-3321.processed.1.zip",filename="Goolam_et_al_2015_count_table.tsv"),header = T,row.names = 1)

##read scRNA-seq expression file if you have downloaded and unzipped the file in advance
#Express=read.table("./Goolam_et_al_2015_count_table.tsv",header = T,row.names = 1)

##data preprocessing
Express=Express[1:41389,]
Express=datapreprocess(Express,log=T)  #log=T is to do log-transformation, log=F is no log-transformation

## do clustering in parallel with 5 cpu cores
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, Express=Express,select_features)
##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.
stopCluster(cl)

##do consensus clustering
cc=consClust() #no parameters if using the predicted number of clusters
#cc=consClust(6) #set the number of clusters = 6

##plot scatter graph with PCA
plotPCA(Express,cc) #  'cc' is label of the predicted clusters

##compute ARI as follows:
library(mclust)
presetlabel=rep(c(2,4,8,16,32,8,2,4),c(6,40,16,6,6,16,10,24)) ## preset 5 cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## 'cc' is predicted label

