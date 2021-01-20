# SCENA
Consensus Clustering of Single-cell RNA-seq Data by Enhancing Network Affinity

## An example using SCENA CPU version (data set: Biase)

```
#The program need to install parallel,SNFtool,apcluster package.
install.packages(parallel)
install.packages(SNFtool)
install.packages(apcluster)
install.packages(mclust)
library(SCENAcpu)
#Import the dataset Express=read.table("E:/Biase3celltypes.txt",header = T, comment.char='!',stringsAsFactors = FALSE,quote = "",sep='\t')
Express=Express[,-1]
len=ncol(Express)
Express=Express[apply(Express,1,function(x) sum(x>1)>len*0.05),]
Express=Express[apply(Express,1,function(x) sum(x>1)<len*0.95),]
Express=apply(Express,2,as.numeric)
#Distribute the number of cores in a computer
library(parallel)
detectCores()
cl <- makeCluster(5)
#Load the file
library(SNFtool)
library(apcluster)
alpha=0.5
clusterExport(cl,"KNN_SMI",envir = environment())
clusterExport(cl,"Express",envir = environment())
clusterExport(cl,"alpha",envir = environment())
#Parallel computation.You can change the parameters
parLapply(cl, c(1,2,3,4,5),K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, select_features1)
stopCluster(cl)
#Read the resulting group file

group1=read.table("./data1.txt",header = T,quote = "",sep=' ')
group2=read.table("./data2.txt",header = T,quote = "",sep=' ')
group3=read.table("./data3.txt",header = T,quote = "",sep=' ')
group4=read.table("./data4.txt",header = T,quote = "",sep=' ')
group5=read.table("./data5.txt",header = T,quote = "",sep=' ')
group1=t(group1)
group2=t(group2)
group3=t(group3)
group4=t(group4)
group5=t(group5)
a1=max(group1)
a2=max(group2)
a3=max(group3)
a4=max(group4)
a5=max(group5)
c1=c(a1,a2,a3,a4,a5)
#Merge the results
tt=result(group1,group2,group3,group4,group5)
#Spectral clustering
#Selection of digital
#Numberis the median in max(group1),max(group2),max(group3),max(group4),max(group5)
res=spectralClustering(tt, c1)
```

## An example using SCENA CPU version (data set: Kolodziejczyk)

```
install.packages(parallel)
install.packages(SNFtool)
install.packages(apcluster)
install.packages(mclust)
library(SCENAcpu)
#Import the dataset
Express=read.table("E:/kolodziejczyk_0.txt",header = T,comment.char='!',stringsAsFactors = FALSE,quote = "",sep='\t')
Express=Express[,-1]
Express=Express[apply(Express,1,function(x) mean(x)>10),]
len=ncol(Express)
Express=Express[apply(Express,1,function(x) sum(x>1)>len0.05),]
Express=Express[apply(Express,1,function(x) sum(x>1)<len0.95),]
dim(Express)
Express=apply(Express,2,as.numeric)
Express=log2(Express+1)
#Distribute the number of cores in a computer
library(parallel)
detectCores()
cl <- makeCluster(5)
#Load the file
library(SNFtool)
library(apcluster)
alpha=0.5
clusterExport(cl,"KNN_SMI",envir = environment())
clusterExport(cl,"Express",envir = environment())
clusterExport(cl,"alpha",envir = environment())
#Parallel computation.You can change the parameters
parLapply(cl, c(1,2,3,4,5),K=20,T=100,X1=50,X2=100,X3=150,X4=200,X5=250, select_features1)
stopCluster(cl)
#Read the resulting group file
group1=read.table("./data1.txt",header = T,quote = "",sep=' ')
group2=read.table("./data2.txt",header = T,quote = "",sep=' ')
group3=read.table("./data3.txt",header = T,quote = "",sep=' ')
group4=read.table("./data4.txt",header = T,quote = "",sep=' ')
group5=read.table("./data5.txt",header = T,quote = "",sep=' ')
group1=t(group1)
group2=t(group2)
group3=t(group3)
group4=t(group4)
group5=t(group5)
a1=max(group1)
a2=max(group2)
a3=max(group3)
a4=max(group4)
a5=max(group5)
c1=c(a1,a2,a3,a4,a5)
#Merge the results
tt=result(group1,group2,group3,group4,group5)
#Spectral clustering
#Selection of digital
#Numberis the median in max(group1),max(group2),max(group3),max(group4),max(group5)
res=spectralClustering(tt, c1)
library(mclust)
#Compute ARI
a<-c(1,1)
for (i in 1:295) {
a[i]=1
}
for (i in 296:454) {
a[i]=2
}
for (i in 455:704) {
a[i]=3
}
a=as.vector(a)
b=as.vector(res)
adjustedRandIndex(a,b)
```
