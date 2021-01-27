# SCENA (Single-cell Clustering by Enhancing Network Affinity)

SCENA is a R package for unsupervised clustering of single cell RNA-Seq data

Version: 1.0.1

Depends: R(>4.0)

Import packages: parallel, SNFtool, gpuR, apcluster, mclust

Citation: Consensus Clustering of Single-cell RNA-seq Data by Enhancing Network Affinity, to be published. 

# 1. Installation
## 1.1 Start by installing the necessary packages  
```
install.packages("parallel")
install.packages("SNFtool")
install.packages("apcluster")
install.packages("mclust") # this package is used to compute ARI
```
## 1.2 Install the SCENA package
```
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")
```
## 1.3 If you computer supports GPU programming (optional)
```
install.packages("gpuR") ##**the GPU computing package
```
#### ***Note: the package "gpuR" was built on Linux x86_64 (https://www.rdocumentation.org/packages/gpuR), and cannot be installed on a Windows system.*

# 2. Usage examples
## 2.1 An example using SCENA CPU version (dataset: Biase)
First, we load the packages
```
pkgs<-c('SNFtool','apcluster','mclust','parallel','SCENAcpu')
lapply(pkgs,library,character.only=TRUE)
```
Second, we load the dataset (rows are genes and columns are cells)

If the dataset is a txt file:
```
Express=read.table("Biase3celltypes.txt",header = T,row.names = 1)
```

If the dataset is a csv file:

```
Express=read.csv("Biase3celltypes.csv",header = T,row.names = 1)
```

If the dataset is a rds file:

```
library(SingleCellExperiment)
biase<-readRDS("E:/biase.rds")
Express=biase@assays$data$normcounts
```
Third, preprocess the input data
```
Express=datapreprocess(Express,lognum = 1)  #log=1 is do log-transformation, log=0 is no log-transformation
```
Fourth, clustering in parallel with 5 cpu cores

```
detectCores()
cl <- makeCluster(5)  # employ 5 cpu cores
clusterExport(cl,"Express",envir = environment())
parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, select_features1) ## parallel clustering with parameter settings
stopCluster(cl)
```



*##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.*

Fifth, do consensus clustering
```
b=consClust()
```

## 2.2 An example using SCENA GPU version


# Contact info
Author: Yaxuan Cui, Shaoqiang Zhang, Yong Chen

Maintainer: Yaxuan Cui, Shaoqiang Zhang (zhangshaoqiang@tjnu.edu.cn)


