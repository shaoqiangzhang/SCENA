# SCENA (Single-cell Clustering by Enhancing Network Affinity)

SCENA is a R package for unsupervised clustering of single cell RNA-Seq data

Version: 1.0.1

Depends: R(>4.0)

Import packages: parallel, SNFtool, gpuR, apcluster, mclust

Citation: Consensus Clustering of Single-cell RNA-seq Data by Enhancing Network Affinity, to be published. 

# 1. Installation
##  Start by installing the necessary packages  
```
install.packages("parallel")
install.packages("SNFtool")
install.packages("apcluster")
install.packages("mclust") # this package is used to compute ARI
```
## Install the SCENA package
```
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA")
```
## If you computer supports GPU programming (optional)
```
install.packages("gpuR") ##**the GPU computing package
```
#### ***Note: the package "gpuR" was built on Linux x86_64 (https://www.rdocumentation.org/packages/gpuR), and cannot be installed on a Windows system.*

# 2. Usage examples
##  An example for the Biase's dataset
First, we load the package
```
library(SCENA)
```
Second, we load the dataset (rows are genes and columns are cells)

If the dataset is a txt file (e.g. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57249):
```
Express=read.table("GSE57249_fpkm.txt",header = T,row.names = 1)
#Express=Express[,1:49] #select a part of cells from the dataset to do clustering
```

If the dataset is a csv file:

```
Express=read.csv("Biase3celltypes.csv",header = T,row.names = 1)
```

If the dataset is a rds file:(e.g.: https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/biase.rds)

```
library(SingleCellExperiment)
biase<-readRDS("biase.rds")
Express=biase@assays$data$normcounts
```
Third, preprocess the input data
```
Express=datapreprocess(Express,lognum = 1)  #log=1 is do log-transformation, log=0 is no log-transformation
```
Fourth, clustering in parallel with 5 GPU cores (or/and GPU cores). In the following commands, only one "parLapply" should be run at the same time.

```
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
clusterExport(cl,"Express",envir = environment())
parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, select_features1) ## the command line only calling CPU cores
#library(gpuR) 
#parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, select_features_gpu) ## the command line calling CPU+GPU cores
stopCluster(cl)
```

*##Note: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.*

Fifth, do consensus clustering
```
b=consClust()
```



# Contact info
Author: Yaxuan Cui, Shaoqiang Zhang, Yong Chen

Maintainer: Yaxuan Cui, Shaoqiang Zhang (zhangshaoqiang@tjnu.edu.cn)


