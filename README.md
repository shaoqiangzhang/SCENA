# SCENA (Single-cell Clustering by Enhancing Network Affinity)

SCENA is a R package for unsupervised clustering of single cell RNA-Seq data

Version: 1.0.1

Depends: R(>3.6)

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
install.packages("gpuR") ## see Note1 
```
*##__Note1__: the package "gpuR" was built on Linux x86_64 (https://www.rdocumentation.org/packages/gpuR), and cannot be installed on a Windows system.*

# 2. Usage examples
**There are several detailed examples in the "examples" folder. You can directly run these examples in R.**
## An example for the Biase's dataset is described as follows.

**First**, load the package, download the dataset (rows are genes and columns are cells) from GEO website, and read file.

```
library(SCENA)
furl<-"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57249&format=file&file=GSE57249%5Ffpkm%2Etxt%2Egz"
download.file(furl,destfile="./GSE57249_fpkm.txt.gz")
Express=read.table(gzfile("GSE57249_fpkm.txt.gz"),header = T,row.names = 1)
```
If you have downloaded and unzipped the compressed file in advance, you can read the txt file directly: 
```
library(SCENA)
Express=read.table("./GSE57249_fpkm.txt",header = T,row.names = 1)
```
**Second**, preprocess the input data as follows.
```
Express=datapreprocess(Express,log=T)  #log=T is to do log-transformation, log=F is no log-transformation
```
**Third**, do clustering in parallel with 5 CPU cores as follows. 

```
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, Express=Express,select_features)##see Note2
stopCluster(cl)
```

*##__Note2__: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.*

Alternatively, if your computer supports GPU computing, you can do clustering in parallel with CPU+GPU as follows.

```
library(gpuR)
source('./ApSpe_GPU.R')## see Note 3
cl <- makeCluster(5)
parLapply(cl,1:5,K=10,T=50,X1=200,X2=400,X3=600,X4=800,X5=1000, Express=Express,select_features_GPU)
stopCluster(cl)
```
*##__Note3__: Because the GPU code cannot be called from the installed SCENA package directly, please copy it to your working path and run it using ’source'.*

**Fourth**, do consensus clustering as follows. 
```
cc=consClust() #no parameters if using the predicted number of clusters
```
or
```
cc=consClust(3) #set the number of clusters = 3
```
### Further analysis
#### Plot a scatter diagram with PCA as follows.
```
plotPCA(Express,cc) #  'cc' is label of the predicted clusters
```
#### Compute Adjusted Rand Index (ARI) between preset and predicted cell types as follows.
```
library(mclust)
presetlabel=rep(c(1:3),c(9,20,27)) ## preset 3 cell types each containing 9 cells, 20 cells, 27 cells,respectively.
adjustedRandIndex(presetlabel,as.vector(cc)) ## 'cc' is predicted label
```


# Contact us
Authors: Yaxuan Cui, Shaoqiang Zhang, Yong Chen

Maintainers: Yaxuan Cui, Shaoqiang Zhang (zhangshaoqiang@tjnu.edu.cn)
