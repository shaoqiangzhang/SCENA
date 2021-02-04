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

**The "data" folder includes the download URLs of some datasets** 

## An example for Yan's dataset is described as follows.

**Step 1**: load the package, download the dataset (rows are genes and columns are cells), and read file.

```
library(SCENA)
furl<-"https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/yan/nsmb.2660-S2.csv"
download.file(furl,destfile="./nsmb.2660-S2.csv")
Express=read.csv("./nsmb.2660-S2.csv", header = T,row.names = 1) ##see Note 2
```
*##__Note2__: If the file suffix is ".txt" or ".tsv" , please use "read.table" instead of "read.csv".*

**Step 2**: preprocess the input data as follows.
```
Express=Express[,2:91] #select 90 cells. You can omit it if you use full data
Express=datapreprocess(Express,log=T)  #log=T is to do log-transformation, log=F is no log-transformation
```
At the end of the preprocessing step, suggestions for setting parameters are provided.

**Step 3**: do clustering in parallel with 5 CPU cores as follows. 

```
detectCores()
cl <- makeCluster(5)  # call 5 cpu cores
parLapply(cl,1:5,K=10,T=50,X1=50,X2=100,X3=150,X4=200,X5=250,Express=Express,select_features)##see Note3
stopCluster(cl)
```

*##__Note3__: K is the number of K-nearest neighbors; T is the number of matrix iterations, X1~X5 are top number of selected features.*

Alternatively, if your computer supports GPU computing, you can do clustering in parallel with CPU+GPU as follows.

```
library(gpuR)
source('./ApSpe_GPU.R')## see Note 4
cl <- makeCluster(5)
parLapply(cl,1:5,K=10,T=50,X1=50,X2=100,X3=150,X4=200,X5=250,Express=Express,select_features_GPU)
stopCluster(cl)
```
*##__Note4__: Because the GPU code cannot be called from the installed SCENA package directly, please copy it to your working path and run it using ’source'.*

**Step 4**: do consensus clustering as follows. 
```
cc=consClust() #no parameters if using the predicted number of clusters
```
or
```
cc=consClust(6) #set the number of clusters =6
```
### For easy of use, Steps 2~4 have been integrated into a new function "scena_cpu" as follows.
```
Express=Express[,2:91] #select 90 cells. You can omit it if you use full data
cc=scena_cpu(Express,log=T) ## only call CPU.  
```

#### Finally, you can plot a scatter graph with PCA as follows.
```
plotPCA(Express,cc) #  'cc' is label of the predicted clusters
```
#### and compute Adjusted Rand Index (ARI) between preset and predicted cell types as follows.
```
library(mclust)
presetlabel=substring(colnames(Express),1,4) ## read the column names as preset cell types
adjustedRandIndex(presetlabel,as.vector(cc)) ## 'cc' is predicted label
```


# Contact us
Authors: Yaxuan Cui, Shaoqiang Zhang, Yong Chen

Maintainers: Yaxuan Cui, Shaoqiang Zhang (zhangshaoqiang@tjnu.edu.cn)
