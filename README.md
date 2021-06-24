# SCENA (Single-cell Clustering by Enhancing Network Affinity)

SCENA is a R package for unsupervised clustering of single cell RNA-Seq data

Version: 1.0.3

Depends: R(>3.6)

Import packages: parallel, doParallel,foreach, SNFtool, gpuR, apcluster, mclust, ggplot2

Citation: Consensus Clustering of Single-cell RNA-seq Data by Enhancing Network Affinity, Briefings in Bioinformatics 2021, DOI: 10.1093/bib/bbab236 

# 1. Installation
##  Start by installing the necessary packages  
```
requiredPackages = c('parallel','doParallel','foreach','ggplot2','SNFtool',"apcluster","mclust")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
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

## To speed up matrix computing, recommend you install the "OpenBLAS" library on your operating system(Linux as an example)
 ```
 git clone https://github.com/xianyi/OpenBLAS.git
 cd OpenBLAS
 make FC=gfortran
 sudo make PREFIX=/usr/local install
 ```
 then,set the number of threads using environment variables as follows
 ```
 export OPENBLAS_NUM_THREADS=4
 ```

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

Express=Express[,2:91] #select 90 cells. You can omit it if you use full data
```
*##__Note2__: If the file suffix is ".txt" or ".tsv" , please use "read.table" instead of "read.csv".*

**Step 2**: preprocess the input data as follows.
```
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
*For a big dataset(number of cells>1000), if you just need an acceptable result quickly, you can set the number of iterations T to be small (e.g. T=20).*

Alternatively, if your computer supports GPU computing, you can do clustering in parallel with CPU+GPU as follows.

```
library(gpuR)
source('./parallelclust_gpu.R')## see Note 4
cl <- makeCluster(5)
parLapply(cl,1:5,K=10,T=50,X1=50,X2=100,X3=150,X4=200,X5=250,Express=Express,select_features_GPU)
stopCluster(cl)
```
*##__Note4__: Because the GPU code cannot be called from the installed SCENA package directly, please copy it to your working path and run it using ’source'.*

**Step 4**: do consensus clustering as follows. 
```
cc=consClust(Express) 
```
or
```
cc=consClust(Express,6) #set the number of clusters =6
```
### For easy of use, Steps 3 and 4 for CPU computing have been integrated into a new function "scena_cpu" as follows.
```
cc=scena_cpu(Express)  ##use default parameters
```
or
```
cc=scena_cpu(Express,T=20) ## "T" is the the number of matrix iterations.  
```
or
```
cc=scena_cpu(Express,T=20,num=6) ##"num" is the number of clusters
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
Authors: Shaoqiang Zhang, Yaxuan Cui, Yong Chen

Maintainer: Shaoqiang Zhang (zhangshaoqiang@tjnu.edu.cn)
