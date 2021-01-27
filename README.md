# SCENA (Single-cell Clustering by Enhancing Network Affinity)

SCENA is a R package for unsupervised clustering of single cell RNA-Seq data

Version: 1.0.1

Depends: R(>3.3)

Import packages: parallel, SNFtool, gpuR, apcluster, mclust

Citation: Consensus Clustering of Single-cell RNA-seq Data by Enhancing Network Affinity, to be published. 

# 1. Installation
## 1.1 Start by installing the necessary packages  
```
install.packages(parallel)
install.packages(SNFtool)
install.packages(apcluster)
install.packages(mclust) 
```
## 1.2 Install the SCENA package

Option1: download the file SCENAcpu_0.1.0.tar.gz, and install it in R

Option2: 
```
install.packages("devtools")
devtools::install_github("shaoqiangzhang/SCENA/SCENAcpu")
```
## 1.3 Note: the GPU version can be directly based on the CPU version (see the example) 

# 2. Usage examples
## 2.1 An example using SCENA CPU version (dataset: Biase)
First, we load the packages
```
library(SNFtool)
library(parallel)
library(apcluster)
library(SCENAcpu)
```
Second, we load the dataset (rows are genes, while columns are cells)

for a txt file:
```
Express=read.table("Biase3celltypes.txt",header = T, comment.char='!',stringsAsFactors = FALSE,quote = "",sep='\t')
```

for a csv file:

```
Express=read.csv
```

for a rds file:

```
Express=
```

## 2.2 An example using SCENA GPU version


# Contact info
Author: Yaxuan Cui, Shaoqiang Zhang

Maintainer: Yaxuan Cui, Shaoqiang Zhang (zhangshaoqiang@tjnu.edu.cn)


