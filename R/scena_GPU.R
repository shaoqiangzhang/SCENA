scena_gpu = function(Express=Express,T=T,num=num){
  library(SNFtool)
  library(apcluster)
  library(parallel)
  

  len=ncol(Express)
  featnum=nrow(Express)
  ##Recommend parameter settings in the next step
  if(len<500){
	K=10
	it=50
  }else if(len<1000){
	K=20
	it=100
  }else{
	K=20
	it=round(len/10, digits = 0)
  }
  if(featnum<=10000){
	X=c(50,100,150,200,250)
  }else if(featnum<=12000){
	X=c(50,100,200,400,800)
  }else{
	X=c(200,400,600,800,1000)
  }
  
  if(missing(T)){
	T=it
  }
  
  cl <- makeCluster(5)  # call 5 cpu cores
  #source("ApSpe_GPU.R")
  parLapply(cl,1:5,K=K,T=T,X1=X[1],X2=X[2],X3=X[3],X4=X[4],X5=X[5],Express=Express,select_features_GPU)
  stopCluster(cl)
  
  group1=read.table("./data1.txt",header = T,quote = "",sep=' ')
  group2=read.table("./data2.txt",header = T,quote = "",sep=' ')
  group3=read.table("./data3.txt",header = T,quote = "",sep=' ')
  group4=read.table("./data4.txt",header = T,quote = "",sep=' ')
  group5=read.table("./data5.txt",header = T,quote = "",sep=' ')
  file.remove("./data1.txt")
  file.remove("./data2.txt")
  file.remove("./data3.txt")
  file.remove("./data4.txt")
  file.remove("./data5.txt")
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
  #Merge the results
  tt=result(Express,group1,group2,group3,group4,group5)
  if(missing(num)==TRUE){
    c1=c(a1,a2,a3,a4,a5)
    c1=median(c1)
    res=spectralClustering(tt, c1)
  }
  else
    res=spectralClustering(tt, num)

  cat("Cluster label:\n")
  print(res)
  return(res)
}