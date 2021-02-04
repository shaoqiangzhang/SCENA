
datapreprocess<-function(Express,log=F) {
  library(SNFtool)
  library(apcluster)
  library(parallel)
  
  Express=apply(Express,2,as.numeric)
  len=ncol(Express)
  Express=Express[apply(Express,1,function(x) sum(x>1)>len*0.05),]
  Express=Express[apply(Express,1,function(x) sum(x>1)<len*0.95),]
  if(log==TRUE){
    Express=log2(Express+1)
  }
  featnum=nrow(Express)
  ##Recommend parameter settings in the next step
  if(len<500){
	K=10
	T=50
  }else if(len<1000){
	K=20
	T=100
  }else{
	K=20
	T=round(len/10, digits = 0)
  }
  if(featnum<=10000){
	X=c(50,100,150,200,250)
  }else if(featnum<=12000){
	X=c(50,100,200,400,800)
  }else{
	X=c(200,400,600,800,1000)
  }
  cat("*Note: Recommend parameter settings in the next clustering step:\n K=", K, ",T=",T,",X1~X5=",X) 
   return(Express)
}


consClust<-function(num=num){
  library(SNFtool)
  group1=read.table("./data1.txt",header = T,quote = "",sep=' ')
  group2=read.table("./data2.txt",header = T,quote = "",sep=' ')
  group3=read.table("./data3.txt",header = T,quote = "",sep=' ')
  group4=read.table("./data4.txt",header = T,quote = "",sep=' ')
  group5=read.table("./data5.txt",header = T,quote = "",sep=' ')
 # file.remove("./data1.txt")
  #file.remove("./data2.txt")
  #file.remove("./data3.txt")
  #file.remove("./data4.txt")
  #file.remove("./data5.txt")
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
  tt=result(group1,group2,group3,group4,group5)
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



plotPCA<-function(Express,b){
  Exp.pca <- princomp(Express, cor=TRUE, scores=TRUE)
  Exp.pca=prcomp(t(Express),center = TRUE,scale. = TRUE)
  plot(Exp.pca, type = "l")
  names(Exp.pca)
  summary(Exp.pca)
  #View(summary(Exp.pca))
  #View(Exp.pca$x)
  plot(Exp.pca$x,col=b,xaxt="n", yaxt="n",pch=16)
}

