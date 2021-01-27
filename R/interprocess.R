
datapreprocess<-function(Express,lognum) {
  Express=apply(Express,2,as.numeric)
  len=ncol(Express)
  Express=Express[apply(Express,1,function(x) sum(x>1)>len*0.05),]
  Express=Express[apply(Express,1,function(x) sum(x>1)<len*0.95),]
  if(lognum==1){
    Express=log2(Express+1)
  }
  return(Express)
}


consClust<-function(){
  library(SNFtool)
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
  print("Cluster label:")
  print(res)
  return(res)
}

