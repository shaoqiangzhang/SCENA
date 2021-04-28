
datapreprocess<-function(Express,log=F) {
  library(SNFtool)
  #library(apcluster)
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
  }else if(len>=5000){
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
