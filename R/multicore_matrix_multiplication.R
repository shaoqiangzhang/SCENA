matmultip <- function(X,C,n_cores) {
	library(doParallel)
	library(foreach)

	ndata = nrow(X)

	mt<-function(X,C,rowlabel){
		X[rowlabel,]%*% C
	}
	
	if((detectCores() - 1)<n_cores){
		n_cores=detectCores() -1
	}
	
	cls <- makeCluster(n_cores)
	registerDoParallel(cls)
	XC <- foreach(rowlabel=1:ndata,.combine='rbind') %dopar% mt(X,C,rowlabel)
	stopCluster(cls)
	return(XC)
}