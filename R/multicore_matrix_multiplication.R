matmultip <- function(X,C,n_cores) {
	library(doParallel)
	library(foreach)

	ndata = nrow(X)

	mt<-function(X,C,rowlabel){
		X[rowlabel,]%*% C
	}
	
	cls <- makeCluster(n_cores)
	registerDoParallel(cls)
	XC <- foreach(rowlabel=1:ndata,.combine='rbind') %dopar% mt(X,C,rowlabel)
	stopCluster(cls)
	return(XC)
}