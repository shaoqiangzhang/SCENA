clust=function(K=10,T=100,features=50, Express=Express){
  library(SNFtool)
  library(apcluster)
  alpha=0.5
  .discretisation <- function(eigenVectors) {

    normalize <- function(x) x / sqrt(sum(x^2))
    eigenVectors = t(apply(eigenVectors,1,normalize))

    n = nrow(eigenVectors)
    k = ncol(eigenVectors)

    R = matrix(0,k,k)
    R[,1] = t(eigenVectors[round(n/2),])

    mini <- function(x) {
      i = which(x == min(x))
      return(i[1])
    }

    c = matrix(0,n,1)
    for (j in 2:k) {
      c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
      i = mini(c)
      R[,j] = t(eigenVectors[i,])
    }

    lastObjectiveValue = 0
    for (i in 1:20) {
      eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)

      svde = svd(t(eigenDiscrete) %*% eigenVectors)
      U = svde[['u']]
      V = svde[['v']]
      S = svde[['d']]

      NcutValue = 2 * (n-sum(S))
      if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
        break

      lastObjectiveValue = NcutValue
      R = V %*% t(U)

    }

    return(list(discrete=eigenDiscrete,continuous =eigenVectors))
  }

  .discretisationEigenVectorData <- function(eigenVector) {

    Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
    maxi <- function(x) {
      i = which(x == max(x))
      return(i[1])
    }
    j = apply(eigenVector,1,maxi)
    Y[cbind(1:nrow(eigenVector),j)] = 1

    return(Y)

  }

  .dominateset <- function(xx,KK=20) {
    ###This function outputs the top KK neighbors.

    zero <- function(x) {
      s = sort(x, index.return=TRUE)
      x[s$ix[1:(length(x)-KK)]] = 0
      return(x)
    }
    normalize <- function(X) X / rowSums(X)
    A = matrix(0,nrow(xx),ncol(xx));
    for(i in 1:nrow(xx)){
      A[i,] = zero(xx[i,]);

    }


    return(normalize(A))
  }
  KNN_SMI_CPU <- function(W,K=10,T=50,n_cores) {
    .discretisation <- function(eigenVectors) {

      normalize <- function(x) x / sqrt(sum(x^2))
      eigenVectors = t(apply(eigenVectors,1,normalize))

      n = nrow(eigenVectors)
      k = ncol(eigenVectors)

      R = matrix(0,k,k)
      R[,1] = t(eigenVectors[round(n/2),])

      mini <- function(x) {
        i = which(x == min(x))
        return(i[1])
      }

      c = matrix(0,n,1)
      for (j in 2:k) {
        c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
        i = mini(c)
        R[,j] = t(eigenVectors[i,])
      }

      lastObjectiveValue = 0
      for (i in 1:20) {
        eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)

        svde = svd(t(eigenDiscrete) %*% eigenVectors)
        U = svde[['u']]
        V = svde[['v']]
        S = svde[['d']]

        NcutValue = 2 * (n-sum(S))
        if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
          break

        lastObjectiveValue = NcutValue
        R = V %*% t(U)

      }

      return(list(discrete=eigenDiscrete,continuous =eigenVectors))
    }

    .discretisationEigenVectorData <- function(eigenVector) {

      Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
      maxi <- function(x) {
        i = which(x == max(x))
        return(i[1])
      }
      j = apply(eigenVector,1,maxi)
      Y[cbind(1:nrow(eigenVector),j)] = 1

      return(Y)

    }

    .dominateset <- function(xx,KK=20) {
      ###This function outputs the top KK neighbors.

      zero <- function(x) {
        s = sort(x, index.return=TRUE)
        x[s$ix[1:(length(x)-KK)]] = 0
        return(x)
      }
      normalize <- function(X) X / rowSums(X)
      A = matrix(0,nrow(xx),ncol(xx));
      for(i in 1:nrow(xx)){
        A[i,] = zero(xx[i,]);

      }


      return(normalize(A))
    }
    normalize <- function(X) X / rowSums(X)
    ##First, normalize the data
    newW =matrix(0,dim(W)[1], dim(W)[2])
    nextW =matrix(0,dim(W)[1], dim(W)[2])
    W = normalize(W);
    W = (W+t(W))/2;
    # Calculate the local transition matrix.
    newW = (.dominateset(W,K))
	if(dim(W)[1] <=10000) {
		for (i in 1:floor(log2(T))) {
		#newW=newW %*% newW;
		newW=matmultip(newW,newW,n_cores)
		}
		#nextW=newW %*% (W) %*% t(newW);
		nextW=matmultip(newW,W,n_cores)
		nextW=matmultip(nextW,t(newW),n_cores)
	}else{
		nextW=(newW + t(newW))/2 ;## enhance W by addition 
	}
	W = nextW + diag(nrow(W));
	W = (W + t(W))/2;
    return(W)
  }

  num=apply(Express, 1, sd)
  Express=as.array(Express)
  num=as.array(num)
  Total<-cbind(num,Express)
  Total1<-Total[order(Total[,1],decreasing = TRUE),]
  #n=ncol(Express)
  
  library(doParallel)
  library(foreach)
  
  variablecores=detectCores() - 1;
  
  if(variablecores>1){
    n_cores=variablecores
  }else{
    n_cores=1
  }

    #alpha=0.5
    Total1=Total1[1:features,]
    #Selecting genes
    Total12=Total1[,-1]
    Total12=apply(Total12,1,as.numeric)
    Data11 = standardNormalization(Total12)
    Dist11 = dist3(as.matrix(Data11),as.matrix(Data11),n_cores)
    W11 = affinityMatrix(Dist11, K, alpha)
    W1 = KNN_SMI_CPU(W11, K, T, n_cores)
    #function reference KNN_SMI
    apresla1<-apcluster(negDistMat(r=7),W1)
    s11<-negDistMat(W1,r=7)
    apreslb1<-apcluster(s11)
    #AP forecast
    g1=length(apreslb1@exemplars)
    group1 = spectralClustering(W1, g1)
    #write.table(group1,file="clust.txt",sep=" ",quote=TRUE)
	
	return(group1)


}
