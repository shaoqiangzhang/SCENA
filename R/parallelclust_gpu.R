select_features_GPU = function(X,K=10,T=100,X1=50,X2=100,X3=150,X4=200,X5=250,Express=Express){
  library(SNFtool)
  library(apcluster)
  library(gpuR)
  alpha=0.5
  
  dist2GPU <- function(X,C) {
    ndata = nrow(X)
    ncentres = nrow(C)

    sumsqX = rowSums(X^2)
    sumsqC = rowSums(C^2)

    XC =matrix(0,ndata, ncentres)
    gpuX=gpuMatrix(X ,type = "float")
    gpuC=gpuMatrix(t(C) ,type = "float")
    gpuXC=gpuMatrix(XC ,type = "float")
    gpuXC = 2 * (gpuX %*% gpuC)
    
    for (i in 1:ndata) {
      for (j in 1:ncentres) {
        XC[i,j]=gpuXC[i,j]
      }
    }
  
    res = matrix(rep(sumsqX,times=ncentres),ndata,ncentres) + t(matrix(rep(sumsqC,times=ndata),ncentres,ndata)) - XC
    res[res < 0] = 0
    return(res)
  }##distance between two matrices

  KNN_SMI_GPU <- function(W,K=20,T=20) {

    library(gpuR)

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

    newW =matrix(0,dim(W)[1], dim(W)[2])
    nextW =matrix(0,dim(W)[1], dim(W)[2])
    W = normalize(W);
    W = (W+t(W))/2;
    newW = (.dominateset(W,K))
	## the following is gpu matrix multiplication
	gpuB=gpuMatrix(newW ,type = "float")
	gpuA=gpuMatrix(W,type = "float")
	if(dim(W)[1]<=5000){
		for (i in 1:floor(log2(T))) {
			gpuB=gpuB %*% gpuB
		}
		gpuA = gpuB %*% (gpuA)%*% t(gpuB); 
		
		for (i in 1:nrow(gpuA)) {
			for (j in 1:ncol(gpuA)) {
				nextW[i,j]=gpuA[i,j]
			}
		}
	}else{
		gpuA=(gpuB+t(gpuB))/2
		for (i in 1:nrow(gpuA)) {
			for (j in 1:ncol(gpuA)) {
				nextW[i,j]=gpuA[i,j]
			}
		}##copy gpu to cpu
	}
    W = nextW + diag(nrow(W));
    W = (W + t(W))/2;


    return(W)
  }
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
  num=apply(Express, 1, sd)
  Express=as.array(Express)
  num=as.array(num)
  Total<-cbind(num,Express)
  Total1<-Total[order(Total[,1],decreasing = TRUE),]
  if(X==1){
    alpha=0.5
    Total1=Total1[1:X1,]
    Total12=Total1[,-1]
    Total12=apply(Total12,1,as.numeric)
    Data11 = standardNormalization(Total12)
    Dist11 = dist2GPU(as.matrix(Data11),as.matrix(Data11))
    W11 = affinityMatrix(Dist11, K, alpha)
    W1 = KNN_SMI_GPU(W11, K, T)
    apresla1<-apcluster(negDistMat(r=7),W1)
    s11<-negDistMat(W1,r=7)
    apreslb1<-apcluster(s11)
    g1=length(apreslb1@exemplars)
    group1 = spectralClustering(W1, g1)
    write.table(group1,file="data1.txt",sep=" ",quote=TRUE)
  }
  if(X==2){
    alpha=0.5
    Total1=Total1[1:X2,]
    Total22=Total1[,-1]
    Total22=apply(Total22,1,as.numeric)
    Data21 = standardNormalization(Total22)
    Dist21 = dist2GPU(as.matrix(Data21),as.matrix(Data21))
    W21 = affinityMatrix(Dist21, K, alpha)
    W2 = KNN_SMI_GPU(W21, K, T)
    apresla2<-apcluster(negDistMat(r=7),W2)
    s21<-negDistMat(W2,r=7)
    apreslb2<-apcluster(s21)
    g2=length(apreslb2@exemplars)
    group2 = spectralClustering(W2, g2)
    write.table(group2,file="data2.txt",sep=" ",quote=TRUE)
  }
  if(X==3){
    alpha=0.5
    Total1=Total1[1:X3,]
    Total32=Total1[,-1]
    Total32=apply(Total32,1,as.numeric)
    Data31 = standardNormalization(Total32)
    Dist31 = dist2GPU(as.matrix(Data31),as.matrix(Data31))
    W31 = affinityMatrix(Dist31, K, alpha)
    W3 = KNN_SMI_GPU(W31, K, T)
    apresla3<-apcluster(negDistMat(r=7),W3)
    s31<-negDistMat(W3,r=7)
    apreslb3<-apcluster(s31)
    g3=length(apreslb3@exemplars)
    group3 = spectralClustering(W3, g3)
    write.table(group3,file="data3.txt",sep=" ",quote=TRUE)
  }
  if(X==4){
    alpha=0.5
    Total1=Total1[1:X4,]
    Total42=Total1[,-1]
    Total42=apply(Total42,1,as.numeric)
    Data41 = standardNormalization(Total42)
    Dist41 = dist2GPU(as.matrix(Data41),as.matrix(Data41))
    W41 = affinityMatrix(Dist41, K, alpha)
    W4 = KNN_SMI_GPU(W41, K, T)
    apresla4<-apcluster(negDistMat(r=7),W4)
    s41<-negDistMat(W4,r=7)
    apreslb4<-apcluster(s41)
    g4=length(apreslb4@exemplars)
    group4 = spectralClustering(W4, g4)
    write.table(group4,file="data4.txt",sep=" ",quote=TRUE)
  }
  if(X==5){
    alpha=0.5
    Total1=Total1[1:X5,]
    Total52=Total1[,-1]
    Total52=apply(Total52,1,as.numeric)
    Data51 = standardNormalization(Total52)
    Dist51 = dist2GPU(as.matrix(Data51),as.matrix(Data51))
    W51 = affinityMatrix(Dist51, K, alpha)
    W5 = KNN_SMI_GPU(W51, K, T)
    apresla5<-apcluster(negDistMat(r=7),W5)
    s51<-negDistMat(W5,r=7)
    apreslb5<-apcluster(s51)
    g5=length(apreslb5@exemplars)
    group5 = spectralClustering(W5, g5)
    write.table(group5,file="data5.txt",sep=" ",quote=TRUE)
  }
}
