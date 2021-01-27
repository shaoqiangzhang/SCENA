KNN_SMI <- function(W,K=10,t=50,n) {
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

  nextW=W;
  for (i in 1:t) {
    prew=nextW;
    nextW=newW %*% (nextW) %*% t(newW);
    lastw=nextW;
    if(max(abs(prew-lastw))<(1/n*0.01))
    {
      break;
    }
  }

  W = nextW + diag(nrow(W));
  W = (W + t(W))/2;


  return(W)
}
