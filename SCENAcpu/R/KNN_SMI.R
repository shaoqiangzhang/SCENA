KNN_SMI <- function(W,K=10,t=50,n) {
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
