result<-function(group1,group2,group3,group4,group5){
  collect <- matrix(rep(0,5*ncol(Express)),5,ncol(Express))
  collect[1,]=group1
  collect[2,]=group2
  collect[3,]=group3
  collect[4,]=group4
  collect[5,]=group5
  #View(collect)


  tran <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y1=collect[1,]
  x1=t(y1)

  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y1[j]==x1[i])
        tran[i,j]=1
      else
        tran[i,j]=0
    }
  }





  tran2 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y2=collect[2,]
  x2=t(y2)

  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y2[j]==x2[i])
        tran2[i,j]=1
      else
        tran2[i,j]=0
    }
  }




  tran3 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y3=collect[3,]
  x3=t(y3)

  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y3[j]==x3[i])
        tran3[i,j]=1
      else
        tran3[i,j]=0
    }
  }

  tran4 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y4=collect[4,]
  x4=t(y4)

  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y4[j]==x4[i])
        tran4[i,j]=1
      else
        tran4[i,j]=0
    }
  }


  tran5 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y5=collect[5,]
  x5=t(y5)

  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y5[j]==x5[i])
        tran5[i,j]=1
      else
        tran5[i,j]=0
    }
  }

  trantotal=tran+tran2+tran3+tran4+tran5
  #View(trantotal)
  for (i in 1:ncol(collect)) {
    for (j in 1:ncol(collect)) {
      if(trantotal[i,j]>2){
        trantotal[i,j]=1
      }
      else{
        trantotal[i,j]=0
      }
    }
  }
  return(trantotal)
}
