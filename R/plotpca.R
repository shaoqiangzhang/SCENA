plotPCA<-function(Express,cc){
  library(ggplot2)
  #library(RColorBrewer)
  #Exp.pca <- princomp(Express, cor=TRUE, scores=TRUE)
  Exp.pca=prcomp(t(Express),center = TRUE,scale. = TRUE)
   
  Clusters=as.factor(cc)
  ggplot(as.data.frame(Exp.pca$x),aes(x=PC1,y=PC2,colour=Clusters)) + 
  geom_point()+theme_bw()
}
