PCA_env_TML<-function(env){
  
  pca.raw<-as.matrix(env)
  c<-1:nrow(pca.raw)
  pca.raw.t<-cbind(c,pca.raw)
  c<-pca.raw.t[is.na(pca.raw.t[,2])==F,1]
  rm(pca.raw.t)

  means<-apply(na.omit(pca.raw),2,mean)
  stds<-apply(na.omit(pca.raw),2,sd)
  
  #Scale transform 
  data.scaled <- data.frame(apply(pca.raw,2,scale))
  
  # Realizar a PCA 
  data.pca <- prcomp(~.,data=data.scaled,retx=TRUE,na.action = na.omit)
  
  #Salvar os coeficientes
  coeficientes<-data.pca$rotation
  write.table(coeficientes,'result_coef.txt',row.names=T,sep='\t')
  
  #Porcentagem da variacaoo explicada
  n.eixos<-length(summary(data.pca)$importance[3,])
  cumulat.var<-summary(data.pca)$importance[3,]
  variacao.explicada<-data.frame(cumulat.var)
  write.table(variacao.explicada,'result_latent.txt',row.names=T,sep='\t')      
  
  #Eixos que explicam 95% da varia??o
  var.95<-cumulat.var<=0.96
  
  #Recuperar os loadings e transformar em um data.frame
  eixos<-as.data.frame(data.pca$x)
  for(x in 1:ncol(eixos)){
    pca.raw[c,x]<-eixos[,x]
  }
  pca.raw<-pca.raw[,var.95]
  xy<-xyFromCell(env[[1]], 1:ncell(env[[1]]) )
  pca.raw<-as.data.frame(cbind(xy,pca.raw))
  gridded(pca.raw)<- ~x+y
  env<-stack(pca.raw)
  return(env)
}