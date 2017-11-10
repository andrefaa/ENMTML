PCAFuturo<-function(Env,
                    Dir,
                    DirP,
                    Save=''){

  #1.Realizar a PCA das variaveis no presente
  DF<-rasterToPoints(Env)
  DF<-na.omit(DF)
  PcaR<-DF[,-c(1:2)]
  
  means<-colMeans(PcaR)
  stds<-apply(PcaR,2,sd)
  
  #Scale transform 
  DScale <- data.frame(apply(PcaR,2,scale))
  
  # Realizar a PCA 
  DPca <- prcomp(DScale,retx=TRUE,center=F,scale=F)

  #Salvar os coeficientes
  Coef<-DPca$rotation
  
  #Porcentagem da variacaoo explicada
  NEixos<-length(summary(DPca)$importance[3,])
  CumVar<-summary(DPca)$importance[3,]
  VarEx<-data.frame(CumVar)
  
  #Recuperar os loadings e transformar em um data.frame
  Eix<-as.data.frame(DPca$x)
  EixXY<-cbind(DF[,(1:2)],Eix)
  gridded(EixXY)<- ~x+y
  PCAPr<-stack(EixXY)
  PCA.95 <- PCAPr[[1:(sum(VarEx<=0.95)+1)]]
  if(Save=="Y"){
    writeRaster(PCA.95,paste(Dir,names(PCA.95),sep="/"),bylayer=T,format="GTiff",overwrite=T)
  }
  
  #2.Projetar a PCA para variaveis do futuro
  
  ProjEX <- unique(file_ext(list.files(DirP)))
  form <- c('bil','asc','txt','tif')
  ProjEX <- ProjEX[ProjEX%in%form]
  
  if(ProjEX == 'bil'){
    ProjT<-brick(stack(file.path(DirP,list.files(DirP,pattern='.bil'))))
  }
  
  if(ProjEX == 'asc'){
    ProjT<-brick(stack(file.path(DirP,list.files(DirP,pattern='.asc'))))
  }
  
  if(ProjEX == 'txt'){
    ProjT<-read.table(file.path(DirP,list.files(DirP,pattern='.txt'),h=T))
    gridded(ProjT)<- ~x+y
    ProjT<-brick(stack(ProjT))
  }
  
  if(ProjEX == 'tif'){
    ProjT<-brick(stack(file.path(DirP,list.files(DirP,pattern='.tif'))))
  }
  
  ProjE<-rasterToPoints(ProjT)
  ProjE<-na.omit(ProjE)
  ProjER<-ProjE[,-c(1:2)]
  
  scale<-sweep(ProjER,2,means)
  scale<-scale %*% diag(1/stds)
  PCAFut<-scale %*% Coef
  colnames(PCAFut) <- colnames(Coef)
  PCAFut <- data.frame(cbind(ProjE[,(1:2)],PCAFut))
  gridded(PCAFut)<- ~x+y
  PCAFut<-stack(PCAFut)
  names(PCAFut) <- names(PCAPr)
  PCAFut.95 <- PCAFut[[1:nlayers(PCA.95)]]

  if(Save=="Y"){
    writeRaster(PCAFut.95,paste(DirP,names(PCAFut.95),sep="/"),bylayer=T,format="GTiff",overwrite=T)
  }
  return(PCAFut.95)
}