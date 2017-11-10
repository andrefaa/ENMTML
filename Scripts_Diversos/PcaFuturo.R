PCAFuturo<-function(Dir,DirP){
  #Parametros
    #dir_in: diretorio com as variaveis ambientais do presente
    #dirp_out: diretorio onde ser? salvo o PCA(presente)
  #dirf_out: diretorio onde ser? salvo o PCA(futuro)
  setwd(Dir)
  
  library(raster)
  library(SDMTools)
  
  #1.Realizar a PCA das variaveis no presente
  env<-stack(file.path(Dir,list.files(Dir,pattern='.asc')))
  
  DF<-rasterToPoints(env)
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
  write.table(Coef,'result_coef.txt',row.names=T,sep='\t')
  
  #Porcentagem da variacaoo explicada
  NEixos<-length(summary(DPca)$importance[3,])
  CumVar<-summary(DPca)$importance[3,]
  VarEx<-data.frame(CumVar)
  write.table(VarEx,'result_latent.txt',row.names=T,sep='\t')

  #Recuperar os loadings e transformar em um data.frame
  Eix<-as.data.frame(DPca$x)
  EixXY<-cbind(DF[,(1:2)],Eix)
  gridded(EixXY)<- ~x+y
  PCA<-stack(EixXY)
  PCA.95 <- PCA[[1:(sum(VarEx<=0.95)+1)]]
  writeRaster(PCA.95,paste(Dir,"PC.tif",sep="/"),bylayer=T,suffix='numbers',format="GTiff",overwrite=T)
  
  #2.Projetar a PCA para variaveis do futuro
  
  Proj<-stack(file.path(DirP,list.files(DirP,pattern='.asc')))

  ProjE<-rasterToPoints(Proj)
  ProjE<-na.omit(ProjE)
  ProjER<-ProjE[,-c(1:2)]
  
  nomes.PC<-NULL
  scale<-sweep(ProjER,2,means)
  scale<-scale %*% diag(1/stds)
  PCAFut<-NULL
  
  for (x in 1:ncol(ProjER)){
    CoefPC<-as.numeric(Coef[,x])
    PC<-scale %*% CoefPC
    PCAFut <- cbind(PCAFut,PC)
  }
  colnames(PCAFut) <- names(environmental)
  PCAFut <- data.frame(cbind(ProjE[,(1:2)],PCAFut))
  gridded(PCAFut)<- ~x+y
  PCAFut<-stack(PCAFut)
  names(PCAFut) <- names(PCA)
  PCAFut.95 <- PCAFut[[1:nlayers(PCA.95)]]
  writeRaster(PCAFut.95,paste(DirP,"PC.tif",sep="/"),bylayer=T,suffix='numbers',format="GTiff",overwrite=T)
}

DirP <- "C:\\OneDrive\\R\\Codigos\\ENM_TheMetaLand\\TestePCAFuturo\\futuro\\mr85"

PCAFuturo(Dir,DirP)
