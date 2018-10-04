PCAFuturo<-function(Env,
                    Dir,
                    DirP,
                    Save=''){

  #0.Create PCA Folder
  DirPCA <- file.path(Dir,"PCA")
  dir.create(DirPCA)
  
  #1.Present PCA
  DF<-rasterToPoints(Env)
  DF<-na.omit(DF)
  PcaR<-DF[,-c(1:2)]
  
  means<-colMeans(PcaR)
  stds<-apply(PcaR,2,sd)
  
  #Scale transform 
  DScale <- data.frame(apply(PcaR,2,scale))
  
  # PCA 
  DPca <- prcomp(DScale,retx=TRUE,center=F,scale=F)

  #Coefficients
  Coef<-DPca$rotation
  Coef2 <- data.frame(cbind(Variable=names(Env),Coef))
  write.table(Coef2,file.path(DirPCA,"Coeficient.txt"),sep="\t",row.names=F)
  
  #Cummulative Variance
  NEixos<-length(summary(DPca)$importance[3,])
  CumVar<-summary(DPca)$importance[3,]
  VarEx<-data.frame(CumVar)
  write.table(VarEx,file.path(DirPCA,"CumulativeVariance.txt"),sep="\t",row.names=F)
  
  #Save PCsthat account for 95% of total variability
  Eix<-as.data.frame(DPca$x)
  EixXY<-cbind(DF[,(1:2)],Eix)
  gridded(EixXY)<- ~x+y
  PCAPr<-stack(EixXY)
  PCA.95 <- PCAPr[[1:(sum(VarEx<=0.95)+1)]]
  if(Save=="Y"){
    writeRaster(PCA.95,paste(DirPCA,names(PCA.95),sep="/"),bylayer=T,format="GTiff",overwrite=T)
  }
  
  #2.Project PCA
  
  ProjEX <- unique(file_ext(list.files(DirP)))
  form <- c('bil','asc','txt','tif')
  ProjEX <- ProjEX[ProjEX%in%form]
  
  if(any(ProjEX == c('asc', 'bil', 'tif'))){
    ProjT<-brick(stack(file.path(DirP,list.files(DirP,paste0('\\.',ProjEX,'$')))))
  }

  if(ProjEX == 'txt'){
    ProjT<-read.table(file.path(DirP,list.files(DirP,pattern='.txt'),h=T))
    gridded(ProjT)<- ~x+y
    ProjT<-brick(stack(ProjT))
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
  DirPCAF <- file.path(DirP,"PCA")
  dir.create(DirPCAF)

  if(Save=="Y"){
    writeRaster(PCAFut.95,paste(DirPCAF,names(PCAFut.95),sep="/"),bylayer=T,format="GTiff",overwrite=T)
  }
  return(PCAFut.95)
}