PCAFuturo<-function(Env,
                    Dir,
                    DirP,
                    Save=''){

  #0.Create PCA Folder
  DirPCA <- file.path(Dir,"PCA")
  dir.create(DirPCA)
  #PCA Tables Folder
  DirPCATab <- file.path(DirPCA,"Tables")
  dir.create(DirPCATab)
  
  #Projection Folders
  
  DirP_PCA <- file.path(dirname(Dir),'Projection_PCA')
  dir.create(DirP_PCA)
  FoldersProj <- list.files(dirname(DirP[[1]]))
  FoldersProj <- as.list(file.path(DirP_PCA,FoldersProj))
  lapply(FoldersProj, function(x) dir.create(x))
  
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
  write.table(Coef2,file.path(DirPCATab,"Coeficient.txt"),sep="\t",row.names=F)
  
  #Cummulative Variance
  NEixos<-length(summary(DPca)$importance[3,])
  CumVar<-summary(DPca)$importance[3,]
  VarEx<-data.frame(CumVar)
  write.table(VarEx,file.path(DirPCATab,"CumulativeVariance.txt"),sep="\t",row.names=F)
  
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
  
  DirP <- as.list(DirP)
  ProjEX <- lapply(DirP, function(x) unique(file_ext(list.files(x))))
  form <- c('bil','asc','txt','tif')
  ProjEX <- unique(unlist(ProjEX)[unlist(ProjEX)%in%form])

  if(any(ProjEX %in% c('asc', 'bil', 'tif'))){
    ProjT<-lapply(DirP, function(x) brick(stack(file.path(x,list.files(x,paste0('\\.',ProjEX,'$'))))))
  }

  if(any(ProjEX %in% 'txt')){
    ProjT <- list()
    for(j in DirP){
      ProjT[[i]]<-read.table(file.path(DirP[[i]],list.files(DirP[[i]],pattern='.txt'),h=T))
      gridded(ProjT[[i]])<- ~x+y
      ProjT[[i]]<-brick(stack(ProjT[[i]]))
    }
  }

  ProjE<-lapply(ProjT, function(x) rasterToPoints(x))
  ProjE<-lapply(ProjE, function(x) na.omit(x))
  ProjER <- lapply(ProjE, function(z) z[,!(colnames(z) %in% c("x", "y"))])
  
  scale<-lapply(ProjER, function(x) sweep(x,2,means))
  scale<-lapply(scale, function(x) x %*% diag(1/stds))
  PCAFut<-lapply(scale, function(x) x %*% Coef)
  PCAFut <- lapply(PCAFut, function(x) data.frame(cbind(ProjE[[1]][,(1:2)],x)))
  PCAFut.95 <- list()
  for(j in 1:length(PCAFut)){
    gridded(PCAFut[[j]])<- ~x+y
    PCAFut[[j]]<-stack(PCAFut[[j]])
    names(PCAFut[[j]]) <- names(PCAPr)
    PCAFut.95[[j]] <- PCAFut[[j]][[1:nlayers(PCA.95)]]
    if(Save=="Y"){
      writeRaster(PCAFut.95[[j]],paste(FoldersProj[[j]],names(PCAFut.95[[j]]),sep="/"),bylayer=T,format="GTiff",overwrite=T)
    }
  }
  return(PCAFut.95)
}
