Env.PCA.Futuro<-function(Dir,DirP,dirf_out){
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
  
  setwd(dirp_out)
  
  #Salvar os coeficientes
  Coef<-DPca$rotation
  #write.table(Coef,'result_coef.txt',row.names=T,sep='\t')
  
  #Porcentagem da variacaoo explicada
  NEixos<-length(summary(DPca)$importance[3,])
  CumVar<-summary(DPca)$importance[3,]
  VarEx<-data.frame(CumVar)
  write.table(variacao.explicada,'result_latent.txt',row.names=T,sep='\t')      

  #Recuperar os loadings e transformar em um data.frame
  Eix<-as.data.frame(DPca$x)
  EixXY<-cbind(DF[,(1:2)],Eix)
  gridded(EixXY)<- ~x+y
  environmental<-stack(EixXY)
  writeRaster(environmental,"PCA.asc",bylayer=T,suffix='numbers',format="ascii")
  
  #2.Projetar a PCA para variaveis do futuro
  
  print("Select the folder with ASC variables for projection:")
  DirP<-choose.dir()
  Proj<-stack(file.path(DirP,list.files(DirP,pattern='.asc')))

  ProjE<-rasterToPoints(Proj)
  ProjE<-na.omit(ProjE)
  ProjER<-ProjE[,-c(1:2)]
  
  project<-stack()
  nomes.PC<-NULL
  scale<-sweep(ProjER,2,means)
  scale<-scale %*% diag(1/stds)
  
  for (x in 1:ncol(project.env.raw)){
    coeficientes.pc<-coeficientes[,x]
    names(coeficientes.pc)<-NULL
    PC<-project.env.raw %*% diag(coeficientes.pc)
    PC<-rowSums(PC)
    PC<-cbind(project.env[,1:2],PC)
    colnames(PC)<-c("x","y","PC")
    PC<-as.data.frame(PC)
    gridded(PC)<- ~x+y
    PC<-raster(PC)
    project<-stack(project,PC)
  }
  names(project)<-names(environmental)
  setwd(dirf_out)
  writeRaster(project,"PCA.asc",bylayer=T,suffix='numbers',format="ascii")
}