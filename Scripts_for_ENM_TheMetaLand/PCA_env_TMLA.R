PCA_env_TMLA<-function(env,Dir){
  #Create PCA Folder
  DirPCA <- file.path(Dir,"PCA")
  dir.create(DirPCA)
  nomes <- names(env)
  
  #Perform PCA
  env <- rasterPCA(env,spca=T)
  coef <- data.frame(cbind(Variables=nomes,env$model$loadings))
  write.table(coef, file.path(DirPCA,"Coefficients.txt"),sep="\t",row.names=F)
  vars<- env$model$sdev^2
  vars<- vars/sum(vars)
  write.table(data.frame(Axis=names(vars),Variance=round(vars,digits=10)), file.path(DirPCA,"CumulativeVariance.txt"),sep="\t",row.names=F)
  env <- env$map[[which(cumsum(vars)<=0.97)]]
  writeRaster(env,file.path(DirPCA,names(env)),format="GTiff",NAflag=-9999,bylayer=T,overwrite=T)
  return(env)
}