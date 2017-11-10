PCA_env_TMLA<-function(env){
  env <- rasterPCA(env,spca=T)
  vars<- env$model$sdev^2
  vars<- vars/sum(vars)
  env <- env$map[[which(cumsum(vars)<=0.97)]]
  writeRaster(env,names(env),format="GTiff",NAflag=-9999,bylayer=T,overwrite=T)
  return(env)
}