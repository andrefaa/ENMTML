PCA_env<-function(env){
  env <- rasterPCA(env,spca=T)
  vars<- env$model$sdev^2
  vars<- vars/sum(vars)
  write.table(cbind(names(env$map),vars),"result_latent.txt",sep="\t",row.names=F)
  load <- env$model$loadings
  write.table(load,"result_coeff.txt",sep="\t")
  env <- env$map[[which(cumsum(vars)<=0.96)]]
  writeRaster(env,names(env),bylayer=T,format="GTiff",overwrite=T)
}
