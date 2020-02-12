#written by Andre Andrade

PCA_ENS_TMLA<-function(BRICK){
  if(raster::nlayers(BRICK)==1){
    ens <- BRICK[[1]]
  }else{
    ens <- rasterPCA(BRICK,spca=T,nComp=1,maskCheck=T)
    ens <- ens$map
    enmean <- calc(BRICK,fun=mean,na.rm=T)
    # enmean <- raster::mean(BRICK)
    co <- raster::layerStats(raster::stack(enmean, ens), 'pearson', na.rm=T)
    co <- co[[1]][1,2]
    if(co<0){
      ens <- ens*-1
    }
    ens <- STANDAR(ens)
  }
  return(ens)
}
