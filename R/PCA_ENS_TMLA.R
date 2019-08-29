#written by Andre Andrade

PCA_ENS_TMLA<-function(BRICK){
  if(nlayers(BRICK)==1){
    ens <- BRICK
    return(ens)
  }else{
    ens <- rasterPCA(BRICK,spca=T,nComp=1)
    ens <- ens$map
    enmean <- mean(BRICK)
    co <- layerStats(stack(enmean, ens), 'pearson', na.rm=T)[[1]][1,2]
    if(co<0){
      ens <- ens*-1
    }
    ens <- STANDAR(ens)
    return(ens)
  }
}
