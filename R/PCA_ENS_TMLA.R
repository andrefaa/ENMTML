#written by Andre Andrade

PCA_ENS_TMLA<-function(BRICK){
  ens <- rasterPCA(BRICK,spca=T,nComp=1)
  ens <- ens$map
  ens <- STANDAR(ens)
  return(ens)
}
