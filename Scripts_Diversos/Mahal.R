# Mahalanobis algorithm
require(raster)
# var=Variable stack
# species=species records
MAHAL<-function(var, species, probability=TRUE){
  envir<-extract(var, species)
  # transformamos os raster em dataframe
  RtoP<-rasterToPoints(var)
  # media das variávies ambientais para cada ponto
  ColMeans<-colMeans(envir)
  # Matriz de variâcia-covar
  varcov <- t(as.matrix(envir)) %*% as.matrix(envir)/nrow(envir)
  # Mahalanobis distance
  map <- mahalanobis(as.matrix(RtoP[,-c(1:2)]), ColMeans, varcov)
  # Probabilidade
  if(probability==TRUE){map <- 1 - pchisq(map, ncol(envir))
  map2<-var[[1]]
  map2[!is.na(map2[])]<-map
  }
  map2<-var[[1]]
  map2[!is.na(map2[])]<-map
  return(map2)
}