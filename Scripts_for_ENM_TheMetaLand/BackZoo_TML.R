BackZoo_TML <- function(dir.in,occ,dir.out){
  
  require(raster)
  
  #Import Zooregions ASCs
  files <- list.files(dir.in,pattern = ".asc")
  zoo <- lapply(paste(dir.in,files,sep="/"),raster)
  
  #Get Species Zoogeographical Locations
  occ.spxy <- occ[,2:3]
  ext <- lapply(zoo,extract,occ.spxy) 
  sp.zoo <- zoo[which(lapply(ext,sum)!=0)]
  if(length(sp.zoo)>1){
    sp.zoo$fun <- max
    sp.zoo <- do.call(mosaic,sp.zoo)
  }else{
    sp.zoo <-sp.zoo[[1]]
  }
  writeRaster(sp.zoo,paste(dir.out,paste(as.character(unique(occ[,1])),".asc",sep=""),sep="/"),format="ascii",overwrite=T)
  return(sp.zoo)
}

 