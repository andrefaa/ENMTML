RemainingAreasMask <- function(dir.in,dir.out){
  setwd(dir.in)
  sps <- stack(list.files(pattern=".tif"))
  nomes <- substr(names(sps),1,nchar(names(sps))-4)
  #Choose a Mask that will be used as base for projection
  msk <- raster(file.choose())
  crs(sps) <- crs(msk)
  sps <- projectRaster(sps,msk)
  sps <- sps>0
  cels <- NULL
  for(x in 1:nlayers(sps)){
    cels[x] <- sum(na.omit(sps[[x]][]))
  }
  area <- cbind(nomes,cels)
  sps <- sps*msk
  cels <- NULL
  for(x in 1:nlayers(sps)){
    cels[x] <- sum(na.omit(sps[[x]][]))
  }
  area <- cbind(area,cels)
  por <- as.numeric(as.character(area[,3]))/as.numeric(as.character(area[,2]))
  area <- cbind(area,por)
  colnames(area) <- c("Species","Original Area","Remaining Area","")
  setwd(dir.out)
  writeRaster(sps,nomes,format="GTiff",bylayer=T,overwrite=T)
  write.table(area,"RemainingArea.txt",sep="\t",row.names=F)
}

RemainingAreasMask("C:\\OneDrive\\Trabajos\\03 - MataAtlanticaDanira\\MSDM\\Crop",
                   "C:\\OneDrive\\Trabajos\\03 - MataAtlanticaDanira\\Remanescentes")
