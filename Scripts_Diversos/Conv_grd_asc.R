Conv_grd_asc<-function(diretorio){
  setwd(diretorio)
  library(raster)
  environmental<-list.files(pattern='tif')
  raster<-stack(environmental)
  names(raster)<-substr(environmental,1,nchar(environmental)-4)
  writeRaster(raster,paste(names(raster),".asc",sep=""),bylayer=T,format="ascii",NAflag=-9999)
}
  
Conv_grd_asc("C:\\OneDrive\\Doutorado\\Bloco1 - ENM\\Cap2 - Virtual Species ENM\\Modelos\\Broad\\Env\\BIO2")
