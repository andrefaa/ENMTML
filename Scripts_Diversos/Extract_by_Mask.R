Extract_by_Mask("C:\\OneDrive\\ProjetoAves\\ENM\\Result\\Richness&BetaDiv_Coleta","C:\\OneDrive\\ProjetoAves\\ENM\\Result\\Richness&BetaDiv_Coleta\\Antropico")

Extract_by_Mask<-function(dir.in,dir.out){
  #Fuction to crop an asc layer by a mask
  
  library(SDMTools)
  library(raster)
  
  env_prep<-list.files(dir.in,pattern='.asc')
  asc_base<-read.asc(paste(dir.in,env_prep[1],sep="/"))
  bloco<-raster(asc_base)
  camadas<-env_prep[-1]
  
  for (a in 1:length(camadas)){
    asc<-read.asc(paste(dir.in,camadas[a],sep="/"))
    raster<-raster(asc)
    bloco<-stack(bloco,raster)
  }
  names(bloco)<-substr(env_prep,1,nchar(env_prep)-4)
  
  
  print("Select Mask file:")
  msk<-shapefile(file.choose())
  msk <- spTransform(msk, CRSobj=crs(bloco))
  asc.msk<-mask(bloco,msk,inverse=T)
  #cat("Choose output filename: ")
  #filename <- as.character(readLines(n = 1))
  #writeRaster(asc.msk,paste(dir.out,paste(filename,".asc",sep=""),sep="/"),format="ascii")
  writeRaster(asc.msk,paste(dir.out,paste(names(bloco),".asc",sep=""),sep="/"),bylayer=T,format="ascii")
}