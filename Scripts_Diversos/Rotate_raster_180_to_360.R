Rotate_180_to_360<-function(dir_input,dir_output){

library(raster)
library(SDMTools)
  
  setwd(dir_input)

  env_prep<-list.files(pattern='.asc')
  asc_base<-read.asc(env_prep[1])
  environmental<-raster(asc_base)
  camadas<-env_prep[-1]
  
  for (a in 1:length(camadas)){
    asc<-read.asc(camadas[a])
    raster<-raster(asc)
    environmental<-stack(environmental,raster)
  }
  names(environmental)<-substr(env_prep,1,nchar(env_prep)-4)  

  e<-extent(0,360,-70,70)
  extent(environmental)<-e  
  rr<-rotate(environmental)
  
  #pacific<-extent(-80,110,-70,70)
  #rr_pacific<-crop(rr,pacific)
  
  setwd(dir_output)

  for(x in 1:nlayers(rr)){
    writeRaster(rr[[x]],paste(names(rr[[x]]),"_Rotated.asc",sep=""),format="ascii")
  }
}

Rotate_180_to_360("D:\\Users\\Andre\\Google Drive\\Mestrado\\Layers\\ETOPO Oceano","C:\\OneDrive\\Poliquetas\\Poliquetas_Pacifico\\Rotated BioOracle")
