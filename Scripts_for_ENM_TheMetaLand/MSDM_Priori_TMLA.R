# Escrito por Poliana Mendes, Santiago Velazco e Andre Andrade

MSDM_Priori_TMLA <- function(Species,
                             var,
                             MSDM,
                             DirMSDM){
  

  #Metodo 1: LatLong----
  if(MSDM=="LatLong"){
    
    DirPRI<-"LatLong"
    if (file.exists(file.path(DirMSDM,DirPRI))){
      DirPRI<-file.path(DirMSDM,DirPRI)
    } else {
      dir.create(file.path(DirMSDM,DirPRI))
      DirPRI<-file.path(DirMSDM,DirPRI)
    }

    if(length(list.files(DirPRI))!=0){
      print("MSDM found! Using already created MSDM")
    }else{
    
      #Carregar Mascara
      rsd<-as.data.frame(var,xy=T,centroids=T)
      long<-rsd
      lat<-rsd
      lins<-which(is.na(long[,3])==F)
      
      #Criar Raster de Longitude
      long[lins,3]<-long[lins,1]
      gridded(long)<-~x+y
      long<-raster(long)
      
      #Criar Raster de Latitude
      lat[lins,3]<-lat[lins,2]
      gridded(lat)<-~x+y
      lat<-raster(lat)
      
      #Create LatLong Stack
      envM <- stack(long,lat)
      names(envM) <- c("Long","Lat")
      rm(lat,long)
      
      writeRaster(envM,file.path(DirPRI,names(envM)),format="GTiff",bylayer=T,overwrite=T,NAflag=-9999)
    }
  }
  
  #Metodo 2- Distancias minimas----
  if(MSDM=="Min"){
  
    DirPRI<-"Min"
    if (file.exists(file.path(DirMSDM,DirPRI))){
      DirPRI<-file.path(DirMSDM,DirPRI)
    } else {
      dir.create(file.path(DirMSDM,DirPRI))
      DirPRI<-file.path(DirMSDM,DirPRI)
    }
    
    if(length(list.files(DirPRI))==length(Species)){
      print("MSDM found! Using already created MSDM")
    }else{
    
      spi<-as(var,'SpatialPixels')@coords
      r <- lapply(Species, function(x) rasterize(x,var,field=1))
      r<-lapply(r, function(x) as(x,'SpatialPixels')@coords)
      distr <- lapply(r, function(x) dist2(spi,x,method = 'euclidean',p=2))
      distr<-lapply(distr,function(x) apply(x, 1,min))
      distr<-lapply(distr, function(x) (x-min(x))/(max(x)-min(x)))
      envM <- list()
      for (b in 1:length(Species)){
        msk<-var[[1]]
        msk[!is.na(msk[])]<-distr[[b]]
        envM[[b]] <- msk
      }
      envM <- stack(envM)
      names(envM) <- names(Species)
      rm(msk)
      if(nlayers(envM)>1){
        writeRaster(envM,file.path(DirPRI,names(Species)),format="GTiff",bylayer=T,overwrite=T,NAflag=-9999)
      }else{
        writeRaster(envM,file.path(DirPRI,names(Species)),format="GTiff",overwrite=T,NAflag=-9999)
      }
    }
  }
  
  #Metodo 3: Distancia cumulativa----
  if(MSDM=="Cum"){  
    
    DirPRI<-"Cum"
    if (file.exists(file.path(DirMSDM,DirPRI))){
      DirPRI<-file.path(DirMSDM,DirPRI)
    } else {
      dir.create(file.path(DirMSDM,DirPRI))
      DirPRI<-file.path(DirMSDM,DirPRI)
    }
    
    if(length(list.files(DirPRI))==length(Species)){
      print("MSDM found! Using already created MSDM")
    }else{
      
      spi<-as(var,'SpatialPixels')@coords
      r <- lapply(Species, function(x) rasterize(x,var,field=1))
      r<-lapply(r, function(x) as(x,'SpatialPixels')@coords)
      distr <- lapply(r, function(x) dist2(spi,x,method = 'euclidean',p=2))
      distr<-lapply(distr, function(x) x+1)
      distr<-lapply(distr, function (x) 1/(1/x^2))
      distr<-lapply(distr, function(x) apply(x,1,sum))
      distr<-lapply(distr, function(x) (x-min(x))/(max(x)-min(x)))
      envM <- list()
      for (b in 1:length(Species)){
        spdist<-var[[1]]
        spdist[!is.na(spdist[])]<-distr[[b]]
        envM[[b]] <- spdist
      }
      envM <- stack(envM)
      names(envM) <- names(Species)
      rm(spdist)
      if(nlayers(envM)>1){
        writeRaster(envM,file.path(DirPRI,names(Species)),format="GTiff",bylayer=T,overwrite=T,NAflag=-9999)
      }else{
        writeRaster(envM,file.path(DirPRI,names(Species)),format="GTiff",overwrite=T,NAflag=-9999)
      }
    }
  }
  
  #Metodo 4: Kernel-Gaussiano----
  if(MSDM=="Kern"){
    
    DirPRI<-"Kern"
    if (file.exists(file.path(DirMSDM,DirPRI))){
      DirPRI<-file.path(DirMSDM,DirPRI)
    } else {
      dir.create(file.path(DirMSDM,DirPRI))
      DirPRI<-file.path(DirMSDM,DirPRI)
    }
    
    if(length(list.files(DirPRI))==length(Species)){
      print("MSDM found! Using already created MSDM")
    }else{
      
      spi<-as(var,'SpatialPixels')@coords
      r <- lapply(Species, function(x) rasterize(x,var,field=1))
      r<-lapply(r, function(x) as(x,'SpatialPixels')@coords)
      distr <- lapply(r, function(x) dist2(spi,x,method = 'euclidean',p=2))
      distp<-lapply(r, function(x) dist2(x,x,method='euclidean',p=2))
      distp1<-lapply(distp, function(x) matrix(0,nrow(x),1))
      envM <- list()
      for (b in 1:length(Species)){
        distsp <- distp[[b]]
        for (c in 1:nrow(distsp)) {
          vec<-distsp[c,]
          distp1[[b]][c]<-min(vec[vec!=min(vec)])
        }
      sd_graus<-max(distp1[[b]])
      distr2<-distr[[b]]
      distr2<-(1/sqrt(2*pi*sd_graus)*exp(-1*(distr[[b]]/(2*sd_graus^2))))
      distr2<-apply(distr2,1,sum)
      spdist<-var[[1]]
      spdist[!is.na(spdist[])]<-distr2
      envM[[b]] <- spdist
      }
      envM <- stack(envM)
      names(envM) <- names(Species)
      rm(spdist)
      if(nlayers(envM)>1){
        writeRaster(envM,file.path(DirPRI,names(Species)),format="GTiff",bylayer=T,overwrite=T,NAflag=-9999)
      }else{
        writeRaster(envM,file.path(DirPRI,names(Species)),format="GTiff",overwrite=T,NAflag=-9999)
      }
    }
  }
  return(DirPRI)
}#fecha a funcao
