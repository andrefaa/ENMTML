MarChegarnaPraia <- function(Dir,DirO){
  
  #Dir <- "C:\\OneDrive\\Variaveis\\Marinhas\\OSCAR_Correntes"
  #CarregarPacotes
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  
  ipak(c("raster","sp","dismo","rgdal","maptools","RStoolbox","tools",
         "SDMTools","maps","flexclust","plyr"))
  
  #Definir o diretorio
  setwd(Dir)
  
  #Carregar Correntes U & V
  u <- stack(list.files(pattern=".nc"),varname="u")
  v <- stack(list.files(pattern=".nc"),varname="v")
  
  #Rotacionar as correntes
  u <- raster::rotate(u)
  v <- raster::rotate(v)
  
  #Importar Mascara
  mask <- raster("Mask.asc")
  
  #Cortar para a area desejada
  uBRA <- crop(u,extent(mask)+5)
  vBRA <- crop(v,extent(mask)+5)
  
  #Tamanho de celula
  uBRA <- disaggregate(uBRA, fact=4)#Factor 4 transforma em celulas de 10km
  vBRA <- disaggregate(vBRA, fact=4)#Factor 4 transforma em celulas de 10km
  
  #Cortar para a area desejada
  uBRA <- crop(uBRA,extent(mask))
  vBRA <- crop(vBRA,extent(mask))
  
  #Fazer o Mar chegar na Praia
  u3 <- is.na(uBRA[[1]])==T & is.na(mask)==F
  u3[u3!=1] <- NA
  xyNA <- rasterToPoints(u3)
  u3E <-  extract(u3,xyNA[,1:2],cellnumbers=T)
  u4 <- uBRA[[1]]
  u4[boundaries(uBRA[[1]])==0] <- NA
  xy4 <- rasterToPoints(u4)
  u4E <-  extract(u4,xy4[,1:2],cellnumbers=T)
  
  v3 <- is.na(vBRA[[1]])==T & is.na(mask)==F
  v3[v3!=1] <- NA
  xyNA <- rasterToPoints(v3)
  v3E <-  extract(v3,xyNA[,1:2],cellnumbers=T)
  v4 <- vBRA[[1]]
  v4[boundaries(vBRA[[1]])==0] <- NA
  xy4 <- rasterToPoints(v4)
  v4E <-  extract(v4,xy4[,1:2],cellnumbers=T)
  
  #MetodoMSDM
  u3C<-as(u3,'SpatialPixels')@coords
  u4C<-as(u4,'SpatialPixels')@coords
  dist_u<-dist2(u3C,u4C,method='euclidean',p=2)
  dist_u[dist_u==0] <- NA
  distmin_u<-apply(dist_u,1,function(x) which(x==min(x, na.rm=T))) 
  distmin_u <- unlist(distmin_u)
  
  v3C<-as(v3,'SpatialPixels')@coords
  v4C<-as(v4,'SpatialPixels')@coords
  dist_v<-dist2(v3C,v4C,method='euclidean',p=2)
  dist_v[dist_v==0] <- NA
  distmin_v<-apply(dist_v,1,function(x) which(x==min(x, na.rm=T))) 
  distmin_v <- unlist(distmin_v)
  
  u5 <- uBRA
  for(x in 1:nlayers(u5)){
    u5[[x]][u3E[,1]] <- u4E[distmin_u,2]
  }  
  
  v5 <- uBRA
  for(x in 1:nlayers(v5)){
    v5[[x]][v3E[,1]] <- v4E[distmin_v,2]
  }  
  
  #Fazer o mar chegar na praia
  u2 <- uBRA
  for(x in 1:nlayers(u2)){
    u2[[x]][u3E[,1]] <- u5[[x]][u3E[,1]]
  } 

  v2 <- vBRA
  for(x in 1:nlayers(v2)){
    v2[[x]][v3E[,1]] <- v5[[x]][v3E[,1]]
  }

  #Salvar as camadas
  writeRaster(u2,names(u2),format="ascii",bylayer=T,NAflag=-9999)
  writeRaster(v2,names(v2),format="ascii",bylayer=T,NAflag=-9999)
}