DummyVariables <- function(mask){

  #Funcao para transformar uma variavel categorica em variaveis dummy
  
  library(rgdal)
  library(raster)
  library(dummies)
  library(SDMTools)

  dir_in<-"C:\\Users\\decoa\\Desktop\\Variaveis_Laila"
  setwd(dir_in)

  lis<-list.files(pattern=".shp")
  shape<- substr(lis[2],1,nchar(lis[2])-4)
  shape <- readOGR(dsn = ".", layer = shape)
  msk.ras <- raster(read.asc(file.choose()))
  
  codes<-cbind(data.frame(shape$CLASS_NAME),data.frame(shape$GRIDCODE))
  codes <- codes[with(codes, order(shape.GRIDCODE)), ]
  clas <- as.character(unique(codes[,1]))[-10]
  
  
  ras <- rasterize(shape,msk.ras,field=shape$GRIDCODE)
  plot(ras)
  
  
  ras.frame <- rasterToPoints(ras)
  dum<-dummy(ras.frame[,3])
  colnames(dum) <- clas
  
  dum <- cbind(ras.frame[,1:2],dum)
  names(dum) <- "x"
  
  for (x in 1:length(clas)){
    temp <- as.data.frame(dum[,c(1,2,(x+2))])
    names(temp[,1]) <- "x"
    gridded (temp) <- ~ x+y
    temp <- raster(temp)
    writeRaster(temp,paste(names(temp),".asc",sep=""),format="ascii")
  }
}

