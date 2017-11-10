
## CAP II - AGRICULTURAL EXPANSION
## TRANSFORM A CSV FILE INTO RASTER

# http://gis.stackexchange.com/questions/20018/how-can-i-convert-data-in-the-form-of-lat-lon-value-into-a-raster-file-using-r



rm(list=ls())

## Escolher o diretório onde estão os arquivos que serão usados:
setwd("C:\\Users\\decoa\\Desktop\\Livia")


pts=read.table("Matriz_Presenca_Ausencia.txt", head=T, sep="\t")
colnames(pts)

soma<-rowSums(pts[,-c(1,2)])
pts<-cbind(pts,soma)

phyP<-pts[,c(1,2,239)] #Coluna 239 é a riqueza

gridded(phyP)<-~x+y
phyP.r<-raster(phyP)
plot(phyP.r)
#gridded(funcP)<- ~ x+y

#richP<-pts[,c(1,2,4)]

# Convert the data frame to a SpatialPointsDataFrame using the sp package and something like:
  
library(sp)
library(rgdal)
library(raster)

coordinates(phyP)=~x+y

# carregar cenários futuros de expansão

setwd("C:\\Users\\LIVIA LAURETO\\Desktop\\Análises de diversidade\\Expansão")

a1b<-raster("a1b_2100_livia_modelo_pronto_016.tif")
a2<-raster("a2_2100_livia_modelo_pronto_016.tif")
b1<-raster("b1_2100_livia_modelo_pronto_016.tif")

# cenário atual:

im16<-raster("image_2010_pronto_016.tif")

# chamar um raster de expansão e copiar informações de projeção (coor.ref.)

a1b
plot(a1b)


# Convert to your regular km system by first telling it what CRS it is, and then spTransform to the destination.

proj4string(phyP)=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")


# Colar inf de projeção após o CRS
phyP = spTransform(phyP,CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

# Tell R that this is gridded:
  
gridded(phyP) = TRUE

# At this point you'll get an error if your coordinates don't lie on a nice regular grid.

# Now use the raster package to convert to a raster and set its CRS:
  
r = raster(phyP)
projection(r) = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

# Now have a look:
  
plot(r)

# cortar o raster de diversidade para o msm extent do raster de expansão


extent???


# Now write it as a geoTIFF file using the raster package:
  
writeRaster(r,"MPDphyP.tif")



## overlay



