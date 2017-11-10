setwd("C:\\OneDrive\\SDM_Aedes\\Variaveis_Matlab\\Pre")

library(SDMTools)
library(raster)

lst<-list.files(pattern=".asc")
asc<-read.asc(lst[1])
ras<-raster(asc)

rsd<-as.data.frame(ras,xy=T,centroids=T)
long<-rsd
lat<-rsd
lins<-which(is.na(long[,3])==F)

long[lins,3]<-long[lins,1]
gridded(long)<-~x+y
long<-raster(long)
writeRaster(long,"Longitude.asc",format="ascii")

lat[lins,3]<-lat[lins,2]
gridded(lat)<-~x+y
lat<-raster(lat)
writeRaster(lat,"Latitude.asc",format="ascii")