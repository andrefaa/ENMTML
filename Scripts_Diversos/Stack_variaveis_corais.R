setwd ("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)")

files<-list.files("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)")

library(raster)

v01 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca1.asc")
v02 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca2.asc")
v03 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca3.asc")
v04 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca4.asc")
v05 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca5.asc")
v06 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca6.asc")
v07 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca7.asc")
v08 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca8.asc")
v09 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/pca9.asc")
v10 <- raster("C:/Users/Andre/Desktop/Andre/Mestrado/Layers/BioOracle/PCA(BioOracle)/bathy5km.asc")

stack <- stack(v01,v02,v03,v04,v05,v06,v07,v08,v09,v10)
e<- extent(-60,-20,-40,10)

costabrasil<- crop(stack,e)

#z <- raster(ncol=4320, nrow=1680)
#zcrop <- crop(z,e)
#res(zcrop) <- 0.02083333
#raster de resolucao 4km

#cor_resPRE <- resample(coraisPRE, zcrop, method='bilinear')
#cor_resFUT <- resample(coraisFUT, zcrop, method='bilinear')
#bathy_res <- resample(bathy, zcrop, method='bilinear')

writeRaster(costabrasil, "var.asc", format="ascii",bylayer=TRUE)


#for (x in files){
  r <- raster(files[x])
  f <- extract(r,e)
  writeRaster(f, ".asc", format="ascii",bylayer=TRUE)
}

