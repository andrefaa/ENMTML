library(raster)

setwd('C:\\MODEL_FUTURO\\esemble_futuro\\SVM\\adeq')
ascc<-read.asc('C:\\MODEL_FUTURO\\esemble_futuro\\SVM\\Thylamys_velutinus.asc')
ras<-raster(ascc)
ras2<-(ras-minValue(ras))/(maxValue(ras)-minValue(ras))
writeRaster(ras2,paste("Thylamys_velutinus",".asc", sep=""),overwrite=TRUE)

setwd('C:\\MODEL_FUTURO\\esemble_futuro\\MXS\\adeq')
ascc<-read.asc('C:\\MODEL_FUTURO\\esemble_futuro\\MXS\\Thylamys_velutinus.asc')
ras<-raster(ascc)
ras2<-(ras-minValue(ras))/(maxValue(ras)-minValue(ras))
writeRaster(ras2,paste("Thylamys_velutinus",".asc", sep=""),overwrite=TRUE)

plot(ras2)

asc<-read.asc('C:\\MODEL_FUTURO\\esemble_futuro\\MXS\\adeq\\Glyphonycteris_behnii.asc')
ras<-raster(asc)
plot(ras)
