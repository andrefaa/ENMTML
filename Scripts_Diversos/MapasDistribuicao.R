library(xlsx)
library(raster)
library(ggmap)
library(maptools)

occ<-read.table(file.choose(),sep="\t",header = T)
occ[,2]<-as.numeric(as.character(occ[,2]))
occ[,3]<-as.numeric(as.character(occ[,3]))
ssp<-unique(occ[,1])
data(wrld_simpl)
bra<-crop(wrld_simpl,extent(-70,-20,-35,10))

for (x in ssp){
  occ.ssp<-occ[occ[,1]==x,]
  tiff(file = paste(as.character(x),".tiff",sep=""), width = 800, height = 800, units = "px", res = 200)
  plot(bra)
  points(occ.ssp[,2:3],cex=0.75,pch=16,col="red")
  title(as.character(x))
  dev.off()
}