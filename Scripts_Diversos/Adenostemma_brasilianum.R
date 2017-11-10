library(rgbif)
library(maps)

e<-extent(-80,-30,-35,0)

adens<-gbif("Adenostemma","brasilianum",download=T,concept=T,ext=e)
gbif("Bidens","brasiliensis", download=F,concept=T)
names(adens)

adens<-data.frame(adens$species,adens$lat,adens$lon)
colnames(adens)<- c("species","lat","lon")

A_brasilianum_completo <- rbind(adens,Adenostemma_brasilianum)

write.table(A_brasilianum_completo, "C:/Users/Andre/Desktop/Andre/Mestrado/Ocorrencias/Ocorrencias CSV/Adenostemma_brasilianum_completo.csv", row.names=F)

#map(xlim=c(-80,-30),ylim=c(-35,0))

