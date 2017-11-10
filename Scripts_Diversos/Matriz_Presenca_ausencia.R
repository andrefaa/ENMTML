presence.absence.raster <- function (mask.raster,species.data) {
  
  # set the background cells in the raster to 0
  mask.raster[!is.na(mask.raster)] <- 0
  
  #set the cells that contain points to 1
  speciesRaster <- rasterize(species.data,mask.raster,field=1)
  speciesRaster <- merge(speciesRaster,mask.raster)
  presenca.ausencia.matriz.ssp <- rasterToPoints(speciesRaster)
  presenca.ausencia.matriz<- cbind(presenca.ausencia.matriz,presenca.ausencia.matriz.ssp[,3])
  
  
  #label the raster
  #names(speciesRaster) <- raster.label
  #return(speciesRaster)
}  

#ALETRAR SOMENTE OS NOMES DOS ARQUIVOS(.ASC E .CSV) DESEJADOS ABAIXO!!!!!!!
setwd("C:/Users/Andre/Google Drive/Mestrado/Poliquetas")
ssp<-list.files(pattern='.csv')

for (i in 1:length(ssp)){
  
  #Import the ascii file that will serve as mask
  asc<-read.asc("sstmean_05.asc")
  myRaster<- raster(asc)
  myRaster<-crop(myRaster,extent)
  plot(myRaster)
  
  #Import the occurrence points (csv format; 2 columns: [1]=longitude, [2]=latitude)
  sp_atual<- read.csv(ssp)
  lista_ssp<-unique(sp_atual[,1])
  presenca.ausencia.matriz<-NULL
  
  for (x in 1:length(lista_ssp)){
    ocorrencias.especie<-sp_atual[sp_atual[,1]==lista_ssp[1],]
    myspecies <- ocorrencias.especie[,c(2:3)]
    
  #Run the function to create the binary map
  presenca.ausencia.matriz<-presence.absence.raster(myRaster,myspecies)
  
  #Save ascii output
  #writeRaster(species.0.1,strcat(ssp[i],"_binary"),format="ascii")
  
}
colnames(presenca.ausencia.matriz)<-lista_ssp
presenca.ausencia.matriz.georef<-cbind(presenca.ausencia.matriz.ssp[,c(1,2)],presenca.ausencia.matriz)
riqueza.celula<-NULL

for (a in 1:nrow(presenca.ausencia.matriz.georef)){
  riqueza.celula[a]<-sum(presenca.ausencia.matriz.georef[a,-c(1,2)])
}


write.table(presenca.ausencia.matriz.georef,"Matriz_Presenca_ausencia_Poliquetas.txt",sep="\t",row.names=F)
