ENM_Stack<-function(diretorio){

setwd(diretorio)
library(raster)
library(SDMTools)

#Cria um bloco com os ENMs das especies
env_prep<-list.files(pattern='.asc')
asc_base<-read.asc(env_prep[1])
bloco<-raster(asc_base)
camadas<-env_prep[-1]

for (a in 1:length(camadas)){
  asc<-read.asc(camadas[a])
  raster<-raster(asc)
  bloco<-stack(bloco,raster)
}
names(bloco)<-substr(env_prep,1,nchar(env_prep)-4)

#Importar a tabela de avaliacao das especies(Feita no Matlab, com o codigo roc_asc5)
avaliacoes<-read.table(file.choose(),h=T,sep='\t')
th_roc<-avaliacoes$th_roc

#Criar os mapas binários

binarios<-bloco

for (b in 1:length(th_roc)){
  binarios[[b]]<-binarios[[b]]>=th_roc[b]
}
names(binarios)<-substr(env_prep,1,nchar(env_prep)-4)

#Stack nos rasters para criar mapas de riqueza
riqueza<-raster(ncol=ncol(binarios),nrow=nrow(binarios),ext=extent(binarios),vals=0)

for (c in 1:nlayers(binarios)){
  riqueza<-riqueza+binarios[[c]]
}

writeRaster(riqueza,"Riqueza_SDM.asc",format='ascii')

}

ENM_Stack('C:\\Users\\Andre\\Google Drive\\Mestrado\\Poliquetas\\ENM\\Modelos\\25km\\Resultados\\Subset-BioOracle\\MAXENT')

###### Z-transform raster ####
#raster<-(raster−minValue(raster))/(maxValue(raster)−minValue(raster))

