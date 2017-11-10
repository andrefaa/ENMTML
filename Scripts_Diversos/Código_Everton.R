Codigo_Everton<-function(diretorio){

setwd(diretorio)

library(raster)
library(SDMTools)

env_prep<-list.files(pattern='.asc')
asc_base<-read.asc(env_prep[1])
environmental<-raster(asc_base)
camadas<-env_prep[-1]

for (a in 1:length(camadas)){
  asc<-read.asc(camadas[a])
  raster<-raster(asc)
  environmental<-stack(environmental,raster)
}
names(environmental)<-substr(env_prep,1,nchar(env_prep)-4)

tabela<-rasterToPoints(environmental)

write.table(tabela,"ASCs_Tabela.txt",row.names=F,sep="\t")
}

Codigo_Everton("C:\\Users\\Andre\\Desktop\\Everton")
