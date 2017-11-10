Binary.Maps("C:\\Users\\decoa\\Desktop\\OdonatasBrasil\\Result\\ENS","C:\\Users\\decoa\\Desktop\\OdonatasBrasil\\Result\\BIN")

Binary.Maps <- function(dir_input,dir_output){

  #Funcao para transformar mapas de adequabilidade em mapas binarios (Pres/AUS):
  #Parametros:
    #dir_input: diretorio com os mapas de adequabilidade (.asc)
    #dir_output: diretorio onde serao salvos os mapas binarios
  #Obs: A tabela com os thresholds deve conter as especies na 1ª linha e os thresholds na segunda
  
library(SDMTools)
library(raster)
  
setwd(dir_input)

print("Selecione a tabela com os thresholds: (1a coluna->ssp / 2a coluna->thresholds")    
thrs<-read.table(file.choose(),header=T,sep="\t")
zeros <- which(thrs[,2]==0)
thrs <- thrs[-zeros,]

modelos<-list.files(pattern = "\\.asc$")
asc.base<-read.asc(modelos[1])
ras<-raster(asc.base)
camadas<-modelos[-1]
    
for (b in 1:length(camadas)){
  asc<-read.asc(paste(dir_input,camadas[b],sep="\\"))
  ras.t<-raster(asc)
  ras<-stack(ras,ras.t)
}
names(ras)<-substr(modelos,1,nchar(modelos)-5)
ras <- ras[[-zeros]]

setwd(dir_output)
    
for (c in 1:nrow(thrs)){
  th.sp<-thrs[c,2]
  bin.sp<-ras[[c]]>=th.sp
  writeRaster(bin.sp,paste(names(ras[[c]]),"binary.asc",sep="_"),format="ascii")
}
}
