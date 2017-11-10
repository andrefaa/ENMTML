Binary.Maps("C:\\OneDrive\\Poliquetas\\Poliquetas_Pacifico\\Resultados\\Generos_AMS","C:\\OneDrive\\Poliquetas\\Poliquetas_Pacifico\\Resultados\\Generos_AMS")

Binary.Maps <- function(dir_input,dir_output){

  #Funcao para transformar mapas de adequabilidade em mapas binarios (Pres/AUS):
  #Parametros:
    #dir_input: diretorio com os mapas de adequabilidade (.asc)
    #dir_output: diretorio onde serao salvos os mapas binarios
  #Obs: A tabela com os thresholds deve conter as especies na 1ª linha e os thresholds na segunda
  
library(SDMTools)
library(raster)

print("Selecione a tabela com os thresholds: (1ª coluna->ssp / Outras coluna->Algoritmos e seus thresholds(ordem alfabetica)")    
thrs<-read.table(file.choose(),header=T,sep="\t")

algor<-list.files(dir_input)

  for (a in 1:length(algor)){
    modelos<-list.files(paste(dir_input,algor[a],sep="\\"),pattern = ".asc")
    asc.base<-read.asc(paste(dir_input,algor[a],modelos[1],sep="\\"))
    ras<-raster(asc.base)
    camadas<-modelos[-1]
    
    for (b in 1:length(camadas)){
      asc<-read.asc(paste(dir_input,algor[a],camadas[b],sep="\\"))
      ras.t<-raster(asc)
      ras<-stack(ras,ras.t)
    }
    names(ras)<-substr(modelos,1,nchar(modelos)-5)
    setwd(dir_output)
    
    for (c in 1:nrow(thrs)){
      th.sp<-thrs[c,a+1]
      bin.sp<-ras[[c]]>=th.sp
      writeRaster(bin.sp,paste(algor[a],names(ras[[c]]),"binary.asc",sep="_"),format="ascii")
    }
  }
}
