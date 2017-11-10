Resample<-function (diretorio_input,diretorio_output){
  
  #Faz o resample no asc base(muda o tamanho da célula e a extensao com base em uma máscara)
  
  #Parametros
    #diretorio_output:diretorio para salvar os arquivos cortados
  
  library(raster)
  library(SDMTools)
  
  #Seleciona odiretorio com os ascs que serao alterados e máscara(que possui a resolucao desejada)
  
  setwd(diretorio_input)
  
  env_prep<-list.files(pattern='.asc')
  asc_base<-read.asc(env_prep[1])
  environmental<-raster(asc_base)
  camadas<-env_prep[-1]
  
  if (length(env_prep)!=1){
  
    for (a in 1:length(camadas)){
      asc<-read.asc(camadas[a])
      raster<-raster(asc)
      environmental<-stack(environmental,raster)
      print(paste("Lendo variaveis ambientais...",camadas[a],sep=""))
    }
    names(environmental)<-substr(env_prep,1,nchar(env_prep)-4)
  }
  
  rm(asc_base)
  rm(asc)
  rm(raster)
    
  print("Selecione o asc(ou raster) que possui a resoluçao desejada:")
  asc_mascara<-read.asc(file.choose())
  asc_mascara<-raster(asc_mascara)
  
  #Executa o resample
  
  asc_resam<-resample(environmental,asc_mascara)
  asc_resam@file@nodatavalue<--9999
  rm(environmental)
  rm(asc_mascara)
  
  #Salva o arquivo com a nova resolucao
  
  setwd(diretorio_output)
  
  for (a in 1:nlayers(asc_resam)){
    writeRaster(asc_resam[[a]],paste(names(asc_resam)[a],".asc",sep=""),format="ascii")
  }
  
}

Resample("C:\\Users\\decoa\\Desktop\\Variaveis_Daniel","C:\\Users\\decoa\\Desktop\\Variaveis_Daniel\\Res")
