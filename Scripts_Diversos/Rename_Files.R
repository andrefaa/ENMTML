Rename_Files <-function(diretorio){
  #Funcao para renomear todos os arquivos em uma pasta
  
  setwd(diretorio)
  
  lista<-list.files(pattern=".mat")
  
  lista.t<-as.numeric(substr(lista,10,11))
  for (x in 1:length(lista)){
    if (lista.t[x]<10){
      lista[x]<-sub( '(?<=.{9})', '0', lista[x], perl=TRUE )  
    }  
  }
  file.rename(list.files(pattern=".mat"), lista)
}

Rename_Files("C:\\OneDrive\\IBM_Tubastraea\\IBM_Layers\\Month_SST\\Rio Grande do Sul(Guaiba)")
